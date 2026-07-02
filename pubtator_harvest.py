#!/usr/bin/env python3
"""
PubTator3 → TTN 致病變異彙整器（獨立測試版）

與主 pipeline（variant-centric：給變異 → 找文獻）相反，本程式是 literature-centric：
    給 PubTator 上 TTN 基因的 journal article 文獻 → 抽出文中提及的 TTN 變異
    → 用本地 LLM 判定致病性並擷取 disease / 影響組織 / 發病年齡 / 遺傳模式
    → 把所有文獻的結果彙整成單一 CSV。

設計重點（皆已實測 PubTator3 API 行為）：
  1. 文獻清單用「基因實體」查詢 @GENE_TTN，而非 free-text `ttn`
     （free-text 會混入 Transient Tachypnea of the Newborn、tensor-network TTNS 等雜訊）。
  2. journal article 過濾：PubTator search 的 filters 參數實測無效，改用 NCBI esummary 的
     pubtype 來判定（順便取得年份）。
  3. 變異來源：PubTator 逐篇 full=true biocjson export（批次 full=true 會回傳錯誤文章，
     故逐篇抓），取 Variant 標註中 CorrespondingGene:7273（TTN）者，含 HGVS / rsID。
  4. 表型欄位 PubTator 沒有，沿用既有 LocalClinicalExtractor 的 schema 由本地 LLM 擷取，
     額外加上 pathogenicity 欄位（由 LLM 依文中語境判定是否致病/disease-causing）。

用法：
    # 只跑 Stage 1-2（不需 GPU），先驗證變異抓取
    conda run -n gm4 python pubtator_harvest.py --max-articles 5 --no-llm

    # 完整流程（需 vLLM / GPU）
    conda run -n gm4 python pubtator_harvest.py --max-articles 100
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import requests

from config import (
    OUTPUT_DIR,
    PUBMED_API_KEY,
    PUBMED_EMAIL,
    LOCAL_LLM_BACKEND,
    LOCAL_LLM_MODEL,
    LOCAL_LLM_TENSOR_PARALLEL,
    LOCAL_LLM_MAX_MODEL_LEN,
    LOCAL_LLM_MAX_CONTEXT_LENGTH,
    TRITON_PTXAS_PATH,
)

# 與 main.py 一致：在載入 vLLM/Triton 前設定 PTXAS 路徑（vLLM 編譯需要）
os.environ.setdefault("TRITON_PTXAS_PATH", TRITON_PTXAS_PATH)

try:
    from config import LLM_SELF_CONSISTENCY_N, LLM_SAMPLING_TEMPERATURE, LLM_USE_GUIDED_JSON
except ImportError:
    LLM_SELF_CONSISTENCY_N = 1
    LLM_SAMPLING_TEMPERATURE = 0.2
    LLM_USE_GUIDED_JSON = False

# 同時在 vLLM 內跑的序列上限（控制 KV cache 峰值，避免長 context OOM）
try:
    from config import LLM_MAX_NUM_SEQS
except ImportError:
    LLM_MAX_NUM_SEQS = 16

# vLLM 啟動要求的 GPU 記憶體比例。unified memory（DGX Spark）下，啟動時可用量會
# 受系統與本行程（一次載入全部文章）影響，0.90 太貪會起不來，預設留多一點。
try:
    from config import LLM_GPU_MEMORY_UTILIZATION
except ImportError:
    LLM_GPU_MEMORY_UTILIZATION = 0.82

# 同時輸出到 stdout 與固定 log 檔（append，不用每次手動 tee）
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE = OUTPUT_DIR / "pubtator_harvest.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_FILE, mode="a", encoding="utf-8"),
    ],
)
logger = logging.getLogger("pubtator_harvest")
logger.info("=" * 70)
logger.info(f"=== 新一次執行開始，log 寫入 {LOG_FILE} ===")

TTN_GENE_ID = "7273"  # human TTN NCBI Gene ID
PUBTATOR_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CACHE_DIR = OUTPUT_DIR / "pubtator_cache"
ABSTRACT_DIR = CACHE_DIR / "abstract"      # 可靠的 abstract 層級 biocjson（id 一定正確）
FULLTEXT_DIR = CACHE_DIR / "fulltext"      # full=true 全文，僅在 id 與請求一致時採用
EUPMC_DIR = CACHE_DIR / "eupmc"            # EuropePMC 全文（純文字）
DEFAULT_OUT = OUTPUT_DIR / "pubtator_ttn_variants.csv"


# ════════════════════════════════════════════════════════════════════════════
# Stage 1 + 2 : PubTator3 / NCBI client
# ════════════════════════════════════════════════════════════════════════════
class PubTatorClient:
    def __init__(self, sleep: float = 0.34):
        self.sleep = sleep  # 禮貌性間隔（無 API key 約 3 req/s）
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "ttn-agent-pubtator-harvest/1.0"})
        for d in (ABSTRACT_DIR, FULLTEXT_DIR, EUPMC_DIR):
            d.mkdir(parents=True, exist_ok=True)

    # ── 通用 GET（含重試） ──────────────────────────────────────────────────
    def _get(self, url: str, params: Optional[dict] = None, timeout: int = 60,
             retries: int = 4) -> Optional[requests.Response]:
        for attempt in range(1, retries + 1):
            try:
                r = self.session.get(url, params=params, timeout=timeout)
                if r.status_code == 200:
                    return r
                logger.warning(f"  HTTP {r.status_code} ({attempt}/{retries}) {url}")
            except requests.RequestException as e:
                logger.warning(f"  request error ({attempt}/{retries}): {e}")
            time.sleep(self.sleep * attempt * 2)
        return None

    # ── Stage 1: 搜尋 @GENE_TTN，取 relevance 排序的 PMID ──────────────────
    def search_ttn_pmids(self, max_candidates: int) -> List[Dict[str, Any]]:
        results: List[Dict[str, Any]] = []
        page = 1
        total_pages = None
        while len(results) < max_candidates:
            r = self._get(
                f"{PUBTATOR_BASE}/search/",
                params={"text": "@GENE_TTN", "sort": "score desc", "page": page},
            )
            if r is None:
                break
            try:
                data = r.json()
            except json.JSONDecodeError:
                logger.error("  search 回傳非 JSON，停止分頁")
                break
            batch = data.get("results", [])
            if not batch:
                break
            if total_pages is None:
                total_pages = data.get("total_pages")
                logger.info(
                    f"@GENE_TTN 命中 {data.get('count')} 篇（{total_pages} 頁），"
                    f"開始收集前 {max_candidates} 個候選 PMID"
                )
            for item in batch:
                results.append({
                    "pmid": str(item["pmid"]),
                    "pmcid": item.get("pmcid", ""),
                    "title": item.get("title", ""),
                    "journal": item.get("journal", ""),
                    "score": item.get("score"),
                })
            logger.info(f"  page {page}: 累積 {len(results)} 候選")
            if total_pages and page >= total_pages:
                break
            page += 1
            time.sleep(self.sleep)
        return results[:max_candidates]

    # ── Stage 1b: esummary 取 pubtype / year，過濾 journal article ─────────
    def annotate_pubtypes(self, pmids: List[str]) -> Dict[str, Dict[str, Any]]:
        out: Dict[str, Dict[str, Any]] = {}
        for i in range(0, len(pmids), 200):
            chunk = pmids[i:i + 200]
            params = {
                "db": "pubmed",
                "id": ",".join(chunk),
                "retmode": "json",
                "email": PUBMED_EMAIL,
            }
            if PUBMED_API_KEY:
                params["api_key"] = PUBMED_API_KEY
            r = self._get(f"{EUTILS_BASE}/esummary.fcgi", params=params)
            if r is None:
                continue
            try:
                res = r.json().get("result", {})
            except json.JSONDecodeError:
                continue
            for pmid in res.get("uids", []):
                rec = res.get(pmid, {})
                pubtypes = rec.get("pubtype", []) or []
                year = ""
                for key in ("sortpubdate", "pubdate", "epubdate"):
                    if rec.get(key):
                        m = re.search(r"\d{4}", rec[key])
                        if m:
                            year = m.group(0)
                            break
                out[pmid] = {"pubtypes": pubtypes, "year": year}
            time.sleep(self.sleep)
        return out

    # ── Stage 2a: 批次 abstract 層級 export（可靠，id 一定正確） ────────────
    def fetch_abstract_docs(self, pmids: List[str], force: bool = False) -> Dict[str, dict]:
        out: Dict[str, dict] = {}
        todo: List[str] = []
        for pmid in pmids:
            cache = ABSTRACT_DIR / f"{pmid}.json"
            if cache.exists() and not force:
                try:
                    out[pmid] = json.loads(cache.read_text(encoding="utf-8"))
                    continue
                except json.JSONDecodeError:
                    pass
            todo.append(pmid)
        for i in range(0, len(todo), 100):
            chunk = todo[i:i + 100]
            r = self._get(
                f"{PUBTATOR_BASE}/publications/export/biocjson",
                params={"pmids": ",".join(chunk)}, timeout=90,
            )
            if r is None:
                continue
            try:
                docs = r.json().get("PubTator3", [])
            except json.JSONDecodeError:
                logger.warning("  abstract export 解析失敗")
                docs = []
            for doc in docs:
                did = str(doc.get("id", ""))
                if did:
                    out[did] = doc
                    (ABSTRACT_DIR / f"{did}.json").write_text(
                        json.dumps(doc, ensure_ascii=False), encoding="utf-8")
            time.sleep(self.sleep)
        return out

    # ── Stage 2b: 批次 full=true 全文 export ──────────────────────────────────
    # 注意：full=true 回傳的 doc['id'] 是 PMCID 號、且批次順序會被打亂，
    # 必須用 passage 的 article-id_pmid 對應回原 PMID（見 doc_pmid()）。
    def fetch_fulltext_docs(self, pmids: List[str], force: bool = False) -> Dict[str, dict]:
        out: Dict[str, dict] = {}
        todo: List[str] = []
        for pmid in pmids:
            cache = FULLTEXT_DIR / f"{pmid}.json"
            if cache.exists() and not force:
                try:
                    out[pmid] = json.loads(cache.read_text(encoding="utf-8"))
                    continue
                except json.JSONDecodeError:
                    pass
            todo.append(pmid)
        for i in range(0, len(todo), 50):
            chunk = todo[i:i + 50]
            r = self._get(
                f"{PUBTATOR_BASE}/publications/export/biocjson",
                params={"pmids": ",".join(chunk), "full": "true"}, timeout=180,
            )
            time.sleep(self.sleep)
            if r is None:
                continue
            try:
                docs = r.json().get("PubTator3", [])
            except json.JSONDecodeError:
                logger.warning("  full=true export 解析失敗")
                continue
            for doc in docs:
                pmid = doc_pmid(doc)
                if not pmid:
                    continue
                out[pmid] = doc
                (FULLTEXT_DIR / f"{pmid}.json").write_text(
                    json.dumps(doc, ensure_ascii=False), encoding="utf-8")
        return out

    # ── Stage 2c: EuropePMC 全文（可靠，僅供 LLM context，不作為變異來源） ──
    def fetch_eupmc_fulltext(self, pmid: str, max_len: int = 200000,
                             force: bool = False) -> Optional[str]:
        cache = EUPMC_DIR / f"{pmid}.txt"
        if cache.exists() and not force:
            txt = cache.read_text(encoding="utf-8")
            return txt or None
        # 直接請求；404 代表無開放全文（正常），不噪音化成 warning
        r = None
        try:
            resp = self.session.get(
                f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmid}/fullTextXML",
                timeout=60)
            if resp.status_code == 200:
                r = resp
        except requests.RequestException:
            pass
        time.sleep(self.sleep)
        if r is None or not r.content:
            cache.write_text("", encoding="utf-8")
            return None
        try:
            from bs4 import BeautifulSoup
            soup = BeautifulSoup(r.content, "xml")
            body = soup.find("body")
            text = body.get_text(" ") if body else ""
            text = re.sub(r"\s+", " ", text).strip()[:max_len]
        except Exception as e:
            logger.debug(f"  PMID {pmid}: EuropePMC 解析失敗 {e}")
            text = ""
        cache.write_text(text, encoding="utf-8")
        return text or None


# ════════════════════════════════════════════════════════════════════════════
# biocjson 解析：抽 TTN 變異 + 組裝全文
# ════════════════════════════════════════════════════════════════════════════
def _parse_variant_infon(identifier: str) -> Dict[str, str]:
    """把 PubTator Variant 的 identifier 字串拆成欄位。
    例： 'tmVar:p|FS|P|12815|T|37;HGVS:p.P12815TfsX37;VariantGroup:0;
          CorrespondingGene:7273;RS#:752101551;CorrespondingSpecies:9606'
    """
    fields: Dict[str, str] = {}
    for part in (identifier or "").split(";"):
        if ":" in part:
            k, _, v = part.partition(":")
            fields[k.strip()] = v.strip()
        elif part.startswith("tmVar"):
            fields["tmVar"] = part
    return fields


def _is_ttn_variant(fields: Dict[str, str]) -> bool:
    for key in ("CorrespondingGene", "OriginalGene", "Gene"):
        val = fields.get(key, "")
        if TTN_GENE_ID in re.split(r"[ ,]+", val):
            return True
    return False


_VARIANT_TYPES = ("variant", "mutation", "dnamutation", "proteinmutation", "snp")


def extract_variants_from_doc(doc: Optional[dict],
                              variants: Optional[Dict[str, Dict[str, Any]]] = None
                              ) -> Dict[str, Dict[str, Any]]:
    """從一份 biocjson doc 抽出 TTN 變異，併入（或建立）variants dict。"""
    if variants is None:
        variants = {}
    if not doc:
        return variants
    for p in doc.get("passages", []):
        for a in p.get("annotations", []):
            ainf = a.get("infons", {})
            if (ainf.get("type") or "").lower() not in _VARIANT_TYPES:
                continue
            fields = _parse_variant_infon(ainf.get("identifier", ""))
            if not _is_ttn_variant(fields):
                continue
            hgvs = fields.get("HGVS", "")
            rsid = fields.get("RS#", "")
            if rsid and not rsid.startswith("rs"):
                rsid = "rs" + rsid
            raw = (a.get("text") or "").strip()
            key = hgvs or rsid or raw.lower()  # 變異身分：HGVS > rsID > 文字
            if not key:
                continue
            v = variants.setdefault(key, {
                "key": key, "hgvs": hgvs, "rsid": rsid,
                "texts": set(), "mention_count": 0,
            })
            if isinstance(v["texts"], list):  # 從快取讀回時可能已是 list
                v["texts"] = set(v["texts"])
            if raw:
                v["texts"].add(raw)
            if hgvs and not v["hgvs"]:
                v["hgvs"] = hgvs
            if rsid and not v["rsid"]:
                v["rsid"] = rsid
            v["mention_count"] += 1
    return variants


def doc_pmid(doc: Optional[dict]) -> Optional[str]:
    """從 biocjson doc 取真正的 PMID。

    full=true 端點回傳的 doc['id'] 其實是 PMCID 號（且批次會重排），
    真正的 PMID 在 passage 的 infons['article-id_pmid']。abstract 端點則 id 即 PMID。
    """
    if not doc:
        return None
    for p in doc.get("passages", []):
        pmid = p.get("infons", {}).get("article-id_pmid")
        if pmid:
            return str(pmid)
    did = str(doc.get("id", ""))
    return did or None


def extract_text_from_doc(doc: Optional[dict]) -> Dict[str, str]:
    """回傳 {pmcid, title, abstract, full_text}。"""
    pmcid, title_text, abstract_text = "", "", ""
    parts: List[str] = []
    if not doc:
        return {"pmcid": "", "title": "", "abstract": "", "full_text": ""}
    for p in doc.get("passages", []):
        infons = p.get("infons", {})
        ptype = (infons.get("type") or "").lower()
        section = (infons.get("section_type") or "").lower()
        ptext = p.get("text", "") or ""
        if infons.get("article-id_pmc") and not pmcid:
            pmcid = "PMC" + str(infons["article-id_pmc"]).replace("PMC", "")
        # title：abstract 端點 type=='title'；full 端點 section_type=='title'(type=='front')
        if (ptype == "title" or section == "title") and not title_text:
            title_text = ptext
        if (ptype == "abstract" or section == "abstract") and not abstract_text:
            abstract_text = ptext
        if ptext:
            parts.append(ptext)
    return {"pmcid": pmcid, "title": title_text,
            "abstract": abstract_text, "full_text": "\n".join(parts)}


# ════════════════════════════════════════════════════════════════════════════
# Stage 3 + 4 : 本地 LLM —— 致病性 + 表型擷取
# ════════════════════════════════════════════════════════════════════════════
PATHOGENICITY_ENUM = [
    "Pathogenic", "Likely pathogenic", "Uncertain significance",
    "Likely benign", "Benign", "Not specified",
]
TISSUE_ENUM = ["Cardiac", "Skeletal", "Both", "Not specified"]
AGE_ENUM = ["Congenital", "Pediatric", "Adult", "Not specified"]
INHERITANCE_ENUM = [
    "Autosomal Dominant", "Autosomal Recessive", "X-linked",
    "De Novo", "Sporadic", "Not specified",
]

EMPTY_RESULT = {
    "pathogenicity": "Not specified",
    "disease": "Not specified",
    "tissue_affected": "Not specified",
    "age_onset": "Not specified",
    "inheritance": "Not specified",
    "patient_count": 0,
    "evidence_sentence": "N/A",
    "reasoning": "",
}


class PathogenicVariantExtractor:
    """獨立的 vLLM 擷取器（沿用 config 的後端/模型設定）。

    與主 pipeline 的 LocalClinicalExtractor 不同：這裡是「給一篇文獻 + 一個變異」，
    判定該變異在文中是否被描述為致病，並抽出 disease / 組織 / 年齡 / 遺傳模式。
    """

    def __init__(self, model_name: Optional[str] = None, tokenizer: Optional[str] = None):
        self.model_name = model_name or LOCAL_LLM_MODEL
        self.tokenizer = tokenizer
        self.is_gemma = "gemma" in self.model_name.lower() or "medgemma" in self.model_name.lower()
        self.max_context_length = LOCAL_LLM_MAX_CONTEXT_LENGTH
        self._init_vllm()

    @staticmethod
    def _resolve_gguf(model_name: str) -> Tuple[str, bool]:
        """支援三種 GGUF 指定法，回傳 (vLLM 可用的 model 路徑, is_gguf)：
          - 'repo_id:file.gguf'  → 從 HF 下載該單檔
          - '/abs/path/x.gguf'   → 本地檔
          - 一般 HF repo / 路徑   → 原樣回傳（非 GGUF）
        """
        if ":" in model_name and model_name.lower().endswith(".gguf"):
            repo, _, fname = model_name.partition(":")
            from huggingface_hub import hf_hub_download
            logger.info(f"下載 GGUF：repo={repo} file={fname}")
            return hf_hub_download(repo_id=repo, filename=fname), True
        if model_name.lower().endswith(".gguf"):
            return model_name, True
        return model_name, False

    def _init_vllm(self):
        from vllm import LLM, SamplingParams
        stop_tokens = ["<end_of_turn>", "\n\n\n"] if self.is_gemma else ["\n\n\n", "```\n\n"]
        model_path, is_gguf = self._resolve_gguf(self.model_name)
        llm_kwargs: Dict[str, Any] = dict(
            model=model_path,
            tensor_parallel_size=LOCAL_LLM_TENSOR_PARALLEL,
            trust_remote_code=True,
            gpu_memory_utilization=LLM_GPU_MEMORY_UTILIZATION,
            max_model_len=LOCAL_LLM_MAX_MODEL_LEN,
            enforce_eager=True,
            enable_prefix_caching=True,
            # 限制同時在跑的序列數，避免長 context × n_samples 把 KV cache 塞爆
            # （上次 503 prompt 全丟進去在 54% 時 OOM: "Cannot get free blocks"）
            max_num_seqs=LLM_MAX_NUM_SEQS,
        )
        if is_gguf:
            # GGUF 內嵌 tokenizer 在 vLLM 支援有限，預設改用原始 HF tokenizer
            tok = self.tokenizer or ("google/gemma-4-31B-it" if self.is_gemma else None)
            if tok:
                llm_kwargs["tokenizer"] = tok
            logger.info(f"GGUF 模式：model={model_path} tokenizer={tok}")
        self.llm = LLM(**llm_kwargs)
        self.n_samples = max(1, int(LLM_SELF_CONSISTENCY_N))
        temperature = float(LLM_SAMPLING_TEMPERATURE) if self.n_samples > 1 else 0.2
        kwargs: Dict[str, Any] = dict(
            n=self.n_samples, temperature=temperature, top_p=0.95, top_k=64,
            max_tokens=2048, stop=stop_tokens,
        )
        if LLM_USE_GUIDED_JSON:
            try:
                from vllm.sampling_params import StructuredOutputsParams
                kwargs["structured_outputs"] = StructuredOutputsParams(json=self._schema())
                logger.info("Guided JSON decoding 已啟用")
            except Exception as e:
                logger.warning(f"無法啟用 guided JSON：{e}")
        self.sampling_params = SamplingParams(**kwargs)
        logger.info(f"vLLM 就緒：model={self.model_name} n={self.n_samples} temp={temperature}")

    @staticmethod
    def _schema() -> dict:
        return {
            "type": "object",
            "properties": {
                "reasoning": {"type": "string"},
                "pathogenicity": {"type": "string", "enum": PATHOGENICITY_ENUM},
                "disease": {"type": "string"},
                "tissue_affected": {"type": "string", "enum": TISSUE_ENUM},
                "age_onset": {"type": "string", "enum": AGE_ENUM},
                "inheritance": {"type": "string", "enum": INHERITANCE_ENUM},
                "patient_count": {"type": "integer", "minimum": 0},
                "evidence_sentence": {"type": "string"},
            },
            "required": ["reasoning", "pathogenicity", "disease", "tissue_affected",
                         "age_onset", "inheritance", "patient_count", "evidence_sentence"],
            "additionalProperties": False,
        }

    # ── context：以變異出現位置為中心開窗 ────────────────────────────────────
    def _build_context(self, article: Dict[str, Any], variant: Dict[str, Any]) -> str:
        ABSTRACT_CAP = 3500
        WINDOW = 4000
        NO_HIT_FALLBACK = 20000

        parts: List[str] = []
        budget = self.max_context_length
        if article.get("title"):
            parts.append(f"[TITLE]\n{article['title']}")
            budget -= len(article["title"]) + 10
        if article.get("abstract"):
            snip = article["abstract"][:ABSTRACT_CAP]
            parts.append(f"[ABSTRACT]\n{snip}")
            budget -= len(snip) + 12

        text = article.get("full_text", "") or ""
        aliases = [a for a in ([variant.get("hgvs"), variant.get("rsid")]
                               + list(variant.get("texts", []))) if a]
        if budget <= 0 or not text:
            return "\n\n".join(parts)

        low = text.lower()
        hits: List[int] = []
        for al in aliases:
            al = al.lower().strip()
            if len(al) < 3:
                continue
            idx = 0
            while True:
                pos = low.find(al, idx)
                if pos == -1:
                    break
                hits.append(pos)
                idx = pos + 1
        hits = sorted(set(hits))

        if not hits:
            parts.append(f"[FULL TEXT START]\n{text[:min(budget, NO_HIT_FALLBACK)]}")
            return "\n\n".join(parts)

        windows: List[Tuple[int, int]] = []
        for pos in hits:
            s, e = max(0, pos - WINDOW), min(len(text), pos + WINDOW)
            if windows and s <= windows[-1][1]:
                windows[-1] = (windows[-1][0], max(windows[-1][1], e))
            else:
                windows.append((s, e))
        chunks = []
        for s, e in windows:
            chunk = text[s:e].strip()
            if chunk:
                chunks.append(("..." if s > 0 else "") + chunk + ("..." if e < len(text) else ""))
        combined = "\n\n---\n\n".join(chunks)[:budget]
        parts.append(f"[RELEVANT SECTIONS (around variant mentions)]\n{combined}")
        return "\n\n".join(parts)

    def _build_prompt(self, context: str, variant: Dict[str, Any]) -> str:
        disp = variant.get("hgvs") or variant.get("rsid") or ", ".join(variant.get("texts", []))
        alias_lines = "\n".join(
            f"  - {a}" for a in dict.fromkeys(
                [variant.get("hgvs"), variant.get("rsid")] + list(variant.get("texts", [])))
            if a
        )
        system_message = f"""You are a clinical genetics expert extracting structured data from a biomedical article about the TTN (titin) gene.

## TARGET TTN VARIANT
Variant: {disp}
Known forms (any of these refer to the SAME variant):
{alias_lines}

## TASK
Read the article text and report, FOR THIS SPECIFIC TTN VARIANT ONLY:
1. Whether the article describes it as disease-causing / pathogenic.
2. The disease/phenotype it is associated with.
3. The affected tissue, age of onset, and inheritance pattern.

### RULES
1. Judge pathogenicity from the article's own language about THIS variant
   (e.g. "pathogenic", "disease-causing", "causative", "loss-of-function",
   segregates with disease, de novo in an affected patient → Pathogenic/Likely pathogenic;
   "benign", "polymorphism", "not associated" → Benign/Likely benign;
   "VUS", "uncertain", control-only → Uncertain significance).
   If the article gives no pathogenicity signal for this variant → "Not specified".
2. The variant may appear only in a table row — that row's cohort disease counts.
3. Use FULL disease names with abbreviation, e.g. "Dilated Cardiomyopathy (DCM)",
   "Tibial Muscular Dystrophy (TMD)", "Limb-Girdle Muscular Dystrophy (LGMD)",
   "Hereditary Myopathy with Early Respiratory Failure (HMERF)".
4. tissue_affected ∈ {{Cardiac, Skeletal, Both, Not specified}}.
5. age_onset ∈ {{Congenital, Pediatric, Adult, Not specified}}.
6. inheritance ∈ {{Autosomal Dominant, Autosomal Recessive, X-linked, De Novo, Sporadic, Not specified}}.
7. evidence_sentence MUST be a direct quote from the text supporting your answer.

## OUTPUT (JSON only)
```json
{{
  "reasoning": "<1-2 sentences>",
  "pathogenicity": "<Pathogenic|Likely pathogenic|Uncertain significance|Likely benign|Benign|Not specified>",
  "disease": "<full disease name, or 'Not specified'>",
  "tissue_affected": "<Cardiac|Skeletal|Both|Not specified>",
  "age_onset": "<Congenital|Pediatric|Adult|Not specified>",
  "inheritance": "<Autosomal Dominant|Autosomal Recessive|X-linked|De Novo|Sporadic|Not specified>",
  "patient_count": <integer>,
  "evidence_sentence": "<direct quote or 'Not specified'>"
}}
```"""
        user_message = f"## ARTICLE TEXT\n\n{context}\n\n---\n\nExtract for the target variant. Output ONLY the JSON object."
        if self.is_gemma:
            return (f"<start_of_turn>user\n{system_message}\n\n{user_message}"
                    f"<end_of_turn>\n<start_of_turn>model\n```json")
        return (f"<|begin_of_text|><|start_header_id|>system<|end_header_id|>\n"
                f"{system_message}<|eot_id|><|start_header_id|>user<|end_header_id|>\n"
                f"{user_message}<|eot_id|><|start_header_id|>assistant<|end_header_id|>\n```json")

    # ── 批次推論 ────────────────────────────────────────────────────────────
    def extract_batch(self, jobs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """jobs: [{article, variant}]，回傳對應的 clinical info dict 清單。"""
        prompts = [self._build_prompt(self._build_context(j["article"], j["variant"]),
                                      j["variant"]) for j in jobs]
        logger.info(f"vLLM 批次推論 {len(prompts)} 個 (文獻×變異) prompt …")
        outputs = self.llm.generate(prompts, self.sampling_params)
        results = []
        for o in outputs:
            samples = [c.text for c in o.outputs]
            parsed = [self._parse_json(s) for s in samples]
            parsed = [p for p in parsed if p]
            results.append(self._vote(parsed) if parsed else dict(EMPTY_RESULT))
        return results

    @staticmethod
    def _parse_json(text: str) -> Optional[dict]:
        clean = text.strip()
        for pre in ("```json", "```"):
            if clean.startswith(pre):
                clean = clean[len(pre):]
        if clean.endswith("```"):
            clean = clean[:-3]
        s, e = clean.find("{"), clean.rfind("}")
        if s == -1 or e == -1:
            return None
        for cand in (clean[s:e + 1], re.sub(r",\s*([}\]])", r"\1", clean[s:e + 1])):
            try:
                return json.loads(cand)
            except json.JSONDecodeError:
                continue
        return None

    @staticmethod
    def _vote(samples: List[dict]) -> dict:
        def pick(field):
            vals = [str(s.get(field, "")).strip() for s in samples
                    if str(s.get(field, "")).strip() not in ("", "Not specified")]
            return Counter(vals).most_common(1)[0][0] if vals else "Not specified"

        counts = []
        for s in samples:
            try:
                counts.append(int(s.get("patient_count", 0) or 0))
            except (ValueError, TypeError):
                pass
        counts.sort()
        out = {
            "pathogenicity": pick("pathogenicity"),
            "disease": pick("disease"),
            "tissue_affected": pick("tissue_affected"),
            "age_onset": pick("age_onset"),
            "inheritance": pick("inheritance"),
            "patient_count": counts[len(counts) // 2] if counts else 0,
            "reasoning": str(samples[0].get("reasoning", "")).strip(),
            "evidence_sentence": str(samples[0].get("evidence_sentence", "")).strip(),
        }
        return out


# ════════════════════════════════════════════════════════════════════════════
# Pipeline orchestration
# ════════════════════════════════════════════════════════════════════════════
def select_journal_articles(candidates: List[Dict[str, Any]],
                            pubtype_map: Dict[str, Dict[str, Any]],
                            max_articles: int,
                            include_reviews: bool) -> List[Dict[str, Any]]:
    selected = []
    for c in candidates:
        meta = pubtype_map.get(c["pmid"], {})
        pubtypes = meta.get("pubtypes", [])
        c["year"] = meta.get("year", "")
        is_journal = "Journal Article" in pubtypes
        is_review = "Review" in pubtypes
        if not is_journal:
            continue
        if is_review and not include_reviews:
            continue
        selected.append(c)
        if len(selected) >= max_articles:
            break
    return selected


def run(args):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    client = PubTatorClient()

    # Stage 1: 候選 PMID（多抓一些以容忍非 journal article 被過濾掉）
    buffer = max(args.max_articles * 3, args.max_articles + 30)
    logger.info("=" * 70)
    logger.info("Stage 1: 搜尋 @GENE_TTN 取候選 PMID")
    candidates = client.search_ttn_pmids(buffer)
    logger.info(f"  收集到 {len(candidates)} 個候選")

    # Stage 1b: 過濾 journal article
    logger.info("Stage 1b: esummary 取得 pubtype / 年份並過濾 journal article")
    pubtype_map = client.annotate_pubtypes([c["pmid"] for c in candidates])
    selected = select_journal_articles(candidates, pubtype_map, args.max_articles,
                                       args.include_reviews)
    logger.info(f"  選出 {len(selected)} 篇 journal article")
    if not selected:
        logger.error("沒有可處理的文獻，結束")
        return 1

    # Stage 2: 取變異 + 全文
    #   - 變異來源：PubTator 標註（優先全文 full=true，缺全文者退回 abstract 層級）
    #   - full=true 的 doc id 是 PMCID、批次會重排，用 article-id_pmid 對應回 PMID
    #   - 全文取不到時，LLM context 再試 EuropePMC，最後退回 abstract
    logger.info("=" * 70)
    logger.info("Stage 2: 抓全文/abstract biocjson + 解析 TTN 變異")
    sel_pmids = [c["pmid"] for c in selected]

    fulltext_docs: Dict[str, dict] = {}
    if not args.no_fulltext:
        logger.info("  批次抓 full=true 全文 …")
        fulltext_docs = client.fetch_fulltext_docs(sel_pmids, force=args.refresh)
        logger.info(f"  取得 {len(fulltext_docs)}/{len(sel_pmids)} 篇全文")

    # 缺全文的才抓 abstract 層級（省請求）
    missing = [p for p in sel_pmids if p not in fulltext_docs]
    abstract_docs = client.fetch_abstract_docs(missing, force=args.refresh) if missing else {}

    articles: List[Dict[str, Any]] = []
    for i, c in enumerate(selected, 1):
        pmid = c["pmid"]
        fdoc = fulltext_docs.get(pmid)
        adoc = abstract_docs.get(pmid)
        primary = fdoc or adoc
        text = extract_text_from_doc(primary)

        variants = extract_variants_from_doc(primary)

        full_text = text["full_text"]
        text_source = "pubtator_full" if fdoc else ("abstract" if adoc else "")
        # 有 abstract doc 但無全文 → 試 EuropePMC 補全文（僅供 LLM context）
        if not fdoc and not args.no_fulltext:
            eup = client.fetch_eupmc_fulltext(pmid, force=args.refresh)
            if eup and len(eup) > len(full_text):
                full_text = eup
                text_source = "europepmc"

        for v in variants.values():  # set → 排序後 list
            v["texts"] = sorted(v["texts"]) if isinstance(v["texts"], set) else v["texts"]

        articles.append({
            "pmid": pmid,
            "pmcid": text["pmcid"] or c.get("pmcid", ""),
            "title": text["title"] or c.get("title", ""),
            "abstract": text["abstract"],
            "full_text": full_text,
            "variants": variants,
            "journal": c.get("journal", ""),
            "year": c.get("year", ""),
            "score": c.get("score"),
            "text_source": text_source,
        })
        logger.info(f"  [{i}/{len(selected)}] PMID {pmid}: "
                    f"{len(variants)} 個 TTN 變異，context={len(full_text)} 字元 "
                    f"(來源={text_source})")

    n_variants = sum(len(a["variants"]) for a in articles)
    logger.info(f"共 {len(articles)} 篇、{n_variants} 個 (文獻×TTN變異) 待處理")

    # 組裝 jobs（每個 = 文獻 × 變異）
    jobs: List[Dict[str, Any]] = []
    for a in articles:
        for v in a["variants"].values():
            jobs.append({"article": a, "variant": v})

    # Stage 3+4: LLM（分批 + 檢查點，可續跑）
    def job_key(job: Dict[str, Any]) -> str:
        return f"{job['article']['pmid']}|{job['variant']['key']}"

    if args.no_llm or not jobs:
        if args.no_llm:
            logger.info("--no-llm：跳過 LLM 擷取，僅輸出變異清單")
        clinical = [dict(EMPTY_RESULT) for _ in jobs]
    else:
        logger.info("=" * 70)
        logger.info("Stage 3+4: 本地 LLM 判定致病性 + 擷取表型")

        checkpoint = CACHE_DIR / "llm_results.jsonl"
        if args.refresh_llm and checkpoint.exists():
            checkpoint.unlink()
        done: Dict[str, Dict[str, Any]] = {}
        if checkpoint.exists():
            for line in checkpoint.read_text(encoding="utf-8").splitlines():
                try:
                    rec = json.loads(line)
                    done[rec["key"]] = rec["info"]
                except (json.JSONDecodeError, KeyError):
                    continue
            logger.info(f"  檢查點已有 {len(done)} 筆結果，續跑剩餘部分")

        todo = [j for j in jobs if job_key(j) not in done]
        logger.info(f"  待推論 {len(todo)}/{len(jobs)}（batch={args.llm_batch_size}）")

        if todo:
            extractor = PathogenicVariantExtractor(model_name=args.llm_model,
                                                   tokenizer=args.llm_tokenizer)
            with open(checkpoint, "a", encoding="utf-8") as ckpt:
                for s in range(0, len(todo), args.llm_batch_size):
                    chunk = todo[s:s + args.llm_batch_size]
                    results = extractor.extract_batch(chunk)
                    for job, info in zip(chunk, results):
                        k = job_key(job)
                        done[k] = info
                        ckpt.write(json.dumps({"key": k, "info": info},
                                              ensure_ascii=False) + "\n")
                    ckpt.flush()
                    logger.info(f"  進度 {min(s + len(chunk), len(todo))}/{len(todo)} "
                                f"已寫入檢查點")

        clinical = [done.get(job_key(j), dict(EMPTY_RESULT)) for j in jobs]

    # Stage 5: 輸出 CSV
    logger.info("=" * 70)
    logger.info("Stage 5: 彙整輸出 CSV")
    rows = []
    for job, info in zip(jobs, clinical):
        a, v = job["article"], job["variant"]
        rows.append({
            "pmid": a["pmid"],
            "pmcid": a.get("pmcid", ""),
            "year": a.get("year", ""),
            "journal": a.get("journal", ""),
            "title": a.get("title", ""),
            "variant_text": "; ".join(v.get("texts", [])),
            "hgvs": v.get("hgvs", ""),
            "rsid": v.get("rsid", ""),
            "gene": "TTN",
            "mention_count": v.get("mention_count", 0),
            "pathogenicity": info.get("pathogenicity", "Not specified"),
            "disease": info.get("disease", "Not specified"),
            "tissue_affected": info.get("tissue_affected", "Not specified"),
            "age_onset": info.get("age_onset", "Not specified"),
            "inheritance": info.get("inheritance", "Not specified"),
            "patient_count": info.get("patient_count", 0),
            "evidence_sentence": info.get("evidence_sentence", ""),
            "reasoning": info.get("reasoning", ""),
            "text_source": a.get("text_source", ""),
            "pubmed_link": f"https://pubmed.ncbi.nlm.nih.gov/{a['pmid']}/",
        })

    import pandas as pd
    df = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False, encoding="utf-8-sig")
    logger.info(f"✓ 已輸出 {len(df)} 列 → {out_path}")
    if not df.empty and not args.no_llm:
        patho = df[df["pathogenicity"].isin(["Pathogenic", "Likely pathogenic"])]
        logger.info(f"  其中致病/可能致病：{len(patho)} 列")
    return 0


def main():
    p = argparse.ArgumentParser(description="PubTator3 → TTN 致病變異彙整（獨立版）")
    p.add_argument("--max-articles", type=int, default=100, help="處理的 journal article 篇數")
    p.add_argument("--out", type=str, default=str(DEFAULT_OUT), help="輸出 CSV 路徑")
    p.add_argument("--no-llm", action="store_true", help="跳過 LLM，只抓變異（不需 GPU）")
    p.add_argument("--no-fulltext", action="store_true",
                   help="不抓全文（full=true / EuropePMC），變異與 context 只用 abstract")
    p.add_argument("--include-reviews", action="store_true", help="一併納入 Review 類型")
    p.add_argument("--refresh", action="store_true", help="忽略快取，重新下載 biocjson")
    p.add_argument("--llm-model", type=str, default=None,
                   help="覆寫 LLM 模型（預設用 config）。GGUF 用 'repo:file.gguf'，"
                        "例：unsloth/gemma-4-31B-it-GGUF:gemma-4-31B-it-Q8_0.gguf")
    p.add_argument("--llm-tokenizer", type=str, default=None,
                   help="覆寫 tokenizer（GGUF 預設用 google/gemma-4-31B-it）")
    p.add_argument("--llm-batch-size", type=int, default=32,
                   help="每批送進 vLLM 的 (文獻×變異) 數量（越小越省 KV cache）")
    p.add_argument("--refresh-llm", action="store_true",
                   help="清掉 LLM 檢查點，重新推論所有變異")
    args = p.parse_args()
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
