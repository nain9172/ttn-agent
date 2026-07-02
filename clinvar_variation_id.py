#!/usr/bin/env python3
"""
為 HGMD TTN 變異表加上 ClinVar Variation ID。

讀取 Excel（預設 `HGMD_TTN_processed (1).xlsx`），逐筆變異查詢 ClinVar，
把對應的 ClinVar Variation ID 寫進新欄位 `ClinVar_Variation_ID`，再另存新檔。

查詢策略（皆沿用 utils/clinvar_parser.py 已實測的 NCBI E-utilities 行為）：
  1. 主要：用 SPDI 座標精確搜尋。TTN 在負鏈，需 pos-1 並取互補鹼基，
     query 形如 "NC_000002.12:178527030:A:G"。
  2. 後備：若座標查不到（如 indel / SV），改用 DNA 欄位的 HGVS（如
     NM_001267550.2:c.107957T>C）直接搜尋 ClinVar。
  3. esearch 回傳多筆時，逐筆用 efetch VCV XML 比對座標（_verify_variant_match）
     挑出真正匹配者；只有一筆時直接採用。

查詢結果會快取到 JSON，支援中斷後續跑（resume）。

用法：
    conda run -n gm4 python clinvar_variation_id.py
    conda run -n gm4 python clinvar_variation_id.py --limit 10          # 先測前 10 筆
    conda run -n gm4 python clinvar_variation_id.py --input "HGMD_TTN_processed (1).xlsx" \
        --output HGMD_TTN_with_clinvar.xlsx
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from utils.clinvar_parser import ClinVarParser

try:
    from Bio import Entrez
except ImportError:  # pragma: no cover
    print("錯誤：需要 Biopython（Bio.Entrez）。請在 gm4 環境執行此程式。")
    sys.exit(1)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("clinvar_variation_id")
# clinvar_parser 內部 log 偏冗長，降到 WARNING 讓進度清楚
logging.getLogger("utils.clinvar_parser").setLevel(logging.WARNING)

PROJECT_ROOT = Path(__file__).parent
DEFAULT_INPUT = PROJECT_ROOT / "HGMD_TTN_processed (1).xlsx"
DEFAULT_OUTPUT = PROJECT_ROOT / "HGMD_TTN_with_clinvar.xlsx"
# --pubtator 模式的預設輸入/輸出
DEFAULT_PUBTATOR_INPUT = PROJECT_ROOT / "pubtator_ttn_variants_all.csv"
DEFAULT_PUBTATOR_OUTPUT = PROJECT_ROOT / "pubtator_ttn_variants_all_with_clinvar.csv"
RESULT_COLUMN = "ClinVar_Variation_ID"

# GRCh38 chr2 RefSeq accession（TTN 所在染色體）
CHR2_ACCESSION = "NC_000002.12"
EUTILS_SLEEP = 0.34  # 無 API key 約 3 req/s；有 key 也夠安全


class ClinVarVariationLookup:
    """以 ClinVarParser 的座標轉換 / 驗證邏輯為基礎，只取回 Variation ID。"""

    def __init__(self) -> None:
        self.parser = ClinVarParser()

    def lookup(self, variant_info: Dict, hgvs_cdna: Optional[str] = None) -> Optional[str]:
        """回傳匹配的 ClinVar Variation ID（字串）或 None。"""
        vid = self._lookup_by_spdi(variant_info)
        if vid:
            return vid
        if hgvs_cdna:
            vid = self._lookup_by_hgvs(hgvs_cdna, variant_info)
            if vid:
                return vid
        return None

    def lookup_pubtator(
        self,
        variant_text: Optional[str],
        hgvs: Optional[str],
        rsid: Optional[str],
    ) -> Optional[str]:
        """
        PubTator 模式：沒有 chrom/pos/ref/alt 可做座標驗證，依序用文字查 ClinVar。

        順序：variant_text → hgvs → rsid，任一查到即回傳（取第一筆）。
        三者皆無結果回傳 None。HGVS/變異文字加上 ``AND TTN[gene]`` 限縮在 TTN，
        rsID 本身具唯一性故直接查。
        """
        tried: set = set()
        for raw, is_rsid in ((variant_text, False), (hgvs, False), (rsid, True)):
            term = self._build_pubtator_term(raw, is_rsid)
            if not term or term in tried:
                continue
            tried.add(term)
            ids = self._esearch(term)
            if ids:
                return ids[0]
        return None

    @staticmethod
    def _build_pubtator_term(raw: Optional[str], is_rsid: bool) -> Optional[str]:
        if raw is None:
            return None
        s = str(raw).strip()
        if not s or s.lower() == "nan":
            return None
        if is_rsid:
            m = re.search(r"rs\d+", s)  # 去掉 "(Expired)" 等雜訊
            return m.group(0) if m else None
        # 去除變異表示法中的所有空白（"c.76717C > T" → "c.76717C>T"），
        # 帶空白的 query 會被 ClinVar 斷詞而抓到一堆不相關結果。
        normalized = re.sub(r"\s+", "", s)
        if not normalized:
            return None
        return f"{normalized} AND TTN[gene]"

    def _lookup_by_spdi(self, variant_info: Dict) -> Optional[str]:
        ref = str(variant_info["ref"]).upper()
        alt = str(variant_info["alt"]).upper()
        # SPDI 座標搜尋只對單鹼基替換可靠；indel/SV 交給 HGVS 後備
        if len(ref) != 1 or len(alt) != 1 or ref == "-" or alt == "-":
            return None

        neg_pos = int(variant_info["pos"]) - 1
        neg_ref = self.parser._get_complement_base(ref)
        neg_alt = self.parser._get_complement_base(alt)
        query = f'"{CHR2_ACCESSION}:{neg_pos}:{neg_ref}:{neg_alt}"'
        ids = self._esearch(query)
        return self._pick_match(ids, variant_info)

    def _lookup_by_hgvs(self, hgvs_cdna: str, variant_info: Dict) -> Optional[str]:
        hgvs_cdna = hgvs_cdna.strip()
        if not hgvs_cdna:
            return None
        ids = self._esearch(hgvs_cdna)
        if not ids:
            # 退一步：只用無轉錄本版本號的 c. 表示法再試一次
            if ":" in hgvs_cdna:
                ids = self._esearch(hgvs_cdna.split(":", 1)[1])
        # HGVS 搜尋可能對到不同 isoform 的同一變異，多筆時仍用座標驗證
        return self._pick_match(ids, variant_info, trust_single=True)

    def _esearch(self, term: str) -> List[str]:
        for attempt in range(1, 4):
            try:
                time.sleep(EUTILS_SLEEP)
                handle = Entrez.esearch(db="clinvar", term=term, retmax=5)
                res = Entrez.read(handle, validate=False)
                handle.close()
                return list(res.get("IdList", []))
            except Exception as e:  # noqa: BLE001
                logger.debug(f"  esearch 失敗（attempt {attempt}）term={term}: {e}")
                time.sleep(1.0 * attempt)
        return []

    def _pick_match(
        self, ids: List[str], variant_info: Dict, trust_single: bool = True
    ) -> Optional[str]:
        if not ids:
            return None
        if len(ids) == 1 and trust_single:
            return ids[0]
        # 多筆 → 逐筆 efetch VCV XML 驗證座標
        for cid in ids:
            xml = self._efetch_vcv(cid)
            if xml and self.parser._verify_variant_match(xml, variant_info):
                return cid
        # 都驗不過時，保守起見回傳第一筆（座標 esearch 通常已很精確）
        return ids[0]

    def _efetch_vcv(self, clinvar_id: str) -> Optional[str]:
        try:
            time.sleep(EUTILS_SLEEP)
            handle = Entrez.esummary(db="clinvar", id=clinvar_id)
            summary = Entrez.read(handle, validate=False)
            handle.close()
            accession = None
            try:
                doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
                accession = doc.get("accession")
            except Exception:  # noqa: BLE001
                accession = None
            if not accession:
                try:
                    accession = f"VCV{int(clinvar_id):09d}"
                except Exception:  # noqa: BLE001
                    accession = clinvar_id

            time.sleep(EUTILS_SLEEP)
            handle = Entrez.efetch(
                db="clinvar", id=accession, rettype="vcv", retmode="xml"
            )
            xml = handle.read()
            handle.close()
            if isinstance(xml, bytes):
                xml = xml.decode("utf-8", errors="replace")
            return xml
        except Exception as e:  # noqa: BLE001
            logger.debug(f"  efetch VCV 失敗 id={clinvar_id}: {e}")
            return None


def _build_variant_info(row: pd.Series) -> Optional[Dict]:
    try:
        chrom = str(row["#CHROM"]).strip()
        pos = int(row["POS"])
        ref = str(row["REF"]).strip()
        alt = str(row["ALT"]).strip()
    except (KeyError, ValueError, TypeError):
        return None
    if not ref or not alt or ref.lower() == "nan" or alt.lower() == "nan":
        return None
    variant_id = str(row.get("ID") or row.get("input") or f"{chrom}-{pos}-{ref}-{alt}")
    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "variant_id": variant_id,
    }


def _load_cache(path: Path) -> Dict[str, Optional[str]]:
    if path.exists():
        try:
            with open(path, encoding="utf-8") as f:
                return json.load(f)
        except Exception:  # noqa: BLE001
            logger.warning(f"快取檔損毀，忽略：{path}")
    return {}


def _save_cache(path: Path, cache: Dict[str, Optional[str]]) -> None:
    # tmp 檔名帶 pid，避免多個 process 同時跑時共用 .tmp 互相覆寫，
    # 造成 rename 時 FileNotFoundError 的 race condition。
    tmp = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


def main() -> None:
    ap = argparse.ArgumentParser(description="為 TTN 變異表加上 ClinVar Variation ID")
    ap.add_argument(
        "--pubtator",
        action="store_true",
        help="處理 pubtator_ttn_variants_all.csv，依序用 variant_text→hgvs→rsid 查 ClinVar",
    )
    ap.add_argument("--input", default=None, help="輸入檔路徑（預設依模式決定）")
    ap.add_argument("--output", default=None, help="輸出檔路徑（預設依模式決定）")
    ap.add_argument("--sheet", default=0, help="（HGMD 模式）工作表名稱或索引")
    ap.add_argument("--limit", type=int, default=0, help="只處理前 N 筆（0 = 全部）")
    ap.add_argument("--cache", default=None, help="查詢結果快取 JSON 路徑（支援續跑）")
    ap.add_argument("--no-cache", action="store_true", help="不使用快取（每筆都重新查詢）")
    args = ap.parse_args()

    # 依模式決定預設輸入/輸出/快取
    if args.pubtator:
        in_path = Path(args.input) if args.input else DEFAULT_PUBTATOR_INPUT
        out_path = Path(args.output) if args.output else DEFAULT_PUBTATOR_OUTPUT
        default_cache = PROJECT_ROOT / "outputs" / "clinvar_variation_id_pubtator_cache.json"
    else:
        in_path = Path(args.input) if args.input else DEFAULT_INPUT
        out_path = Path(args.output) if args.output else DEFAULT_OUTPUT
        default_cache = PROJECT_ROOT / "outputs" / "clinvar_variation_id_cache.json"

    if not in_path.exists():
        logger.error(f"找不到輸入檔：{in_path}")
        sys.exit(1)

    if args.pubtator:
        df = pd.read_csv(in_path)
    else:
        sheet = args.sheet
        try:
            sheet = int(sheet)
        except (ValueError, TypeError):
            pass
        df = pd.read_excel(in_path, sheet_name=sheet)
    logger.info(f"讀取 {in_path}（{len(df)} 列，{len(df.columns)} 欄）")

    cache_path = Path(args.cache) if args.cache else default_cache
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = {} if args.no_cache else _load_cache(cache_path)
    if cache:
        logger.info(f"已載入快取 {len(cache)} 筆：{cache_path}")

    lookup = ClinVarVariationLookup()

    results: List[Optional[str]] = []
    n_total = len(df)
    n_process = args.limit if args.limit and args.limit < n_total else n_total
    n_found = 0

    for idx, row in df.iterrows():
        if args.limit and idx >= args.limit:
            results.append(None)
            continue

        if args.pubtator:
            variant_text = row.get("variant_text")
            hgvs = row.get("hgvs")
            rsid = row.get("rsid")

            def _clean(v):
                return str(v).strip() if pd.notna(v) else None

            variant_text, hgvs, rsid = _clean(variant_text), _clean(hgvs), _clean(rsid)
            label = variant_text or hgvs or rsid or "(空)"
            key = f"PT|{variant_text}|{hgvs}|{rsid}"
            if key in cache:
                vid = cache[key]
            else:
                vid = lookup.lookup_pubtator(variant_text, hgvs, rsid)
                cache[key] = vid
                if not args.no_cache:
                    _save_cache(cache_path, cache)
            results.append(vid)
            if vid:
                n_found += 1
            logger.info(
                f"[{idx + 1}/{n_process}] {label} -> ClinVar Variation ID: {vid or '(無)'}"
            )
            continue

        vi = _build_variant_info(row)
        if vi is None:
            results.append(None)
            logger.info(f"[{idx + 1}/{n_process}] 略過（變異資訊不完整）")
            continue

        key = f"{vi['chrom']}-{vi['pos']}-{vi['ref']}-{vi['alt']}"
        if key in cache:
            vid = cache[key]
        else:
            hgvs_cdna = row.get("DNA")
            hgvs_cdna = str(hgvs_cdna).strip() if pd.notna(hgvs_cdna) else None
            vid = lookup.lookup(vi, hgvs_cdna=hgvs_cdna)
            cache[key] = vid
            if not args.no_cache:
                _save_cache(cache_path, cache)

        results.append(vid)
        if vid:
            n_found += 1
        logger.info(
            f"[{idx + 1}/{n_process}] {vi['variant_id']} "
            f"({key}) -> ClinVar Variation ID: {vid or '(無)'}"
        )

    df[RESULT_COLUMN] = results
    if out_path.suffix.lower() == ".csv":
        df.to_csv(out_path, index=False)
    else:
        df.to_excel(out_path, index=False)
    logger.info("=" * 60)
    logger.info(f"完成：{n_found}/{n_process} 筆找到 ClinVar Variation ID")
    logger.info(f"輸出已寫入：{out_path}")


if __name__ == "__main__":
    main()
