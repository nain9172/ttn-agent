"""把輸入變異對回 HGMD / PubTator / PubMind 三個知識庫，取出「這顆變異在該資料庫裡登記的疾病」。

跟 Track B 的類比推理（mechanistic_inference）不同：這裡不做任何推論，
只做查表，命中就把資料庫原本寫的疾病照實列出來。

比對層級分兩種：
  exact        — chrom-pos-ref-alt 完全相同（同一顆變異）
  same_position — 同一個 genomic 位置但 allele 不同 / 或該位置落在 indel 區間內

HGMD 原始 CSV 只有 startCoord/endCoord 沒有 REF/ALT，所以 exact 比對靠
acc_num 接到已處理過的 KB（HGMD_TTN_processed）拿 REF/ALT；沒有 KB 對應的
（多為 frameshift / gross indel）只能做 same_position。

PubTator 原本完全沒有座標，靠 utils.hgvs_resolver 離線把 c./p. HGVS 換成座標，
結果快取在 outputs/ 以免每次重算。data/pubtator/ 這一版另外帶了 ClinVar 標註
欄位，唯一命中且落在 chr2 的列直接有權威座標，優先採用；其餘仍走 resolver。
"""

from __future__ import annotations

import json
import logging
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from utils.pubmind_lookup import DEFAULT_PUBMIND_JSON, PubMindLookup
from utils.pubmind_lookup import disease_counts as pubmind_disease_counts

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA = PROJECT_ROOT / "data"

DEFAULT_HGMD_CSV = DATA / "HGMD_latest.csv"
DEFAULT_HGMD_KB = DATA / "HGMD_TTN_processed.xlsx"
DEFAULT_PUBTATOR_CSV = DATA / "pubtator" / "pubtator_ttn_variants_clinvar_annotated.csv"
DEFAULT_CACHE = PROJECT_ROOT / "outputs" / "pubtator_hgvs_coord_cache.json"

# indel 區間超過這個長度就不做 same_position 比對（gross deletion 動輒數 kb，
# 落在區間內不代表跟輸入變異有關）
MAX_SPAN_FOR_POSITION_MATCH = 50

_NULLS = {"", "null", "none", "nan", "n/a", "na", "not specified", "-"}


def _clean(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if text.lower() in _NULLS:
        return None
    return text


def _as_int(value: Any) -> Optional[int]:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _clean_id(value: Any) -> Optional[str]:
    """pandas 會把 PMID / ClinVar ID 這類整數欄位讀成 float（29382405.0），去掉小數尾。"""
    text = _clean(value)
    if text is None:
        return None
    if text.endswith(".0") and text[:-2].isdigit():
        return text[:-2]
    return text


@dataclass
class DatabaseHits:
    """單一變異的查表結果。"""

    hgmd_exact: List[Dict[str, Any]] = field(default_factory=list)
    hgmd_same_position: List[Dict[str, Any]] = field(default_factory=list)
    pubtator_exact: List[Dict[str, Any]] = field(default_factory=list)
    pubtator_same_position: List[Dict[str, Any]] = field(default_factory=list)
    pubmind_exact: List[Dict[str, Any]] = field(default_factory=list)
    # PubMind 資料沒載入時為 False，報告就不顯示 PubMind 卡片
    pubmind_checked: bool = False

    @property
    def has_exact(self) -> bool:
        return bool(self.hgmd_exact or self.pubtator_exact)

    def hgmd_diseases(self) -> List[Tuple[str, int]]:
        return _disease_counts(self.hgmd_exact, "disease")

    def pubtator_diseases(self) -> List[Tuple[str, int]]:
        return _disease_counts(self.pubtator_exact, "disease")

    def pubmind_diseases(self) -> List[Tuple[str, int]]:
        return pubmind_disease_counts(self.pubmind_exact)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "hgmd_exact": self.hgmd_exact,
            "hgmd_same_position": self.hgmd_same_position,
            "pubtator_exact": self.pubtator_exact,
            "pubtator_same_position": self.pubtator_same_position,
            "hgmd_diseases": self.hgmd_diseases(),
            "pubtator_diseases": self.pubtator_diseases(),
            "pubmind_exact": self.pubmind_exact,
            "pubmind_checked": self.pubmind_checked,
            "pubmind_diseases": self.pubmind_diseases(),
        }


def _disease_counts(records: List[Dict[str, Any]], key: str) -> List[Tuple[str, int]]:
    """把命中紀錄裡的疾病名稱去重計數，多的排前面。"""
    counter: Counter = Counter()
    for rec in records:
        disease = _clean(rec.get(key))
        if disease:
            counter[disease] += 1
    return counter.most_common()


class VariantDatabaseLookup:
    """HGMD + PubTator 的變異→疾病查表器（建一次索引，之後 O(1) 查詢）。"""

    def __init__(
        self,
        hgmd_csv: Path = DEFAULT_HGMD_CSV,
        hgmd_kb: Path = DEFAULT_HGMD_KB,
        pubtator_csv: Path = DEFAULT_PUBTATOR_CSV,
        cache_path: Optional[Path] = DEFAULT_CACHE,
        resolver: Any = None,
        pubmind_json: Optional[Path] = DEFAULT_PUBMIND_JSON,
    ) -> None:
        self.hgmd_csv = Path(hgmd_csv)
        self.hgmd_kb = Path(hgmd_kb)
        self.pubtator_csv = Path(pubtator_csv)
        self.cache_path = Path(cache_path) if cache_path else None
        self._resolver = resolver

        self.hgmd_by_variant: Dict[str, List[Dict[str, Any]]] = {}
        self.hgmd_by_position: Dict[int, List[Dict[str, Any]]] = {}
        self.pubtator_by_variant: Dict[str, List[Dict[str, Any]]] = {}
        self.pubtator_by_position: Dict[int, List[Dict[str, Any]]] = {}

        self.stats: Dict[str, Any] = {}

        self._load_hgmd()
        self._load_pubtator()
        self.pubmind = PubMindLookup(pubmind_json) if pubmind_json else None

    # ── HGMD ──────────────────────────────────────────────────────────────

    def _load_hgmd(self) -> None:
        if not self.hgmd_csv.exists():
            logger.warning("HGMD CSV not found: %s", self.hgmd_csv)
            return

        df = pd.read_csv(self.hgmd_csv, low_memory=False)

        # acc_num → variant_id（REF/ALT 只有已處理過的 KB 有）
        acc_to_variant: Dict[str, str] = {}
        if self.hgmd_kb.exists():
            try:
                kb = pd.read_excel(self.hgmd_kb)
                for row in kb.itertuples():
                    acc = _clean(getattr(row, "ID", None))
                    vid = _clean(getattr(row, "input", None))
                    if acc and vid:
                        acc_to_variant[acc] = vid
            except Exception as e:  # noqa: BLE001
                logger.warning("Could not read HGMD KB for REF/ALT mapping: %s", e)

        n_exact = 0
        for row in df.itertuples():
            acc = _clean(getattr(row, "acc_num", None))
            start = _as_int(getattr(row, "startCoord", None))
            end = _as_int(getattr(row, "endCoord", None)) or start

            record = {
                "disease": _clean(getattr(row, "disease", None)) or "Not specified",
                "mutype": _clean(getattr(row, "mutype", None)),
                "hgvs_c": _clean(getattr(row, "hgvs", None)),
                "refseq": _clean(getattr(row, "refseq", None)),
                "amino": _clean(getattr(row, "amino", None)),
                "codon": _clean(getattr(row, "codon", None)),
                "descr": _clean(getattr(row, "descr", None)),
                "inheritance": _clean(getattr(row, "expected_inheritance", None)),
                "tag": _clean(getattr(row, "tag", None)),
                "pmid": _clean_id(getattr(row, "pmid", None)),
                "author": _clean(getattr(row, "author", None)),
                "year": _clean_id(getattr(row, "year", None)),
                "journal": _clean(getattr(row, "fullname", None)),
                "title": _clean(getattr(row, "title", None)),
                "dbsnp": _clean(getattr(row, "dbsnp", None)),
                "clinvar_id": _clean_id(getattr(row, "clinvarID", None)),
                "clinvar_clnsig": _clean(getattr(row, "clinvar_clnsig", None)),
                "acc_num": acc,
                "start": start,
                "end": end,
                "comments": _clean(getattr(row, "comments", None)),
            }

            variant_id = acc_to_variant.get(acc) if acc else None
            if variant_id:
                record["variant_id"] = variant_id
                self.hgmd_by_variant.setdefault(variant_id, []).append(record)
                n_exact += 1

            if start is None:
                continue
            span = (end or start) - start
            if 0 <= span <= MAX_SPAN_FOR_POSITION_MATCH:
                for pos in range(start, (end or start) + 1):
                    self.hgmd_by_position.setdefault(pos, []).append(record)
            else:
                self.hgmd_by_position.setdefault(start, []).append(record)

        self.stats["hgmd_rows"] = len(df)
        self.stats["hgmd_with_alleles"] = n_exact
        self.stats["hgmd_variants"] = len(self.hgmd_by_variant)
        logger.info(
            "HGMD loaded: %d entries (%d with REF/ALT → %d unique variants)",
            len(df),
            n_exact,
            len(self.hgmd_by_variant),
        )

    # ── PubTator ──────────────────────────────────────────────────────────

    @property
    def resolver(self) -> Any:
        """延後建立 HGVS resolver（載入轉錄本模型有成本，沒用到就別建）。"""
        if self._resolver is None:
            from utils.hgvs_resolver import HgvsResolver

            import config

            self._resolver = HgvsResolver(
                model_path=config.TTN_TRANSCRIPT_MODEL_PATH,
                fasta_path=getattr(config, "REFERENCE_GENOME_PATH", None),
                allow_network=False,
            )
        return self._resolver

    def _load_coord_cache(self) -> Dict[str, Any]:
        if self.cache_path and self.cache_path.exists():
            try:
                with open(self.cache_path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception as e:  # noqa: BLE001
                logger.warning("Could not read coord cache: %s", e)
        return {}

    def _save_coord_cache(self, cache: Dict[str, Any]) -> None:
        if not self.cache_path:
            return
        try:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.cache_path, "w", encoding="utf-8") as f:
                json.dump(cache, f)
        except Exception as e:  # noqa: BLE001
            logger.warning("Could not write coord cache: %s", e)

    def _load_pubtator(self) -> None:
        if not self.pubtator_csv.exists():
            logger.warning("PubTator CSV not found: %s", self.pubtator_csv)
            return

        df = pd.read_csv(self.pubtator_csv, low_memory=False)
        cache = self._load_coord_cache()
        cache_dirty = False

        # 快取沒蓋到的列才需要 resolver，先確認它能建起來（建不起來就只用快取）
        try:
            _ = self.resolver
        except Exception as e:  # noqa: BLE001
            logger.warning("HGVS resolver unavailable, PubTator lookup limited to cache: %s", e)

        # 同一篇論文裡的 p.HGVS 可以幫 c.HGVS 判斷該用哪個轉錄本編號
        p_hint: Dict[Tuple[Any, str], str] = {}
        for row in df.itertuples():
            hgvs = _clean(getattr(row, "hgvs", None))
            rsid = _clean(getattr(row, "rsid", None))
            if hgvs and hgvs.lower().startswith("p.") and rsid:
                p_hint.setdefault((getattr(row, "pmid", None), rsid.split("(")[0]), hgvs)

        n_resolved = 0
        n_failed = 0
        for row in df.itertuples():
            hgvs = _clean(getattr(row, "hgvs", None))
            variant_text = _clean(getattr(row, "variant_text", None))
            rsid = _clean(getattr(row, "rsid", None))
            rsid = rsid.split("(")[0] if rsid else None

            queries = [q for q in (hgvs, variant_text) if q]
            if not queries:
                continue

            hint = p_hint.get((getattr(row, "pmid", None), rsid)) if rsid else None
            if hint and hgvs and hgvs.lower().startswith("p."):
                hint = None

            cache_key = "|".join(queries) + f"#{hint or ''}"
            if cache_key in cache:
                coord = cache[cache_key]
            else:
                coord = self._resolve_queries(queries, hint)
                cache[cache_key] = coord
                cache_dirty = True

            # ClinVar 唯一命中的座標比離線 resolver 權威，優先採用。限定 chr2：
            # 用 rsID 對到 ClinVar 時偶爾會連到別條染色體的變異，那是錯配。
            cv_pos = None
            if str(_clean(getattr(row, "chromosome", None)) or "") == "2":
                cv_pos = _as_int(getattr(row, "position", None))

            if not coord and cv_pos is None:
                n_failed += 1
                continue
            n_resolved += 1

            # variant_id 需要 REF/ALT，只有 resolver 給得出來；ClinVar 欄位只有位置。
            pos = cv_pos if cv_pos is not None else (coord or {}).get("pos")
            pos_source = "clinvar_annotation" if cv_pos is not None else "hgvs_resolver"

            record = {
                "variant_id": (coord or {}).get("variant_id"),
                "hgvs_c": hgvs,
                "variant_text": variant_text,
                "rsid": rsid,
                "disease": _clean(getattr(row, "disease", None)) or "Not specified",
                "pathogenicity": _clean(getattr(row, "pathogenicity", None)),
                "tissue": _clean(getattr(row, "tissue_affected", None)),
                "age_onset": _clean(getattr(row, "age_onset", None)),
                "inheritance": _clean(getattr(row, "inheritance", None)),
                "patient_count": _clean(getattr(row, "patient_count", None)),
                "evidence_sentence": _clean(getattr(row, "evidence_sentence", None)),
                "pmid": _clean_id(getattr(row, "pmid", None)),
                "year": _clean_id(getattr(row, "year", None)),
                "journal": _clean(getattr(row, "journal", None)),
                "title": _clean(getattr(row, "title", None)),
                "pubmed_link": _clean(getattr(row, "pubmed_link", None)),
                "resolved_by": (coord or {}).get("resolved_by"),
                "resolve_confidence": (coord or {}).get("confidence"),
                "pos_source": pos_source,
                "clinvar_variation_id": _clean_id(getattr(row, "ClinVar_Variation_ID", None)),
            }

            vid = (coord or {}).get("variant_id")
            if vid:
                self.pubtator_by_variant.setdefault(vid, []).append(record)
            if pos:
                self.pubtator_by_position.setdefault(int(pos), []).append(record)

        if cache_dirty:
            self._save_coord_cache(cache)

        self.stats["pubtator_rows"] = len(df)
        self.stats["pubtator_resolved_rows"] = n_resolved
        self.stats["pubtator_unresolved_rows"] = n_failed
        self.stats["pubtator_variants"] = len(self.pubtator_by_variant)
        logger.info(
            "PubTator loaded: %d rows (%d resolved → %d unique variants, %d unresolvable)",
            len(df),
            n_resolved,
            len(self.pubtator_by_variant),
            n_failed,
        )

    def _resolve_queries(self, queries: List[str], hint: Optional[str]) -> Optional[Dict[str, Any]]:
        for q in queries:
            try:
                resolved = self.resolver.resolve(q, hgvs_p=hint)
            except Exception:  # noqa: BLE001  解析失敗的原因不影響查表，跳過即可
                continue
            if resolved is None:
                continue
            return {
                "variant_id": resolved.variant_id,
                "pos": resolved.pos,
                "ref": resolved.ref,
                "alt": resolved.alt,
                "resolved_by": resolved.resolved_by,
                "confidence": resolved.confidence,
            }
        return None

    # ── 查詢 ──────────────────────────────────────────────────────────────

    def lookup(
        self,
        variant_id: Optional[str] = None,
        pos: Optional[int] = None,
        chrom: str = "2",
        ref: Optional[str] = None,
        alt: Optional[str] = None,
    ) -> DatabaseHits:
        """查一顆變異。給 variant_id 或 (pos, ref, alt) 皆可。"""
        if not variant_id and pos is not None and ref and alt:
            variant_id = f"{chrom}-{pos}-{ref}-{alt}"
        if variant_id and pos is None:
            parts = variant_id.split("-")
            if len(parts) >= 2:
                pos = _as_int(parts[1])

        hits = DatabaseHits()

        if variant_id:
            hits.hgmd_exact = list(self.hgmd_by_variant.get(variant_id, []))
            hits.pubtator_exact = list(self.pubtator_by_variant.get(variant_id, []))
            parts = variant_id.split("-")
            if self.pubmind is not None and self.pubmind.loaded and len(parts) == 4:
                hits.pubmind_checked = True
                hits.pubmind_exact = self.pubmind.lookup(int(parts[1]), parts[2], parts[3])

        if pos is not None:
            exact_accs = {r.get("acc_num") for r in hits.hgmd_exact}
            hits.hgmd_same_position = [
                r for r in self.hgmd_by_position.get(pos, [])
                if r.get("acc_num") not in exact_accs
            ]
            exact_keys = {(r.get("pmid"), r.get("hgvs_c")) for r in hits.pubtator_exact}
            hits.pubtator_same_position = [
                r for r in self.pubtator_by_position.get(pos, [])
                if (r.get("pmid"), r.get("hgvs_c")) not in exact_keys
            ]

        return hits


_SHARED: Optional[VariantDatabaseLookup] = None


def get_shared_lookup(**kwargs: Any) -> VariantDatabaseLookup:
    """單例，讓 main.py 每次執行只建一次索引。"""
    global _SHARED
    if _SHARED is None:
        _SHARED = VariantDatabaseLookup(**kwargs)
    return _SHARED


if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    lookup = VariantDatabaseLookup()
    print(json.dumps(lookup.stats, indent=2))

    for vid in sys.argv[1:]:
        hits = lookup.lookup(variant_id=vid)
        print(f"\n=== {vid}")
        print(f"  HGMD exact       : {hits.hgmd_diseases()}")
        print(f"  HGMD same pos    : {len(hits.hgmd_same_position)} entries")
        print(f"  PubTator exact   : {hits.pubtator_diseases()}")
        print(f"  PubTator same pos: {len(hits.pubtator_same_position)} entries")
