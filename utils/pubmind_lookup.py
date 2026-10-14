"""把變異對回 PubMind 匯出的 TTN 文獻變異（data/pubmind_ttn.json）。

跟 variant_db_lookup 的 HGMD / PubTator 一樣只做查表，命中就照 PubMind 原本
登記的疾病、致病性、文獻列出來。

座標慣例：報告的 variant ID ref/alt 是 TTN 轉錄方向（負股），PubMind 的
parsed_variants 是正股 genomic allele，所以比對前要取互補。只有 aa_change 的紀錄
會列出所有可能的 codon 改變（ref/alt 為 3 bp），其中只差一個鹼基的就是能造成該
胺基酸改變的 SNV，標成 protein 層級；dna_change 直接解出的單鹼基則是 dna 層級。
沒有座標的紀錄（PubMind 自己也解不出 transcript）不納入比對。
"""

from __future__ import annotations

import ast
import json
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PUBMIND_JSON = PROJECT_ROOT / "data" / "pubmind_ttn.json"

COMPLEMENT = str.maketrans("ACGT", "TGCA")
RE_COUNTED = re.compile(r"\s*([^,]+?)\s*\((\d+)\)")
RE_PMID = re.compile(r"PMID:(\d+)")
RE_REFERENCE = re.compile(r"^(?:(PMC\d+))?\(?(?:PMID:(\d+))?\)?(?::(.*))?$")


def _blank(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text in ("-", "nan", "None") else text


def _counted(field: str) -> Counter:
    """'dilated cardiomyopathy (3), heart failure (1)' → Counter。"""
    counted: Counter = Counter()
    for name, n in RE_COUNTED.findall(field or ""):
        if name.strip() != "-":
            counted[name.strip()] += int(n)
    return counted


def _diseases(entry: Dict) -> Counter:
    counted = _counted(entry.get("MONDO_name_counted") or "")
    if not counted:
        # 沒有 MONDO 正規化的紀錄就退回原始 disease 欄位
        for name in (entry.get("disease") or "").split("|"):
            name = name.strip()
            if name and name != "-":
                counted[name] += 1
    return counted


def _split(entry: Dict, key: str, n: int) -> List[str]:
    parts = [_blank(p) for p in str(entry.get(key) or "").split(" | ")]
    return parts if len(parts) == n else [""] * n


def _publications(entry: Dict) -> List[Dict[str, str]]:
    """PubMind 的多值欄位以 ' | ' 串接、每一段對應一筆文獻紀錄（實測 840 筆全部對齊），
    拆開後每列就是「一篇文獻對這顆變異的描述」，可以跟 HGMD / PubTator 一樣逐列列出 PMID。"""
    n = int(entry.get("Num_of_record_used") or 0)
    refs = _split(entry, "formatted_reference", n)
    disease = _split(entry, "disease", n)
    mondo = _split(entry, "MONDO_name_09", n)
    pathogenicity = _split(entry, "pathogenicity", n)
    reasoning = _split(entry, "LLM_reasoning", n)
    phenotype = _split(entry, "phenotype", n)
    rows = []
    for j in range(n):
        m = RE_REFERENCE.match(refs[j].strip())
        pmc, pmid, section = (m.groups() if m else (None, None, None))
        rows.append({
            "pmid": pmid or "",
            "pmcid": pmc or "",
            "section": _blank(section),
            "disease": disease[j] or mondo[j],
            "pathogenicity": pathogenicity[j],
            "evidence": reasoning[j] or phenotype[j],
        })
    return rows


def _record(entry: Dict, level: str) -> Dict[str, Any]:
    dna = _blank(entry.get("dna_change"))
    aa = _blank(entry.get("aa_change"))
    try:
        score: Optional[float] = float(entry.get("pathogenicity_score"))
    except (TypeError, ValueError):
        score = None
    return {
        "pvid": _blank(entry.get("PVID")),
        "label": " / ".join(filter(None, [dna and f"c.{dna}", aa and f"p.{aa}"])),
        "level": level,
        "sig": _blank(entry.get("pathogenicity_sum")),
        "score": score,
        # 公開下載版（/download_json）沒有 Num_of_paper_used，改數 PMCID_PMID_counted 的文獻數
        "papers": int(entry.get("Num_of_paper_used") or 0)
        or len(_counted(entry.get("PMCID_PMID_counted") or "")),
        "records": int(entry.get("Num_of_record_used") or 0),
        "diseases": _diseases(entry).most_common(),
        "phenotypes": _counted(entry.get("HPO_term_counted") or "").most_common(),
        "pmids": list(dict.fromkeys(RE_PMID.findall(entry.get("PMCID_PMID_counted") or ""))),
        "confidence": _blank(entry.get("confidence_criteria")).replace(";", " · "),
        "publications": _publications(entry),
    }


class PubMindLookup:
    """PubMind TTN 變異查表器，key 為報告格式的 pos-ref-alt（轉錄方向 allele）。"""

    def __init__(self, path: Path = DEFAULT_PUBMIND_JSON) -> None:
        self.path = Path(path)
        self.by_variant: Dict[str, List[Dict[str, Any]]] = {}
        self.loaded = False
        if not self.path.exists():
            logger.warning("PubMind data not found: %s", self.path)
            return
        entries = json.loads(self.path.read_text(encoding="utf-8"))
        for entry in entries:
            raw = entry.get("parsed_variants")
            if not raw or raw == "-":
                continue
            try:
                variants = ast.literal_eval(raw)
            except (ValueError, SyntaxError):
                continue
            seen = set()
            for v in variants:
                ref, alt = str(v.get("ref", "")), str(v.get("alt", ""))
                if str(v.get("chr")) != "2" or len(ref) != len(alt) or "-" in ref + alt:
                    continue
                diff = [i for i in range(len(ref)) if ref[i] != alt[i]]
                if len(diff) != 1:
                    continue
                i = diff[0]
                key = f"{int(v['pos']) + i}-{ref[i].translate(COMPLEMENT)}-{alt[i].translate(COMPLEMENT)}"
                if key in seen:
                    continue
                seen.add(key)
                level = "dna" if len(ref) == 1 else "protein"
                self.by_variant.setdefault(key, []).append(_record(entry, level))
        self.loaded = True
        logger.info("loaded PubMind: %d entries -> %d SNV keys", len(entries), len(self.by_variant))

    def lookup(self, pos: int, ref: str, alt: str) -> List[Dict[str, Any]]:
        return list(self.by_variant.get(f"{int(pos)}-{ref.upper()}-{alt.upper()}", []))


def disease_counts(records: List[Dict[str, Any]]) -> List[tuple]:
    counter: Counter = Counter()
    for rec in records:
        for name, n in rec.get("diseases") or []:
            counter[name] += n
    return counter.most_common()
