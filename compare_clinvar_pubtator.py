#!/usr/bin/env python3
"""
以 ClinVar Variation ID 比對 HGMD 變異表與 PubTator 變異表，統計每個 HGMD 變異
在 PubTator 文獻中被標註的組織分佈。

流程：
  1. 讀 pubtator_ttn_variants_all_with_clinvar.csv，依 ClinVar Variation ID 分組，
     統計每個 ID 的 tissue_affected 次數（同一變異可能被多篇文獻記錄多次）：
        Cardiac  → pub_cardiac
        Skeletal → pub_skeletal
        Both     → pub_both
     （Not specified / Others 不計入這三欄）
  2. 以 HGMD_TTN_with_clinvar.xlsx 為基礎，依其 ClinVar Variation ID 查表，
     新增 pub_cardiac / pub_skeletal / pub_both 三欄後另存。

ClinVar ID 在兩邊讀進來可能是 int / float（帶 .0）/ str，統一正規化成整數字串再比對。

用法：
    conda run -n gm4 python compare_clinvar_pubtator.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("compare_clinvar_pubtator")

PROJECT_ROOT = Path(__file__).parent
DEFAULT_HGMD = PROJECT_ROOT / "HGMD_TTN_with_clinvar.xlsx"
DEFAULT_PUBTATOR = PROJECT_ROOT / "pubtator_ttn_variants_all_with_clinvar.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "HGMD_TTN_with_clinvar_pubcount.xlsx"

CLINVAR_COL = "ClinVar_Variation_ID"
TISSUE_COL = "tissue_affected"

# tissue_affected 值 → 輸出欄位
TISSUE_TO_COLUMN = {
    "cardiac": "pub_cardiac",
    "skeletal": "pub_skeletal",
    "both": "pub_both",
}
OUTPUT_COLUMNS = ["pub_cardiac", "pub_skeletal", "pub_both"]
OVERLAP_COLUMN = "overlap"


def normalize_clinvar_id(value) -> Optional[str]:
    """把 ClinVar ID 統一成整數字串（去掉浮點 .0、空白），無效值回 None。"""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return None
    try:
        return str(int(float(s)))
    except (ValueError, TypeError):
        return s


def build_tissue_counts(pub: pd.DataFrame) -> Dict[str, Counter]:
    """依 ClinVar ID 分組，統計各組織標註次數。"""
    counts: Dict[str, Counter] = defaultdict(Counter)
    for _, row in pub.iterrows():
        cid = normalize_clinvar_id(row.get(CLINVAR_COL))
        if not cid:
            continue
        tissue = str(row.get(TISSUE_COL, "")).strip().lower()
        col = TISSUE_TO_COLUMN.get(tissue)
        if col:
            counts[cid][col] += 1
    return counts


def main() -> None:
    ap = argparse.ArgumentParser(
        description="以 ClinVar ID 比對 HGMD 與 PubTator，統計組織分佈"
    )
    ap.add_argument("--hgmd", default=str(DEFAULT_HGMD), help="HGMD（含 ClinVar ID）Excel")
    ap.add_argument(
        "--pubtator", default=str(DEFAULT_PUBTATOR), help="PubTator（含 ClinVar ID）CSV"
    )
    ap.add_argument("--output", default=str(DEFAULT_OUTPUT), help="輸出 Excel 路徑")
    args = ap.parse_args()

    hgmd_path = Path(args.hgmd)
    pub_path = Path(args.pubtator)
    out_path = Path(args.output)
    for p in (hgmd_path, pub_path):
        if not p.exists():
            logger.error(f"找不到輸入檔：{p}")
            sys.exit(1)

    pub = pd.read_csv(pub_path)
    hgmd = pd.read_excel(hgmd_path)
    logger.info(f"PubTator：{len(pub)} 列；HGMD：{len(hgmd)} 列")

    for name, df in (("PubTator", pub), ("HGMD", hgmd)):
        if CLINVAR_COL not in df.columns:
            logger.error(f"{name} 缺少欄位 {CLINVAR_COL}")
            sys.exit(1)
    if TISSUE_COL not in pub.columns:
        logger.error(f"PubTator 缺少欄位 {TISSUE_COL}")
        sys.exit(1)

    counts = build_tissue_counts(pub)
    logger.info(f"PubTator 中有 {len(counts)} 個不同 ClinVar ID 帶有組織標註")

    # overlap 以「ID 是否出現在 PubTator」為準（含 tissue 為 Not specified/Others 者），
    # 無法用 ID 比對（HGMD 無 ClinVar ID，或該 ID 不在 PubTator）計為 False。
    pub_ids = {
        cid for cid in (normalize_clinvar_id(v) for v in pub[CLINVAR_COL]) if cid
    }
    logger.info(f"PubTator 中共有 {len(pub_ids)} 個不同 ClinVar ID")

    col_data = {c: [] for c in OUTPUT_COLUMNS}
    overlap_data = []
    n_matched = 0
    for _, row in hgmd.iterrows():
        cid = normalize_clinvar_id(row.get(CLINVAR_COL))
        is_overlap = bool(cid) and cid in pub_ids
        overlap_data.append(is_overlap)
        if is_overlap:
            n_matched += 1
        c = counts.get(cid) if cid else None
        for col in OUTPUT_COLUMNS:
            col_data[col].append(int(c[col]) if c else 0)

    for col in OUTPUT_COLUMNS:
        hgmd[col] = col_data[col]
    hgmd[OVERLAP_COLUMN] = overlap_data

    hgmd.to_excel(out_path, index=False)

    totals = {col: int(sum(col_data[col])) for col in OUTPUT_COLUMNS}
    logger.info("=" * 60)
    logger.info(
        f"HGMD 變異中有 {n_matched}/{len(hgmd)} 筆 overlap=True"
        f"（其餘 {len(hgmd) - n_matched} 筆無法用 ID 比對，overlap=False）"
    )
    logger.info(
        f"累計標註次數 → pub_cardiac={totals['pub_cardiac']}, "
        f"pub_skeletal={totals['pub_skeletal']}, pub_both={totals['pub_both']}"
    )
    logger.info(f"輸出已寫入：{out_path}")


if __name__ == "__main__":
    main()
