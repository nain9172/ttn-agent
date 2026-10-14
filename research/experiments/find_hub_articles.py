#!/usr/bin/env python3
"""從 outputs/ 的報告統計「一篇文獻提及了幾顆變異」，挑出 hub 文獻與它們覆蓋的變異。

為什麼要這樣挑評估集：第一輪的評估集是隨機分層抽樣，alias gate 只放行 14% 的文獻，
六個機制設定裡有四個結果完全相同——因為真正進到 LLM 的文獻太少，機制沒有作用空間。
提及大量變異的文獻（通常是帶補充表的大型 cohort 研究）正好是機制設計要對付的情境：
同一張表裡有幾十顆位置相近的變異，最容易張冠李戴。用這些文獻覆蓋的變異當評估集，
alias gate、evidence grounding 這些機制才有東西可擋。

判定「提及」直接沿用報告裡的 badge（html_report._article_mentions_variant 的結果），
與 Tissue Involvement Statistics 的計數同源。

用法（conda 環境 gm4）：
    python research/experiments/find_hub_articles.py                    # 只看統計
    python research/experiments/find_hub_articles.py --top 5 --out research/experiments/llm_eval_set_hub.csv
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = PROJECT_ROOT / "outputs"
GT_XLSX = PROJECT_ROOT / "HGMD_TTN_with_clinvar_pubcount.xlsx"

# html_report.py 產生的文獻條目：<a href=".../pubmed/PMID/">標題</a> …
# <span class="mention-badge mention-true">Mention input variant: True</span>
ARTICLE_RE = re.compile(
    r'article-title">.*?pubmed\.ncbi\.nlm\.nih\.gov/(\d+)'   # PMID
    r'.*?mention-badge mention-(true|false)',                 # 有沒有提到輸入變異
    re.DOTALL,
)
# 報告檔名：<variant>_<model tag>.html；no-llm 的報告沒有跑 LLM，mention 一律 False
REPORT_RE = re.compile(r"^(\w+-\d+-[ACGT]+-[ACGT]+)_(.+)\.html$", re.IGNORECASE)


def scan_reports(skip_no_llm: bool = True) -> pd.DataFrame:
    """回傳 (variant, pmid, mention) 三欄的表，一列一個 (變異, 文獻) 配對。"""
    rows = []
    files = sorted(OUTPUT_DIR.glob("*.html"))
    for i, path in enumerate(files, 1):
        m = REPORT_RE.match(path.name)
        if not m:
            continue
        variant, tag = m.group(1), m.group(2)
        if skip_no_llm and tag == "no-llm":
            continue
        try:
            html = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for pmid, mention in ARTICLE_RE.findall(html):
            rows.append((variant, pmid, mention == "true", tag))
        if i % 500 == 0:
            print(f"  掃描 {i}/{len(files)}…", flush=True)

    df = pd.DataFrame(rows, columns=["variant", "pmid", "mention", "model_tag"])
    # 同一顆變異可能有多份報告（不同模型 tag）；同一篇文獻只要任一份報告判定提及就算提及
    return (df.groupby(["variant", "pmid"], as_index=False)["mention"].max())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=5, help="取前幾篇 hub 文獻（預設 5）")
    ap.add_argument("--max-variants", type=int, default=60,
                    help="評估集上限；超過時取被最多 hub 文獻提及的變異（預設 60）")
    ap.add_argument("--out", type=Path, default=None, help="評估集輸出 CSV")
    ap.add_argument("--pairs-out", type=Path,
                    default=PROJECT_ROOT / "research" / "experiments" / "report_article_variant_pairs.csv",
                    help="所有 (變異, 文獻) 配對的輸出，之後不必重掃")
    args = ap.parse_args()

    if args.pairs_out.exists():
        print(f"[info] 直接讀既有的 {args.pairs_out}（要重掃請先刪掉）")
        pairs = pd.read_csv(args.pairs_out)
    else:
        print(f"[info] 掃描 {OUTPUT_DIR} 的報告…")
        pairs = scan_reports()
        pairs.to_csv(args.pairs_out, index=False)
        print(f"[info] 寫出 {args.pairs_out}")

    print(f"\n共 {len(pairs)} 個 (變異, 文獻) 配對，"
          f"{pairs.variant.nunique()} 顆變異，{pairs.pmid.nunique()} 篇文獻")

    mentioned = pairs[pairs["mention"]]
    counts = (mentioned.groupby("pmid")
              .agg(n_variants=("variant", "nunique"))
              .sort_values("n_variants", ascending=False))
    print(f"\n提及最多變異的文獻（共 {len(counts)} 篇曾提及至少一顆）：")
    print(counts.head(15).to_string())

    # ground truth：只有 HGMD 有表型標籤的變異才進得了評估集
    gt = pd.read_excel(GT_XLSX)
    gt["gt_label"] = gt["phenotype"].map({"cardiac": "cardiac", "skeletal": "skeletal",
                                          "both": "both", "unsure": "others"})
    gt = gt.dropna(subset=["input", "gt_label"])

    hubs = counts.head(args.top).index.tolist()
    hub_variants = mentioned[mentioned["pmid"].isin(hubs)]
    per_variant = (hub_variants.groupby("variant")
                   .agg(n_hub_articles=("pmid", "nunique"),
                        hub_pmids=("pmid", lambda s: ";".join(sorted({str(v) for v in s})))))

    eval_set = (gt.merge(per_variant, left_on="input", right_index=True, how="inner")
                .sort_values(["n_hub_articles", "input"], ascending=[False, True]))

    # 同一顆變異在 HGMD 可能有多筆記錄。表型一致的去重即可；互相衝突的（例如同時
    # 登記成 skeletal 與 cardiac）真值無法判定，整顆排除，不能讓它被計分兩次。
    labels_per_variant = eval_set.groupby("input")["gt_label"].nunique()
    conflicted = labels_per_variant[labels_per_variant > 1].index.tolist()
    if conflicted:
        print(f"\n[warn] {len(conflicted)} 顆變異的 HGMD 表型互相衝突，已排除：")
        for vid in conflicted:
            rows = eval_set[eval_set["input"] == vid]
            print(f"  {vid}: " + " / ".join(f"{r.ID}={r.PHEN}" for r in rows.itertuples()))
        eval_set = eval_set[~eval_set["input"].isin(conflicted)]
    eval_set = eval_set.drop_duplicates("input")
    print(f"\n前 {args.top} 篇 hub 文獻共提及 {len(per_variant)} 顆變異，"
          f"其中 {len(eval_set)} 顆在 HGMD 有表型標籤")
    if eval_set.empty:
        print("[warn] 沒有可用的變異，試著加大 --top")
        return 1
    print(eval_set["gt_label"].value_counts().to_string())

    if args.max_variants and len(eval_set) > args.max_variants:
        eval_set = eval_set.head(args.max_variants)
        print(f"\n取前 {args.max_variants} 顆（被最多 hub 文獻提及的）")
        print(eval_set["gt_label"].value_counts().to_string())

    cols = ["input", "#CHROM", "POS", "REF", "ALT", "ID", "DNA", "PROT", "PHEN",
            "phenotype", "gt_label", "n_hub_articles", "hub_pmids", "ClinVar_Variation_ID"]
    eval_set = eval_set[[c for c in cols if c in eval_set.columns]].reset_index(drop=True)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        eval_set.to_csv(args.out, index=False)
        print(f"\n評估集 {len(eval_set)} 顆 -> {args.out}")
    else:
        print("\n（沒有給 --out，只印統計）")
        print(eval_set.head(10).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
