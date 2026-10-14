#!/usr/bin/env python3
"""把缺少 Evo2 Prediction 區塊的報告補上分數。

為什麼會缺：這些報告是用 `--skip-evo2` 跑出來的（批次為了省時間跳過序列層級
計分），所以報告裡根本沒有那一段，不是算出來失敗。

分數來源，依可信度排序：
  1. 同一顆變異的另一份報告  同一套流程算出來的，數值與渲染完全一致
  2. outputs/evo2_scores_*.csv（score_evo2_for_variants.py；新的檔案優先）
  2b. outputs/evo2_scoring_checkpoint.csv      7,179 顆，Evo2-7B / 8,192 bp（早期在 gm3 算的）
     （evo2_ttn_validation_results.csv 是它的子集，數值相同）
  3. experiments/evo2_window_scores.csv (window=8192)  1,800 顆，實驗一產生

注意：2 和 3 是另一次執行的結果。Evo2 的 Hyena FFT 卷積不是位元級決定性的，
同一顆變異重跑的 delta 大約會差 1e-4；門檻是 -7.14e-4，所以**貼近門檻的變異，
補上的標籤有可能與這份報告自己跑會得到的結果不同**。每個補上的區塊都會標註
來源，並記錄在 manifest 裡，可以還原。

用法（conda 環境 gm4）：
    python backfill_evo2_scores.py --dry-run          # 只統計，不改檔
    python backfill_evo2_scores.py --limit 3          # 先改 3 份看看
    python backfill_evo2_scores.py                    # 全部
    python backfill_evo2_scores.py --revert           # 依 manifest 還原
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

ROOT = Path(__file__).parent.resolve()
sys.path.insert(0, str(ROOT))

import pandas as pd

from config import PATHOGENIC_THRESHOLD

OUTPUT_DIR = ROOT / "outputs"
CHECKPOINT_CSV = OUTPUT_DIR / "evo2_scoring_checkpoint.csv"
# score_evo2_for_variants.py 補算出來的分數（欄位與 checkpoint 相同，多一欄 vid）
EXTRA_CSVS = sorted(OUTPUT_DIR.glob("evo2_scores_*.csv"))
WINDOW_CSV = ROOT / "research" / "experiments" / "evo2_window_scores.csv"
MANIFEST = ROOT / "outputs" / "evo2_backfill_manifest.json"

EVO2_HEADING = "<h2>Evo2 Prediction</h2>"
REPORT_RE = re.compile(r"^(\w+-\d+-[ACGT]+-[ACGT]+)_(.+)\.html$", re.IGNORECASE)
# 從既有報告讀回三個分數（順序固定：delta / ref / var）
SCORE_RE = re.compile(
    r'<h2>Evo2 Prediction</h2>.*?'
    r'prediction-result">(\w+)</div>.*?'
    r'score-value">(-?[\d.]+)</div><div>Delta Score.*?'
    r'score-value">(-?[\d.]+)</div><div>Ref Score.*?'
    r'score-value">(-?[\d.]+)</div><div>Var Score',
    re.DOTALL,
)
# 插入點：Evo2 區塊原本就排在這幾段之前，依序嘗試
ANCHORS = [
    '<section class="section"><h2>Protein Domain Localization</h2>',
    '<section class="section"><h2>Literature Review',
    '<section class="section"><h2>Mechanistic Inference',
    '<section class="section"><h2>Curated Database Annotations',
]


def evo2_section(delta: float, ref: float, var: float, prediction: str) -> str:
    """與 html_report._generate_evo2_section 完全相同的標記，不加任何額外註解。

    prediction 一律沿用來源的判定（既有報告的標題、CSV 的
    pred_by_default_threshold），不在這裡重算門檻。
    """
    cls = "prediction-pathogenic" if prediction == "pathogenic" else "prediction-benign"
    return (
        f'<section class="section"><h2>Evo2 Prediction</h2>\n'
        f'            <div class="prediction-box {cls}">'
        f'<div class="prediction-result">{prediction.upper()}</div>\n'
        f'            <div class="score-display">\n'
        f'                <div><div class="score-value">{delta:.4f}</div><div>Delta Score</div></div>\n'
        f'                <div><div class="score-value">{ref:.4f}</div><div>Ref Score</div></div>\n'
        f'                <div><div class="score-value">{var:.4f}</div><div>Var Score</div></div>\n'
        f'            </div></div></section>'
    )


def build_score_index() -> Dict[str, Tuple[float, float, float, str, str]]:
    """variant_id -> (delta, ref, var, prediction, 來源標籤)，依可信度由高到低填入。"""
    scores: Dict[str, Tuple[float, float, float, str, str]] = {}

    # 1) 既有報告（同一套流程）
    for path in OUTPUT_DIR.glob("*.html"):
        m = REPORT_RE.match(path.name)
        if not m:
            continue
        vid = m.group(1)
        if vid in scores:
            continue
        head = path.read_text(encoding="utf-8", errors="replace")
        hit = SCORE_RE.search(head)
        if hit:
            scores[vid] = (float(hit.group(2)), float(hit.group(3)), float(hit.group(4)),
                           hit.group(1).lower(), "report")

    # 2) score_evo2_for_variants.py 算出的分數（setup_and_run.sh 這次跑的），新的優先；
    #    排在舊 checkpoint 之前，確保用的是目前環境算出來的分數
    for extra in sorted(EXTRA_CSVS, key=lambda f: f.stat().st_mtime, reverse=True):
        e = pd.read_csv(extra)
        e = e[e["vid"] != "vid"]           # 逐筆 append 時可能混入重複表頭
        for row in e.itertuples():
            if row.vid not in scores and pd.notna(row.delta_score):
                scores[row.vid] = (float(row.delta_score), float(row.ref_score),
                                   float(row.var_score),
                                   str(row.pred_by_default_threshold).strip().lower(),
                                   extra.stem)

    # 2b) 早期全量計分的 checkpoint（7,179 顆，涵蓋 ClinVar 之外的 HGMD 變異）
    if CHECKPOINT_CSV.exists():
        v = pd.read_csv(CHECKPOINT_CSV)
        v["vid"] = (v["chr"].astype(str) + "-" + v["position"].astype(str)
                    + "-" + v["ref"].astype(str) + "-" + v["alt"].astype(str))
        for row in v.itertuples():
            if row.vid not in scores and pd.notna(row.delta_score):
                scores[row.vid] = (float(row.delta_score), float(row.ref_score),
                                   float(row.var_score),
                                   str(row.pred_by_default_threshold).strip().lower(),
                                   "checkpoint_csv")

    # 3) 實驗一的 8,192 bp 結果
    if WINDOW_CSV.exists():
        w = pd.read_csv(WINDOW_CSV)
        w = w[w["window"] == 8192].drop_duplicates("variant_id")
        for row in w.itertuples():
            if row.variant_id not in scores:
                delta = float(row.delta_mean_window)
                scores[row.variant_id] = (
                    delta, float(row.ref_score), float(row.var_score),
                    "pathogenic" if delta < PATHOGENIC_THRESHOLD else "benign", "window_csv")

    return scores


def revert() -> int:
    if not MANIFEST.exists():
        print(f"[error] 找不到 {MANIFEST}", file=sys.stderr)
        return 1
    entries = json.loads(MANIFEST.read_text(encoding="utf-8"))
    reverted = 0
    for entry in entries:
        path = Path(entry["file"])
        if not path.exists() or not entry.get("block"):
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        new = text.replace(entry["block"], "", 1)
        if new != text:
            path.write_text(new, encoding="utf-8")
            reverted += 1
    MANIFEST.unlink()
    print(f"已還原 {reverted} 份報告，manifest 已刪除")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="只統計，不改檔")
    ap.add_argument("--limit", type=int, default=None, help="只處理前 N 份")
    ap.add_argument("--revert", action="store_true", help="依 manifest 移除補上的區塊")
    args = ap.parse_args()

    if args.revert:
        return revert()

    print("[1/3] 建立分數索引…", flush=True)
    scores = build_score_index()
    by_source: Dict[str, int] = {}
    for *_, src in scores.values():
        by_source[src] = by_source.get(src, 0) + 1
    print(f"      {len(scores)} 顆變異有分數 {by_source}")

    print("[2/3] 掃描缺少 Evo2 區塊的報告…", flush=True)
    todo = []
    no_score = []
    for path in sorted(OUTPUT_DIR.glob("*.html")):
        m = REPORT_RE.match(path.name)
        if not m:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if EVO2_HEADING in text:
            continue
        vid = m.group(1)
        if vid not in scores:
            no_score.append(path.name)
            continue
        todo.append((path, vid, text))
    print(f"      缺分數 {len(todo) + len(no_score)} 份；可補 {len(todo)}；無來源 {len(no_score)}")
    if no_score:
        print(f"      無來源的: {no_score[:10]}")

    if args.limit:
        todo = todo[: args.limit]
    if args.dry_run:
        print("[dry-run] 不改檔")
        return 0

    print(f"[3/3] 寫入 {len(todo)} 份…", flush=True)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8")) if MANIFEST.exists() else []
    written = skipped = 0
    for i, (path, vid, text) in enumerate(todo, 1):
        delta, ref, var, prediction, src = scores[vid]
        anchor = next((a for a in ANCHORS if a in text), None)
        if anchor is None:
            skipped += 1
            continue
        section = evo2_section(delta, ref, var, prediction)
        path.write_text(text.replace(anchor, section + anchor, 1), encoding="utf-8")
        manifest.append({"file": str(path), "variant": vid, "source": src,
                         "delta": round(delta, 6), "prediction": prediction,
                         "block": section})
        written += 1
        if i % 500 == 0:
            print(f"      {i}/{len(todo)}", flush=True)

    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n完成：補上 {written} 份；找不到插入點而略過 {skipped} 份")
    print(f"manifest: {MANIFEST}（--revert 可還原）")
    if no_score:
        print(f"\n仍然缺分數的 {len(no_score)} 份：{no_score}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
