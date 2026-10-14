#!/usr/bin/env python3
"""對指定的變異跑 Evo2 計分，輸出成與 evo2_scoring_checkpoint.csv 同格式的 CSV。

用途：某些報告是用 `--skip-evo2` 跑出來的，缺 Evo2 Prediction 區塊。這支腳本把
分數補算出來，再用 backfill_evo2_scores.py 寫回 HTML。

計分方式與正式流程完全相同（utils/evo2_predictor.Evo2Predictor）：
以變異為中心取 8,192 bp 窗口、ref 與 var 各算一次整窗平均 log-likelihood，
delta = var - ref；致病性用 config.PATHOGENIC_THRESHOLD 判定。

用法（conda 環境 gm4）：
    python score_evo2_for_variants.py --input "TTN variants test.xlsx" \
        --out outputs/evo2_scores_ttn_test.csv
    python score_evo2_for_variants.py --input data/foo.csv --only-missing-report
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import time
from pathlib import Path
from typing import List

ROOT = Path(__file__).parent.resolve()
sys.path.insert(0, str(ROOT))
os.environ.setdefault("TRITON_PTXAS_PATH", "/usr/local/cuda/bin/ptxas")

import pandas as pd

from config import PATHOGENIC_THRESHOLD

VARIANT_RE = re.compile(r"^(\w+)-(\d+)-([ACGT]+)-([ACGT]+)$", re.IGNORECASE)
EVO2_HEADING = "<h2>Evo2 Prediction</h2>"


def read_inputs(path: Path) -> List[str]:
    df = pd.read_excel(path) if path.suffix.lower() in (".xlsx", ".xls") else pd.read_csv(path)
    if "input" not in df.columns:
        raise SystemExit(f"{path} 沒有 input 欄")
    out, seen = [], set()
    for raw in df["input"].dropna():
        v = str(raw).strip()
        if VARIANT_RE.match(v) and v not in seen:
            seen.add(v)
            out.append(v)
    return out


def report_has_evo2(variant: str) -> bool:
    files = glob.glob(str(ROOT / "outputs" / f"{variant}_*.html"))
    return any(EVO2_HEADING in Path(f).read_text(encoding="utf-8", errors="ignore") for f in files)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, required=True, help="含 input 欄的 xlsx / csv")
    ap.add_argument("--out", type=Path, required=True, help="輸出的分數 CSV")
    ap.add_argument("--only-missing-report", action="store_true", default=True,
                    help="只算報告裡還沒有 Evo2 區塊的變異（預設開啟）")
    ap.add_argument("--all", dest="only_missing_report", action="store_false",
                    help="不管報告狀態，全部重算")
    args = ap.parse_args()

    variants = read_inputs(args.input)
    print(f"[info] {args.input} 共 {len(variants)} 顆合法變異")
    if args.only_missing_report:
        variants = [v for v in variants if not report_has_evo2(v)]
        print(f"[info] 其中報告缺 Evo2 區塊的 {len(variants)} 顆")

    done = set()
    if args.out.exists():
        done = set(pd.read_csv(args.out)["vid"])
        print(f"[info] 續跑：{len(done)} 顆已計分")
    todo = [v for v in variants if v not in done]
    if not todo:
        print("[info] 沒有待計分的變異")
        return 0

    from utils.evo2_predictor import Evo2Predictor

    print(f"[info] 載入 Evo2（{len(todo)} 顆待計分）…", flush=True)
    predictor = Evo2Predictor()
    predictor._load_model()
    predictor._load_reference_sequence()

    rows, failed = [], []
    for i, vid in enumerate(todo, 1):
        chrom, pos, ref, alt = VARIANT_RE.match(vid).groups()
        t0 = time.time()
        result = predictor.predict({"variant_id": vid, "chrom": chrom, "pos": int(pos),
                                    "ref": ref.upper(), "alt": alt.upper()})
        if not result.get("success"):
            failed.append((vid, result.get("error", "unknown")))
            print(f"[{i}/{len(todo)}] {vid} 失敗：{result.get('error')}", flush=True)
            continue
        delta = float(result["delta_score"])
        rows.append({
            "vid": vid, "chr": chrom, "position": int(pos), "ref": ref.upper(), "alt": alt.upper(),
            "ref_score": float(result["ref_score"]), "var_score": float(result["var_score"]),
            "delta_score": delta,
            "pred_by_default_threshold": "pathogenic" if delta < PATHOGENIC_THRESHOLD else "benign",
        })
        print(f"[{i}/{len(todo)}] {vid} delta={delta:+.5f} "
              f"-> {rows[-1]['pred_by_default_threshold']} ({time.time() - t0:.1f}s)", flush=True)
        # 每筆都存檔，中斷可續跑
        pd.DataFrame(rows).to_csv(args.out, mode="a", index=False,
                                  header=not args.out.exists())
        rows = []

    print(f"\n完成，輸出：{args.out}")
    if failed:
        print(f"失敗 {len(failed)} 顆：{failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
