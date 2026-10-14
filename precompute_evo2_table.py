#!/usr/bin/env python3
"""預先計算整條 TTN 上所有 SNV 的 Evo2 分數，給 site 的「未收錄變異」查表用。

計分方式與 utils/evo2_predictor.Evo2Predictor.predict 完全相同：以變異為中心取
EVO2_WINDOW_SIZE（8,192 bp）窗口，整窗平均 log-likelihood，delta = var - ref。
差別只在同一個位置的 ref 與三個 alt 放進同一個 forward（batch 4），ref 只算一次；
實測與逐顆呼叫 score_sequences 的結果逐位元相同，速度約快 1.8 倍。

輸出（--out-dir，預設 outputs/evo2_table/）：
  scores.f64   float64 memmap，shape (L, 4)；[i, b] = 序列索引 i 換成鹼基 "ACGT"[b]
               後的窗口平均 log-likelihood（ref 那一格就是 ref 分數）。
               delta = scores[i, alt] - scores[i, ref]。
  done.u8      uint8 memmap，shape (L,)；1 = 該位置已算完。
  meta.json    座標慣例與參數。
  progress.log 進度。
  evo2_table.js  給 site 查表用的精簡版（--export；每次執行結束也會自動更新）：
               float16 delta，shape (L, 4)，ref 那一格為 0、未算的位置為 NaN，
               base64 包在 window.EVO2_TABLE 裡，從 file:// 開啟也能用 <script> 載入。

座標：序列索引 i = TTN_SEQUENCE_START - pos（data/sequence.fasta 是負股、已反向互補），
鹼基是轉錄方向，與報告 variant ID 的 ref/alt 一致。

順序：先算 7 條轉錄本的外顯子 + 剪接區 ±10 bp，再算內含子；可隨時中斷續跑。
停止方式：Ctrl-C / kill（SIGTERM），或建立 <out-dir>/STOP 檔，會在目前 batch 結束後存檔離開。

用法（conda 環境 gm4，GPU 需空出來）：一般透過 init_evo2_table.sh 啟動；直接呼叫：
    python precompute_evo2_table.py            # 計算（可續跑）
    python precompute_evo2_table.py --status   # 進度
    python precompute_evo2_table.py --export   # 只重新匯出 evo2_table.js
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent.resolve()
sys.path.insert(0, str(ROOT))
os.environ.setdefault("TRITON_PTXAS_PATH", "/usr/local/cuda/bin/ptxas")

import numpy as np

import config

BASES = "ACGT"
SPLICE_FLANK = 10


def hot_mask(n: int) -> np.ndarray:
    """7 條轉錄本的外顯子 + 剪接區 ±SPLICE_FLANK，以序列索引表示。

    沒有轉錄本模型時全部視為同等優先（照序列順序算），不影響結果，只影響先後。
    """
    hot = np.zeros(n, dtype=bool)
    if not config.TTN_TRANSCRIPT_MODEL_PATH.exists():
        return hot
    model = json.loads(config.TTN_TRANSCRIPT_MODEL_PATH.read_text())
    for tx in model["transcripts"].values():
        for e in tx["exons"]:
            lo, hi = int(e["g_start"]) - SPLICE_FLANK, int(e["g_end"]) + SPLICE_FLANK
            # genomic → 序列索引（負股，pos 越小索引越大）
            a = max(0, config.TTN_SEQUENCE_START - hi)
            b = min(n - 1, config.TTN_SEQUENCE_START - lo)
            if a <= b:
                hot[a:b + 1] = True
    return hot


def priority_order(n: int) -> np.ndarray:
    """外顯子 + 剪接區排前面，內含子接在後面。"""
    hot = hot_mask(n)
    idx = np.arange(n)
    return np.concatenate([idx[hot], idx[~hot]])


def open_arrays(out: Path, n: int):
    out.mkdir(parents=True, exist_ok=True)
    scores_path, done_path = out / "scores.f64", out / "done.u8"
    fresh = not scores_path.exists()
    scores = np.memmap(scores_path, dtype=np.float64, mode="w+" if fresh else "r+", shape=(n, 4))
    done = np.memmap(done_path, dtype=np.uint8, mode="w+" if fresh else "r+", shape=(n,))
    if fresh:
        scores[:] = np.nan
        done[:] = 0
    return scores, done


def export_table(out: Path, n: int, dest: Path) -> None:
    """把 scores 轉成 site 用的 float16 delta 表。"""
    import base64

    if not (out / "done.u8").exists():
        print("還沒有任何計算結果，無法匯出")
        return
    scores = np.memmap(out / "scores.f64", dtype=np.float64, mode="r", shape=(n, 4))
    done = np.memmap(out / "done.u8", dtype=np.uint8, mode="r", shape=(n,)).astype(bool)
    seq = _load_sequence()
    ref_idx = np.array([BASES.find(b) for b in seq])
    ref_scores = scores[np.arange(n), np.clip(ref_idx, 0, 3)]
    delta = (scores - ref_scores[:, None]).astype(np.float16)
    delta[~done] = np.nan
    meta = json.loads((out / "meta.json").read_text(encoding="utf-8"))
    meta.update({
        "encoding": "base64 little-endian float16, row-major (length, 4); "
                    "delta of base 'ACGT'[b] vs reference; NaN = not computed",
        "n_done": int(done.sum()),
    })
    payload = base64.b64encode(delta.astype("<f2").tobytes()).decode("ascii")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".tmp")
    tmp.write_text(
        "window.EVO2_TABLE = " + json.dumps({**meta, "data": payload}, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    tmp.replace(dest)
    print(f"已匯出 {dest}（{int(done.sum()):,} / {n:,} 個位置，{dest.stat().st_size / 1e6:.1f} MB）")


def _load_sequence() -> str:
    lines = config.REFERENCE_GENOME_PATH.read_text().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">")).upper()


def _env_info() -> dict:
    info = {}
    try:
        from importlib.metadata import version

        for pkg in ("evo2", "torch", "transformer_engine", "vtx"):
            try:
                info[pkg] = version(pkg)
            except Exception:
                pass
    except ImportError:
        pass
    return info


def status(out: Path, n: int) -> None:
    if not (out / "done.u8").exists():
        print("尚未開始")
        return
    done = np.memmap(out / "done.u8", dtype=np.uint8, mode="r", shape=(n,)).astype(bool)
    hot = hot_mask(n)
    print(f"完成 {int(done.sum()):,} / {n:,} 個位置（{done.mean():.1%}）")
    print(f"  外顯子 + 剪接區 {int((done & hot).sum()):,} / {int(hot.sum()):,}")
    print(f"  內含子          {int((done & ~hot).sum()):,} / {int((~hot).sum()):,}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=config.OUTPUT_DIR / "evo2_table")
    ap.add_argument("--positions-per-batch", type=int, default=2,
                    help="一次 forward 放幾個位置（每個位置 4 條序列）；實測 2 比 1 快約 3%%")
    ap.add_argument("--limit", type=int, default=None, help="本次最多算幾個位置（試跑用）")
    ap.add_argument("--status", action="store_true", help="只顯示進度")
    ap.add_argument("--export", action="store_true", help="只把目前結果匯出成 evo2_table.js")
    ap.add_argument("--export-path", type=Path, default=None,
                    help="evo2_table.js 輸出位置（預設 <out-dir>/evo2_table.js）")
    ap.add_argument("--validate", type=int, default=0, metavar="N",
                    help="從 evo2_scoring_checkpoint.csv 抽 N 顆變異，算它們的位置後比對既有 delta")
    args = ap.parse_args()

    out: Path = args.out_dir
    export_path = args.export_path or out / "evo2_table.js"
    seq = _load_sequence()
    n = len(seq)

    if args.status:
        status(out, n)
        return 0
    if args.export:
        export_table(out, n, export_path)
        return 0

    from utils.evo2_predictor import Evo2Predictor

    predictor = Evo2Predictor()

    scores, done = open_arrays(out, n)
    (out / "meta.json").write_text(json.dumps({
        "model": config.EVO2_MODEL,
        "window_size": predictor.window_size,
        "ttn_sequence_start": config.TTN_SEQUENCE_START,
        "ttn_sequence_end": config.TTN_SEQUENCE_END,
        "length": n,
        "bases": BASES,
        "index": "i = ttn_sequence_start - pos (minus strand; bases in transcript orientation)",
        "score": "mean log-likelihood of the window; delta = scores[i, alt] - scores[i, ref]",
        "pathogenic_threshold": config.PATHOGENIC_THRESHOLD,
        # 不同環境（如 gm3 的 FP8 vs gm4）算出的分數會差約 1e-4，記下來方便追查
        "environment": _env_info(),
    }, indent=2), encoding="utf-8")

    checkpoint = None
    if args.validate:
        import pandas as pd

        checkpoint = pd.read_csv(config.OUTPUT_DIR / "evo2_scoring_checkpoint.csv")
        checkpoint = checkpoint.sample(n=min(args.validate, len(checkpoint)), random_state=0)
        checkpoint["i"] = config.TTN_SEQUENCE_START - checkpoint["position"].astype(int)
        todo = sorted({int(i) for i in checkpoint["i"] if 0 <= i < n and not done[i]})
    else:
        todo = [int(i) for i in priority_order(n) if not done[i]]
    if args.limit:
        todo = todo[: args.limit]
    log = open(out / "progress.log", "a", encoding="utf-8")

    def say(msg: str) -> None:
        line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
        print(line, flush=True)
        log.write(line + "\n")
        log.flush()

    def report_validation() -> None:
        if checkpoint is None:
            return
        diffs, flips = [], 0
        for row in checkpoint.itertuples():
            i = int(row.i)
            delta = scores[i, BASES.index(str(row.alt).upper())] - scores[i, BASES.index(str(row.ref).upper())]
            diffs.append(abs(delta - float(row.delta_score)))
            flips += (delta < config.PATHOGENIC_THRESHOLD) != (float(row.delta_score) < config.PATHOGENIC_THRESHOLD)
        say(f"驗證 {len(diffs)} 顆：|Δ| 最大 {max(diffs):.2e}、中位數 {float(np.median(diffs)):.2e}，"
            f"致病性判定不同 {flips} 顆")

    if not todo:
        say("全部位置都已算完")
        report_validation()
        return 0

    stop = {"flag": False}

    def on_signal(signum, _frame):
        stop["flag"] = True
        say(f"收到訊號 {signum}，目前 batch 結束後存檔離開")

    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGINT, on_signal)

    import torch
    from evo2.scoring import logits_to_logprobs, prepare_batch

    say(f"載入 {config.EVO2_MODEL}；本次 {len(todo):,} 個位置，已完成 {int(done.sum()):,} / {n:,}")
    predictor._load_model()
    model, tokenizer = predictor.model.model, predictor.model.tokenizer
    half = predictor.window_size // 2

    t_start = time.time()
    finished = 0
    step = max(1, args.positions_per_batch)
    for b0 in range(0, len(todo), step):
        if stop["flag"] or (out / "STOP").exists():
            break
        group = todo[b0:b0 + step]
        # 視窗在序列兩端會被截短，長度不同就不能放同一個 batch
        batches, cur, cur_len = [], [], None
        for i in group:
            s, e = max(0, i - half), min(n, i + half)
            if cur and e - s != cur_len:
                batches.append(cur)
                cur = []
            cur.append(i)
            cur_len = e - s
        batches.append(cur)

        for positions in batches:
            seqs, slots = [], []
            for i in positions:
                s, e = max(0, i - half), min(n, i + half)
                ref = seq[s:e]
                k = i - s
                for b in BASES:
                    seqs.append(ref if b == ref[k] else ref[:k] + b + ref[k + 1:])
                    slots.append((i, BASES.index(b)))
            input_ids, lengths = prepare_batch(seqs, tokenizer, device="cuda:0")
            with torch.inference_mode():
                logits, _ = model(input_ids)
            lp = logits_to_logprobs(logits, input_ids).float().cpu().numpy()
            for row, (i, b) in enumerate(slots):
                scores[i, b] = float(np.mean(lp[row][: lengths[row]]))
            for i in positions:
                done[i] = 1
            finished += len(positions)

        if finished % 200 < step or b0 + step >= len(todo):
            scores.flush()
            done.flush()
            rate = (time.time() - t_start) / finished
            remaining = (n - int(done.sum())) * rate
            say(f"{int(done.sum()):,} / {n:,}（{rate:.2f} s/位置，全部剩約 {remaining / 86400:.1f} 天）")

    scores.flush()
    done.flush()
    say(f"本次結束：算了 {finished:,} 個位置，累計 {int(done.sum()):,} / {n:,}")
    report_validation()
    if checkpoint is None:
        export_table(out, n, export_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
