#!/usr/bin/env python3
"""
批次跑 main.py，從 HGMD_TTN_processed.xlsx 讀取 input 欄位（chr-pos-ref-alt 格式），
為每個變異產生 HTML 報告，解析其中「Disease / tissue extracted from the article
for the target variant」區塊，把 8 個分類數字寫回 xlsx。

8 個欄位（依 mention input variant 區分 + 4 個 tissue 類別）：
  - mention_Cardiac
  - mention_Skeletal
  - mention_Both
  - mention_Not_Specified
  - no_mention_Cardiac
  - no_mention_Skeletal
  - no_mention_Both
  - no_mention_Not_Specified

特性：
  - 每跑完一個變異就立刻 save xlsx，避免 SSH 中斷後遺失進度
  - 已經有結果（8 個欄位都填好）的變異會跳過
  - 若 outputs/ 目錄下已存在對應 gemma-4 HTML 報告，直接解析、不重跑 main.py
  - Ctrl+C 會等本筆完成並儲存再退出
"""

from __future__ import annotations

import argparse
import re
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
from bs4 import BeautifulSoup

PROJECT_ROOT = Path(__file__).parent.resolve()
XLSX_PATH = PROJECT_ROOT / "HGMD_TTN_processed.xlsx"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
MAIN_PY = PROJECT_ROOT / "main.py"

MODEL_TAG = "gemma-4-31B-it"

INPUT_COL = "input"

RESULT_COLS = [
    "mention_Cardiac",
    "mention_Skeletal",
    "mention_Both",
    "mention_Not_Specified",
    "no_mention_Cardiac",
    "no_mention_Skeletal",
    "no_mention_Both",
    "no_mention_Not_Specified",
]

TISSUE_LABEL_TO_KEY = {
    "Cardiac": "Cardiac",
    "Skeletal": "Skeletal",
    "Both": "Both",
    "Not Specified": "Not_Specified",
}

VARIANT_RE = re.compile(r"^\d+\-\d+\-[ACGT]+\-[ACGT]+$", re.IGNORECASE)

_stop_requested = False


def _signal_handler(signum, frame):  # noqa: ANN001
    global _stop_requested
    _stop_requested = True
    print("\n[!] 收到中斷信號，跑完目前這筆並儲存後就結束 (再按一次 Ctrl+C 強制離開)…",
          flush=True)
    signal.signal(signal.SIGINT, signal.default_int_handler)


def find_existing_report(variant_id: str) -> Optional[Path]:
    """找已存在的 gemma-4 報告。"""
    candidate = OUTPUT_DIR / f"{variant_id}_{MODEL_TAG}.html"
    if candidate.exists():
        return candidate
    return None


def parse_report(html_path: Path) -> Optional[dict]:
    """解析 HTML 報告，取出 8 個 Disease/tissue 計數。

    回傳格式為 {RESULT_COLS[i]: int}，若整份報告沒有 Tissue Involvement Statistics
    區塊則回傳 None（代表這份報告不完整，要求重跑）。
    """
    try:
        html = html_path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        print(f"  [warn] 讀檔失敗 {html_path.name}: {exc}", flush=True)
        return None

    if "Tissue Involvement Statistics" not in html:
        return None

    soup = BeautifulSoup(html, "html.parser")

    counts = {col: 0 for col in RESULT_COLS}

    headers = soup.find_all("h4")
    for h4 in headers:
        title = h4.get_text(" ", strip=True)
        if title.startswith("Articles that directly mention the input variant"):
            prefix = "mention_"
        elif title.startswith("Articles that do NOT directly mention the input variant"):
            prefix = "no_mention_"
        else:
            continue

        container = h4.parent
        if container is None:
            continue

        no_articles_p = container.find(
            "p", string=lambda s: s and "No articles in this category" in s
        )
        if no_articles_p is not None:
            for label, key in TISSUE_LABEL_TO_KEY.items():
                counts[f"{prefix}{key}"] = 0
            continue

        for card in container.select(".info-card"):
            label_el = card.select_one(".label")
            value_el = card.select_one(".value")
            if not label_el or not value_el:
                continue
            label = label_el.get_text(strip=True)
            if label not in TISSUE_LABEL_TO_KEY:
                continue
            m = re.match(r"\s*(\d+)", value_el.get_text(strip=True))
            if not m:
                continue
            counts[f"{prefix}{TISSUE_LABEL_TO_KEY[label]}"] = int(m.group(1))

    return counts


def run_main_py(variant_id: str, timeout: Optional[int] = None) -> bool:
    """執行 python main.py <variant>。回傳是否成功（return code 0）。"""
    cmd = [sys.executable, str(MAIN_PY), variant_id]
    print(f"  [run] {' '.join(cmd)}", flush=True)
    try:
        result = subprocess.run(
            cmd,
            cwd=str(PROJECT_ROOT),
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        print(f"  [error] {variant_id} 執行超時 ({timeout}s)", flush=True)
        return False
    except KeyboardInterrupt:
        raise
    except Exception as exc:  # noqa: BLE001
        print(f"  [error] {variant_id} 執行失敗: {exc}", flush=True)
        return False
    return result.returncode == 0


def ensure_columns(df: pd.DataFrame) -> pd.DataFrame:
    """確保 8 個結果欄位存在，並使用 nullable Int64 型別。"""
    for col in RESULT_COLS:
        if col not in df.columns:
            df[col] = pd.NA
        df[col] = df[col].astype("Int64")
    return df


def save_xlsx(df: pd.DataFrame, path: Path) -> None:
    """以暫存檔 + rename 方式安全儲存，避免寫到一半中斷造成損壞。"""
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_excel(tmp, index=False)
    tmp.replace(path)


def row_is_done(row: pd.Series) -> bool:
    return all(pd.notna(row.get(c)) for c in RESULT_COLS)


def write_counts(df: pd.DataFrame, idx: int, counts: dict) -> None:
    for col in RESULT_COLS:
        df.at[idx, col] = int(counts.get(col, 0))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xlsx",
        type=Path,
        default=XLSX_PATH,
        help=f"輸入/輸出的 xlsx 路徑 (預設: {XLSX_PATH})",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="從第幾列開始 (0-based, 不含表頭)",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="跑到第幾列為止 (exclusive)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="最多跑幾筆 (在 skip 後計算)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="每筆 main.py 執行的逾時秒數 (預設不限)",
    )
    parser.add_argument(
        "--parse-only",
        action="store_true",
        help="只解析已存在的 HTML 報告，不跑 main.py（用來補回先前產生過的結果）",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="忽略既有結果，全部重跑",
    )
    parser.add_argument(
        "--variant",
        type=str,
        default=None,
        help="只跑單一變異 (例如 2-178527031-T-C)，會更新對應的列",
    )
    args = parser.parse_args()

    xlsx_path: Path = args.xlsx.resolve()
    if not xlsx_path.exists():
        print(f"[error] 找不到 xlsx: {xlsx_path}", file=sys.stderr)
        return 1
    if not MAIN_PY.exists():
        print(f"[error] 找不到 main.py: {MAIN_PY}", file=sys.stderr)
        return 1

    print(f"[info] 讀取 {xlsx_path}", flush=True)
    df = pd.read_excel(xlsx_path)
    if INPUT_COL not in df.columns:
        print(f"[error] xlsx 缺少 '{INPUT_COL}' 欄位", file=sys.stderr)
        return 1

    df = ensure_columns(df)
    save_xlsx(df, xlsx_path)

    signal.signal(signal.SIGINT, _signal_handler)

    total_rows = len(df)
    start = max(0, args.start)
    end = total_rows if args.end is None else min(args.end, total_rows)

    indices = list(range(start, end))
    if args.variant is not None:
        target = args.variant.strip()
        match_idx = df.index[df[INPUT_COL].astype(str) == target].tolist()
        if not match_idx:
            print(f"[error] 在 xlsx 中找不到變異 {target}", file=sys.stderr)
            return 1
        indices = match_idx

    processed = 0
    skipped = 0
    parsed_existing = 0
    ran_main = 0
    failed = 0

    for idx in indices:
        if _stop_requested:
            break

        if args.limit is not None and processed >= args.limit:
            break

        row = df.loc[idx]
        variant_id = str(row[INPUT_COL]).strip()
        if not variant_id or variant_id.lower() == "nan":
            continue
        if not VARIANT_RE.match(variant_id):
            print(f"[skip] row {idx}: '{variant_id}' 不符合 chr-pos-ref-alt 格式",
                  flush=True)
            continue

        if not args.force and row_is_done(row):
            skipped += 1
            continue

        print(f"\n[{idx + 1}/{total_rows}] variant = {variant_id}", flush=True)
        t0 = time.time()

        report_path = find_existing_report(variant_id)
        counts: Optional[dict] = None

        if report_path is not None and not args.force:
            print(f"  [reuse] 找到既有報告 {report_path.name}", flush=True)
            counts = parse_report(report_path)
            if counts is not None:
                parsed_existing += 1
            else:
                print("  [warn] 既有報告沒有 Tissue Involvement Statistics 區塊，將重跑",
                      flush=True)

        if counts is None:
            if args.parse_only:
                print("  [skip] --parse-only 模式且沒有可用既有報告", flush=True)
                continue

            ok = run_main_py(variant_id, timeout=args.timeout)
            if not ok:
                failed += 1
                print(f"  [error] main.py 失敗，跳過 {variant_id}", flush=True)
                continue

            report_path = find_existing_report(variant_id)
            if report_path is None:
                failed += 1
                print(f"  [error] 找不到產生的報告 {variant_id}_{MODEL_TAG}.html",
                      flush=True)
                continue

            counts = parse_report(report_path)
            if counts is None:
                failed += 1
                print("  [error] 解析新產生的報告失敗（沒有 Tissue Involvement Statistics）",
                      flush=True)
                continue

            ran_main += 1

        write_counts(df, idx, counts)
        save_xlsx(df, xlsx_path)
        processed += 1

        dt = time.time() - t0
        summary = ", ".join(f"{c.replace('mention_', 'M:').replace('no_M:', 'N:')}={counts[c]}"
                            for c in RESULT_COLS)
        print(f"  [ok] {variant_id} 完成 ({dt:.1f}s) -> {summary}", flush=True)

    print("\n========== 完成 ==========", flush=True)
    print(f"已處理: {processed} 筆", flush=True)
    print(f"  - 解析既有報告: {parsed_existing}", flush=True)
    print(f"  - 新跑 main.py: {ran_main}", flush=True)
    print(f"略過 (已有結果): {skipped}", flush=True)
    print(f"失敗: {failed}", flush=True)
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
