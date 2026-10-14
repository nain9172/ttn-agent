#!/usr/bin/env python3
"""檢查 PubTator 上是否有新的 TTN 文獻，決定要不要重跑 pubtator_harvest。

比對 PubTator3 目前 @GENE_TTN 的命中篇數與上次 harvest 時記下的篇數
（outputs/pubtator_cache/harvest_state.json）。還沒有狀態檔時，改讀
outputs/pubtator_harvest.log 裡最後一次的「@GENE_TTN 命中 N 篇」當基準。

  python scripts/data_prep/pubtator_check_new.py            # exit 0 = 有新文獻（或從未 harvest），10 = 已是最新
  python scripts/data_prep/pubtator_check_new.py --record   # harvest 完成後記下目前篇數
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STATE = ROOT / "outputs" / "pubtator_cache" / "harvest_state.json"
HARVEST_LOG = ROOT / "outputs" / "pubtator_harvest.log"
SEARCH = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text=@GENE_TTN&page=1"
UP_TO_DATE = 10


def current_count() -> int:
    with urllib.request.urlopen(SEARCH, timeout=60) as resp:
        return int(json.load(resp)["count"])


def last_count():
    if STATE.exists():
        state = json.loads(STATE.read_text())
        return int(state["gene_ttn_count"]), state.get("date", "?")
    if HARVEST_LOG.exists():
        hits = re.findall(r"^(\S+) .*@GENE_TTN 命中 (\d+) 篇", HARVEST_LOG.read_text(errors="replace"), re.M)
        if hits:
            return int(hits[-1][1]), hits[-1][0]
    return None, None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--record", action="store_true", help="記下目前的篇數（harvest 完成後呼叫）")
    args = ap.parse_args()

    now = current_count()
    if args.record:
        STATE.parent.mkdir(parents=True, exist_ok=True)
        STATE.write_text(json.dumps({"gene_ttn_count": now, "date": time.strftime("%Y-%m-%d")}, indent=2))
        print(f"[pubtator] 已記錄：@GENE_TTN {now:,} 篇")
        return 0

    before, when = last_count()
    if before is None:
        print(f"[pubtator] 尚未 harvest 過；PubTator 目前 @GENE_TTN {now:,} 篇")
        return 0
    if now > before:
        print(f"[pubtator] 有新文獻：{when} 為 {before:,} 篇，現在 {now:,} 篇（+{now - before:,}）")
        return 0
    print(f"[pubtator] 已是最新：{when} 為 {before:,} 篇，現在 {now:,} 篇")
    return UP_TO_DATE


if __name__ == "__main__":
    sys.exit(main())
