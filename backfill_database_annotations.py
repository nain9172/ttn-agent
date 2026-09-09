"""把 Curated Database Annotations 區塊回填進 outputs/ 既有的 HTML 報告。

報告本身已經有 ClinVar / Evo2 / 文獻 / 推理的內容，這裡只是多加一段
「這顆變異在 HGMD 與 PubTator 裡登記的疾病」，所以不需要重跑整個 pipeline
（不必載入 LLM、不必重算 Evo2），直接查表 + 注入 HTML 即可。

用法：
    python backfill_database_annotations.py                # 全部回填
    python backfill_database_annotations.py --dry-run      # 只看會改哪些、命中多少
    python backfill_database_annotations.py --limit 3      # 先試 3 份
    python backfill_database_annotations.py --outputs DIR  # 換目錄

重跑是安全的：已存在的區塊會被換成最新版本，不會重複堆疊。
"""

from __future__ import annotations

import argparse
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import config
from utils.html_report import HTMLReportGenerator
from utils.variant_db_lookup import VariantDatabaseLookup

logger = logging.getLogger("backfill_db")

RE_STYLE = re.compile(r"<style>.*?</style>", re.S)
RE_EXISTING = re.compile(r'\s*<section class="section dbann-section.*?</section>', re.S)
RE_VARIANT_ID = re.compile(r"^([\dXYM]+)-(\d+)-([ACGTN]+)-([ACGTN]+)$", re.I)
SUMMARY_MARKER = "<h2>Variant Summary</h2>"


def parse_variant_id(filename: str) -> Optional[str]:
    """報告檔名格式為 {chrom}-{pos}-{ref}-{alt}_{model}.html。"""
    stem = filename.split("_", 1)[0]
    return stem if RE_VARIANT_ID.match(stem) else None


def inject(text: str, section_html: str, css: str) -> str:
    """把區塊插在 Variant Summary 之後，並把 CSS 換成目前版本。"""
    text = RE_EXISTING.sub("", text, count=1)  # 重跑時先移除舊版本

    idx = text.find(SUMMARY_MARKER)
    if idx < 0:
        raise ValueError("找不到 Variant Summary 區塊，報告格式與預期不符")
    close = text.find("</section>", idx)
    if close < 0:
        raise ValueError("Variant Summary 區塊沒有結尾標籤")
    insert_at = close + len("</section>")

    text = RE_STYLE.sub(lambda _: f"<style>\n{css}\n    </style>", text, count=1)
    # 換 CSS 會改變字串長度，插入點要重算
    idx = text.find(SUMMARY_MARKER)
    close = text.find("</section>", idx)
    insert_at = close + len("</section>")

    return text[:insert_at] + "\n        " + section_html.strip() + text[insert_at:]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outputs", type=Path, default=config.OUTPUT_DIR, help="報告目錄")
    parser.add_argument("--limit", type=int, default=None, help="最多處理幾份（試跑用）")
    parser.add_argument("--dry-run", action="store_true", help="只統計，不寫檔")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    reports = sorted(args.outputs.glob("*.html"))
    if not reports:
        print(f"{args.outputs} 裡沒有 HTML 報告")
        return
    print(f"找到 {len(reports)} 份報告")

    lookup = VariantDatabaseLookup(
        hgmd_csv=config.HGMD_FULL_CSV_PATH,
        hgmd_kb=config.HGMD_KB_PATH,
        pubtator_csv=config.PUBTATOR_VARIANTS_CSV_PATH,
        cache_path=config.PUBTATOR_COORD_CACHE_PATH,
    )
    generator = HTMLReportGenerator()
    css = generator._get_css()

    stats: Counter = Counter()
    failures: List[Tuple[str, str]] = []
    disease_examples: List[str] = []

    todo = reports[: args.limit] if args.limit else reports
    for path in todo:
        variant_id = parse_variant_id(path.name)
        if not variant_id:
            stats["檔名無法解析變異"] += 1
            failures.append((path.name, "檔名不是 chrom-pos-ref-alt_model.html"))
            continue

        hits = lookup.lookup(variant_id=variant_id)
        if hits.hgmd_exact:
            stats["HGMD 精確命中"] += 1
        if hits.pubtator_exact:
            stats["PubTator 精確命中"] += 1
        if not hits.has_exact:
            stats["兩邊都沒精確命中"] += 1

        if len(disease_examples) < 5 and hits.has_exact:
            hgmd_d = ", ".join(d for d, _ in hits.hgmd_diseases()) or "—"
            pt_d = ", ".join(d for d, _ in hits.pubtator_diseases()) or "—"
            disease_examples.append(f"{variant_id}: HGMD [{hgmd_d}] · PubTator [{pt_d}]")

        text = path.read_text(encoding="utf-8", errors="ignore")
        section = generator._generate_database_section(hits)

        try:
            out_text = inject(text, section, css)
        except ValueError as e:
            stats["注入失敗"] += 1
            failures.append((path.name, str(e)))
            continue

        if not args.dry_run:
            tmp = path.with_suffix(".html.tmp")
            tmp.write_text(out_text, encoding="utf-8")
            tmp.replace(path)
        stats["已更新" if not args.dry_run else "可更新"] += 1

    print("\n=== 結果 ===")
    for key, n in stats.most_common():
        print(f"  {key:<20} {n:>5}")

    if disease_examples:
        print("\n命中範例：")
        for line in disease_examples:
            print(f"  - {line}")

    if failures:
        print(f"\n未處理 {len(failures)} 份：")
        for name, reason in failures[:10]:
            print(f"  - {name}: {reason}")


if __name__ == "__main__":
    main()
