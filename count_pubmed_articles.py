#!/usr/bin/env python3
"""以「文章」為單位統計 outputs 內 gemma-4-31B-it.html 報告。

每份報告對應一個輸入變異 (檔名，例如 2-178527031-T-C)。
每篇文章對應一個 <div class="article-title"> 區塊，內含：
  - 連結中的 PMID (pubmed.ncbi.nlm.nih.gov/<PMID>/)
  - source-badge：文章來源 (PubMed 或 LitVar2，皆來自 PubMed，合併計算)
  - mention-badge：Mention input variant: True / False

輸出：每篇 PubMed 文章 (PMID) 提及了幾個心臟/肌肉變異
      (即在多少份報告中被標記為 Mention input variant: True)。
"""

import csv
import glob
import os
import re

OUTPUTS_DIR = "/home/ryan910702/ttn-agent/outputs"
PATTERN = os.path.join(OUTPUTS_DIR, "*_gemma-4-31B-it.html")
CSV_OUT = "/home/ryan910702/ttn-agent/pubmed_article_stats.csv"

# 抓每篇文章標題行的 PMID、標題、來源與是否提及輸入變異
ARTICLE_RE = re.compile(
    r'<a href="https://pubmed\.ncbi\.nlm\.nih\.gov/(?P<pmid>\d+)/?"[^>]*>(?P<title>.*?)</a>'
    r'<span class="source-badge">(?P<source>[^<]*)</span>'
    r'<span class="mention-badge[^"]*">Mention input variant:\s*(?P<mention>True|False)</span>'
)


def variant_from_filename(path):
    name = os.path.basename(path)
    return name.replace("_gemma-4-31B-it.html", "")


def main():
    files = sorted(glob.glob(PATTERN))

    # pmid -> 資訊
    articles = {}
    # 至少被一篇文獻提及 (Mention input variant: True) 的變異
    variants_mentioned_by_lit = set()

    for path in files:
        variant = variant_from_filename(path)
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            html = f.read()

        for m in ARTICLE_RE.finditer(html):
            pmid = m.group("pmid")
            title = m.group("title").strip()
            mention_true = m.group("mention") == "True"

            info = articles.setdefault(
                pmid,
                {"title": title, "reports": set(), "mention_variants": set()},
            )
            if not info["title"]:
                info["title"] = title
            info["reports"].add(variant)
            if mention_true:
                info["mention_variants"].add(variant)
                variants_mentioned_by_lit.add(variant)

    # 依「提及的變異數量」由多到少排序
    rows = sorted(
        articles.items(),
        key=lambda kv: (len(kv[1]["mention_variants"]), len(kv[1]["reports"])),
        reverse=True,
    )

    with open(CSV_OUT, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "pmid",
            "title",
            "num_variants_mentioned",       # Mention input variant: True 的變異數
            "num_reports_appeared",          # 曾出現在幾份報告中
            "variants_mentioned",            # 被提及的變異清單
        ])
        for pmid, info in rows:
            writer.writerow([
                pmid,
                info["title"],
                len(info["mention_variants"]),
                len(info["reports"]),
                ";".join(sorted(info["mention_variants"])),
            ])

    total_articles = len(articles)
    mentioned_articles = sum(1 for _, i in rows if i["mention_variants"])
    print(f"報告數量 (gemma-4-31B-it.html)：{len(files)}")
    print(f"不重複 PubMed 文章數 (依 PMID)：{total_articles}")
    print(f"至少提及過 1 個心臟/肌肉變異的文章數：{mentioned_articles}")
    print(f"至少被文獻提及過 1 次的變異數：{len(variants_mentioned_by_lit)} / {len(files)}")
    print(f"\n每篇文章明細 (依提及變異數排序) 已輸出：{CSV_OUT}")
    print("\nTop 10 提及最多變異的文章：")
    for pmid, info in rows[:10]:
        print(f"  PMID {pmid}: {len(info['mention_variants'])} 個變異 - {info['title'][:70]}")


if __name__ == "__main__":
    main()
