#!/usr/bin/env python3
"""實驗二／三的第一階段：把評估集每顆變異的文獻抓下來存成 JSON（不用 GPU）。

為什麼要單獨一個階段：實驗二要比兩三個模型，實驗三要比五六種機制組合，如果每
個設定都從 main.py 從頭跑，同一批文獻會被重複下載與 Docling 轉檔好幾次，而且
兩個設定看到的輸入還不保證一樣（PubMed / LitVar 回傳會變）。先把語料凍結成檔
案，之後每個設定都吃同一份輸入，差異才只來自模型或機制本身。

抓取流程與 main.run_pipeline 的 Step 2 / 4 / 4a 相同：
  ClinVar（取 PMID、rsID、HGVS 別名）→ EnhancedPubMedSearcher（全文 + Docling）
  → LitVarSearcher（rsID 找更多文獻）→ 依 PMID 去重

用法（先 conda activate gm4）：
  python research/experiments/build_llm_corpus.py                      # 用 llm_eval_set.csv
  python research/experiments/build_llm_corpus.py --limit 5 --no-full-text   # 快速試跑
  python research/experiments/build_llm_corpus.py --variants my.csv --force

輸出：research/experiments/corpus/<variant>.json，一顆變異一檔，已存在就跳過（可續跑）。

注意：Docling 轉 PDF 很吃 host 記憶體，不要跟 vLLM 批次（run_batch_inprocess.sh）
同時跑，128 GB unified memory 會被兩邊搶。
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

from config import (
    ENABLE_DOCLING_PDF,
    ENABLE_LITVAR_SEARCH,
    LITVAR_MAX_RESULTS,
    MAX_TEXT_LENGTH,
)
from utils.clinvar_parser import ClinVarParser
from utils.enhanced_pubmed_search import EnhancedPubMedSearcher
from utils.litvar_search import LitVarSearcher
from utils.variant_parser import parse_variant
from utils.variant_utils import get_variant_aliases

CORPUS_DIR = PROJECT_ROOT / "research" / "experiments" / "corpus"
DEFAULT_VARIANTS = PROJECT_ROOT / "research" / "experiments" / "llm_eval_set.csv"

logger = logging.getLogger("build_llm_corpus")

# 只留抽取階段真的會用到的欄位，JSON 才不會爆大
ARTICLE_KEYS = [
    "pmid", "title", "abstract", "journal", "year", "authors", "pubmed_link",
    "full_text", "text_for_llm", "gate_text", "source", "has_full_text",
    "supplementary_files", "docling_markdown_path",
]


def fetch_one(variant_id: str, use_full_text: bool, use_docling: bool) -> Dict:
    """回傳 {variant_info, clinvar_info, articles}；與 run_pipeline 的抓法一致。"""
    variant_info = parse_variant(variant_id)

    clinvar_info = ClinVarParser().parse_variant(variant_info) or None

    searcher = EnhancedPubMedSearcher(
        try_full_text=use_full_text,
        max_text_length=MAX_TEXT_LENGTH,
        use_docling=use_docling and use_full_text,
    )
    vid = f"{variant_info['chrom']}-{variant_info['pos']}-{variant_info['ref']}-{variant_info['alt']}"
    aliases = get_variant_aliases(vid, clinvar_info)

    articles: List[Dict] = []
    pmid_list = (clinvar_info or {}).get("pmid_list") or []
    if pmid_list:
        articles = searcher.search(variant_info, pmid_list=pmid_list, variant_aliases=aliases)

    if ENABLE_LITVAR_SEARCH:
        litvar = LitVarSearcher(use_docling=use_docling and use_full_text).search_multiple_formats(
            variant_info=variant_info,
            clinvar_info=clinvar_info,
            variant_aliases=aliases,
            max_results=LITVAR_MAX_RESULTS,
        )
        seen = {a.get("pmid") for a in articles}
        articles.extend(a for a in litvar if a.get("pmid") not in seen)

    slim = [{k: a.get(k) for k in ARTICLE_KEYS if k in a} for a in articles]
    return {
        "variant": vid,
        "variant_info": variant_info,
        "clinvar_info": clinvar_info,
        "aliases": aliases,
        "articles": slim,
        "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "full_text": use_full_text,
        "docling": use_docling and use_full_text,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS,
                    help="評估集 CSV，需有 input 欄（預設 research/experiments/llm_eval_set.csv）")
    ap.add_argument("--out-dir", type=Path, default=CORPUS_DIR)
    ap.add_argument("--limit", type=int, default=None, help="只抓前 N 顆（試跑用）")
    ap.add_argument("--force", action="store_true", help="已存在的也重抓")
    ap.add_argument("--no-full-text", action="store_true",
                    help="只抓摘要，不下載 PDF、不跑 Docling（快速煙霧測試用）")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
                        stream=sys.stdout)

    if not args.variants.exists():
        print(f"[error] 找不到 {args.variants}；先在 notebook 裡用 "
              f"llm_eval_common.build_eval_set() 產生", file=sys.stderr)
        return 1

    df = pd.read_csv(args.variants)
    if "input" not in df.columns:
        print(f"[error] {args.variants} 沒有 input 欄", file=sys.stderr)
        return 1
    variants = [str(v).strip() for v in df["input"].dropna()]
    if args.limit:
        variants = variants[: args.limit]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    use_full_text = not args.no_full_text

    done = skipped = failed = 0
    for i, vid in enumerate(variants, 1):
        out_path = args.out_dir / f"{vid}.json"
        if out_path.exists() and not args.force:
            skipped += 1
            continue
        print(f"[{i}/{len(variants)}] {vid}", flush=True)
        t0 = time.time()
        try:
            payload = fetch_one(vid, use_full_text, ENABLE_DOCLING_PDF)
        except Exception as exc:  # noqa: BLE001 — 單顆失敗不該停掉整批
            failed += 1
            logger.error(f"{vid} 抓取失敗: {exc}", exc_info=True)
            continue
        out_path.write_text(json.dumps(payload, ensure_ascii=False, default=str),
                            encoding="utf-8")
        done += 1
        print(f"  -> {len(payload['articles'])} 篇，{time.time() - t0:.1f}s", flush=True)

    print(f"\n完成 {done} 顆；略過（已存在）{skipped} 顆；失敗 {failed} 顆")
    print(f"語料目錄：{args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
