#!/usr/bin/env python3
"""實驗二／三的第二階段：用「一種設定」跑完整個評估集的臨床抽取（需要 GPU）。

輸入是 build_llm_corpus.py 凍結好的語料，所以每個設定看到的文獻完全相同，差異
只來自這裡的旗標。一次執行只載入一個模型、只用一種設定，結果逐顆變異 append 到
JSONL，中斷後可續跑。

實驗二（換模型）：
  python research/experiments/run_extraction_config.py --label gemma4-31b-nvfp4 \
      --model nvidia/Gemma-4-31B-IT-NVFP4
  python research/experiments/run_extraction_config.py --label medgemma-27b \
      --model google/medgemma-27b-it

實驗三（關掉防幻覺機制，一次關一個）：
  python research/experiments/run_extraction_config.py --label all-guards          # 全開（基準）
  python research/experiments/run_extraction_config.py --label no-alias-gate  --no-alias-gate
  python research/experiments/run_extraction_config.py --label no-guided-json --no-guided-json
  python research/experiments/run_extraction_config.py --label no-evidence    --no-verify-evidence
  python research/experiments/run_extraction_config.py --label n1             --self-consistency 1
  python research/experiments/run_extraction_config.py --label none           \
      --no-alias-gate --no-guided-json --no-verify-evidence --self-consistency 1

每個機制的開關都是改 utils.local_clinical_extractor 的模組層常數（那裡是 from
config import 進去的名字，改 config 物件沒有用），而且一定要在建立 extractor
之前改：self-consistency 的 n 與 guided JSON 在 __init__ 就寫進 SamplingParams。

不論旗標怎麼設，每篇文獻都會另外算兩個「事後探針」，機制關掉時才能量出它本來
會擋掉什麼：
  alias_present            alias gate 的判定（gate 關掉時也算）
  evidence_grounded_probe  evidence 與原文的 5-gram 重疊是否過關（同上）
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from copy import deepcopy
from pathlib import Path
from typing import Dict, List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

CORPUS_DIR = PROJECT_ROOT / "research" / "experiments" / "corpus"
RUNS_DIR = PROJECT_ROOT / "research" / "experiments" / "round1_random60" / "runs"
DEFAULT_VARIANTS = PROJECT_ROOT / "research" / "experiments" / "round1_random60" / "llm_eval_set.csv"


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", required=True, help="這個設定的名字，也是輸出檔名")
    ap.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    ap.add_argument("--corpus-dir", type=Path, default=CORPUS_DIR)
    ap.add_argument("--out-dir", type=Path, default=RUNS_DIR)
    ap.add_argument("--limit", type=int, default=None, help="只跑前 N 顆")
    ap.add_argument("--force", action="store_true", help="重跑，覆蓋既有輸出")
    ap.add_argument("--pmid-column", default=None,
                    help=("評估集裡列出「這顆變異要用哪些文獻」的欄位（分號分隔的 PMID，"
                          "例如 find_hub_articles.py 產生的 hub_pmids）。給了就只餵這些文獻，"
                          "語料裡其他文獻一律不送進 LLM"))

    ap.add_argument("--backend", default=None, help="預設取 config.LOCAL_LLM_BACKEND")
    ap.add_argument("--model", default=None, help="預設取 config.LOCAL_LLM_MODEL")
    ap.add_argument("--max-model-len", type=int, default=None)
    ap.add_argument("--max-context-length", type=int, default=None)
    ap.add_argument("--gpu-memory-utilization", type=float, default=None,
                    help=("覆寫 config.LLM_GPU_MEMORY_UTILIZATION。預設 0.85 是「總記憶體的 85%%」，"
                          "機器上還有別的東西佔著時會直接 OOM；所有設定要用同一個值才公平"))

    # 防幻覺機制：預設沿用 config.py 的值，旗標只用來關掉／改數字
    ap.add_argument("--self-consistency", type=int, default=None,
                    help="LLM_SELF_CONSISTENCY_N（1 = 關掉 self-consistency）")
    ap.add_argument("--temperature", type=float, default=None)
    ap.add_argument("--no-guided-json", dest="guided_json", action="store_false", default=None)
    ap.add_argument("--guided-json", dest="guided_json", action="store_true", default=None)
    ap.add_argument("--no-verify-evidence", dest="verify_evidence", action="store_false", default=None)
    ap.add_argument("--verify-evidence", dest="verify_evidence", action="store_true", default=None)
    ap.add_argument("--no-alias-gate", dest="alias_gate", action="store_false", default=None)
    ap.add_argument("--alias-gate", dest="alias_gate", action="store_true", default=None)
    ap.add_argument("--evidence-threshold", type=float, default=None,
                    help="LLM_EVIDENCE_NGRAM_OVERLAP_THRESHOLD")
    return ap


def apply_knobs(args) -> Dict[str, object]:
    """把旗標寫進 utils.local_clinical_extractor 的模組常數，回傳實際生效的設定。"""
    import utils.local_clinical_extractor as lce

    if args.self_consistency is not None:
        lce.LLM_SELF_CONSISTENCY_N = max(1, args.self_consistency)
    if args.temperature is not None:
        lce.LLM_SAMPLING_TEMPERATURE = args.temperature
    if args.guided_json is not None:
        lce.LLM_USE_GUIDED_JSON = args.guided_json
    if args.verify_evidence is not None:
        lce.LLM_VERIFY_EVIDENCE = args.verify_evidence
    if args.alias_gate is not None:
        lce.LLM_REQUIRE_DETERMINISTIC_ALIAS = args.alias_gate
    if args.evidence_threshold is not None:
        lce.LLM_EVIDENCE_NGRAM_OVERLAP_THRESHOLD = args.evidence_threshold

    return {
        "self_consistency_n": lce.LLM_SELF_CONSISTENCY_N,
        "temperature": lce.LLM_SAMPLING_TEMPERATURE,
        "guided_json": lce.LLM_USE_GUIDED_JSON,
        "verify_evidence": lce.LLM_VERIFY_EVIDENCE,
        "alias_gate": lce.LLM_REQUIRE_DETERMINISTIC_ALIAS,
        "evidence_threshold": lce.LLM_EVIDENCE_NGRAM_OVERLAP_THRESHOLD,
    }


class JsonFailureCounter:
    """包住 _robust_json_parse，數 schema / JSON 破損次數（guided JSON 的效果）。"""

    def __init__(self, extractor):
        self.extractor = extractor
        self.failures = 0
        self.parsed = 0
        self._orig = extractor._robust_json_parse

        def wrapped(text: str):
            self.parsed += 1
            result = self._orig(text)
            if result.get("reasoning") == "JSON Parsing Failed":
                self.failures += 1
            return result

        extractor._robust_json_parse = wrapped  # type: ignore[method-assign]

    def take(self) -> Dict[str, int]:
        out = {"json_parse_failures": self.failures, "samples_parsed": self.parsed}
        self.failures = self.parsed = 0
        return out


def raw_text_of(article: Dict) -> str:
    """與 batch_extract 相同的取文順序。"""
    return (article.get("text_for_llm") or article.get("full_text")
            or article.get("abstract") or "")


def article_record(extractor, article: Dict, aliases: List[str]) -> Dict:
    """一篇文獻的結果 + 兩個事後探針。"""
    info = article.get("clinical_info") or {}
    raw_text = raw_text_of(article)
    gate_text = article.get("gate_text") or raw_text
    source = str(info.get("extraction_source", ""))

    try:
        alias_present, matched = extractor._alias_present_in_text(aliases, gate_text, article)
    except Exception:  # noqa: BLE001 — 探針失敗不該讓整批停下來
        alias_present, matched = None, None

    evidence = str(info.get("evidence_sentence", "") or "")
    grounded_probe: Optional[bool] = None
    overlap_probe: Optional[float] = None
    if evidence and evidence.strip().lower() not in ("n/a", "not specified", ""):
        try:
            grounded_probe, score = extractor._verify_evidence_grounding(evidence, raw_text, article)
            overlap_probe = round(float(score), 3)
        except Exception:  # noqa: BLE001
            pass

    # 這篇有沒有真的送進 LLM（gate 擋掉或文本太短就沒有）
    prompted = source != "deterministic_no_alias" and info.get("reasoning") != "Text too short"

    return {
        "pmid": article.get("pmid"),
        "disease": info.get("disease"),
        "tissue_affected": info.get("tissue_affected"),
        "patient_count": info.get("patient_count"),
        "evidence_sentence": evidence[:400],
        "extraction_source": source,
        "evidence_grounded": info.get("evidence_grounded"),
        "evidence_overlap_score": info.get("evidence_overlap_score"),
        "alias_present": alias_present,
        "alias_matched": matched,
        "evidence_grounded_probe": grounded_probe,
        "evidence_overlap_probe": overlap_probe,
        "prompted": prompted,
        "has_full_text": bool(article.get("has_full_text")),
        "text_chars": len(raw_text),
    }


def main() -> int:
    args = build_parser().parse_args()

    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
                        stream=sys.stdout)

    import pandas as pd
    from config import (LOCAL_LLM_BACKEND, LOCAL_LLM_MAX_CONTEXT_LENGTH,
                        LOCAL_LLM_MAX_MODEL_LEN, LOCAL_LLM_MODEL,
                        LOCAL_LLM_TENSOR_PARALLEL)

    if not args.variants.exists():
        print(f"[error] 找不到評估集 {args.variants}", file=sys.stderr)
        return 1
    variants_df = pd.read_csv(args.variants)
    variants = [str(v).strip() for v in variants_df["input"].dropna()]
    if args.limit:
        variants = variants[: args.limit]

    # 每顆變異的文獻白名單：只餵指定的那幾篇（hub 文獻實驗用）
    pmid_allow: Dict[str, set] = {}
    if args.pmid_column:
        if args.pmid_column not in variants_df.columns:
            print(f"[error] 評估集沒有 {args.pmid_column} 欄", file=sys.stderr)
            return 1
        for _, row in variants_df.iterrows():
            raw = str(row.get(args.pmid_column) or "").strip()
            if raw and raw.lower() != "nan":
                pmid_allow[str(row["input"]).strip()] = {
                    p.strip() for p in raw.replace(",", ";").split(";") if p.strip()
                }
        print(f"[info] 只餵 {args.pmid_column} 指定的文獻："
              f"{sum(len(v) for v in pmid_allow.values())} 個 (變異, 文獻) 配對", flush=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / f"{args.label}.jsonl"
    if out_path.exists() and args.force:
        out_path.unlink()
    done = set()
    if out_path.exists():
        with open(out_path, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    done.add(json.loads(line)["variant"])
        print(f"[info] 續跑：{len(done)} 顆已完成", flush=True)

    todo = [v for v in variants if v not in done]
    missing = [v for v in todo if not (args.corpus_dir / f"{v}.json").exists()]
    if missing:
        print(f"[warn] {len(missing)} 顆沒有語料檔，會被跳過；先跑 build_llm_corpus.py"
              f"（例：{missing[:3]}）", flush=True)
        todo = [v for v in todo if v not in set(missing)]
    if not todo:
        print("[info] 沒有待處理的變異")
        return 0

    knobs = apply_knobs(args)          # 一定要在建立 extractor 之前
    if args.gpu_memory_utilization is not None:
        # _init_vllm 是在函式裡 from config import，所以改模組屬性就會生效
        import config as _config
        _config.LLM_GPU_MEMORY_UTILIZATION = args.gpu_memory_utilization
        print(f"[info] gpu_memory_utilization={args.gpu_memory_utilization}")
    backend = args.backend or LOCAL_LLM_BACKEND
    model = args.model or LOCAL_LLM_MODEL
    print(f"[info] label={args.label} model={model} backend={backend}")
    print(f"[info] knobs={knobs}")
    print(f"[info] 待處理 {len(todo)} 顆 -> {out_path}", flush=True)

    from utils.local_clinical_extractor import LocalClinicalExtractor

    extractor = LocalClinicalExtractor(
        backend=backend,
        model_name=model,
        tensor_parallel_size=LOCAL_LLM_TENSOR_PARALLEL,
        max_model_len=args.max_model_len or LOCAL_LLM_MAX_MODEL_LEN,
        max_context_length=args.max_context_length or LOCAL_LLM_MAX_CONTEXT_LENGTH,
    )
    counter = JsonFailureCounter(extractor)

    sys.path.insert(0, str(PROJECT_ROOT / "research" / "experiments"))
    from llm_eval_common import tissue_counts

    for i, vid in enumerate(todo, 1):
        payload = json.loads((args.corpus_dir / f"{vid}.json").read_text(encoding="utf-8"))
        articles = deepcopy(payload["articles"])          # batch_extract 會就地改寫
        if vid in pmid_allow:
            keep = pmid_allow[vid]
            before = len(articles)
            articles = [a for a in articles if str(a.get("pmid")) in keep]
            if not articles:
                print(f"  [skip] {vid}: 語料裡找不到指定的文獻 {sorted(keep)}", flush=True)
                continue
            print(f"  [filter] {before} 篇 -> {len(articles)} 篇（只留 {sorted(keep)}）", flush=True)
        variant_info = payload["variant_info"]
        clinvar_info = payload.get("clinvar_info")
        aliases = payload.get("aliases") or []

        print(f"[{i}/{len(todo)}] {vid}（{len(articles)} 篇）", flush=True)
        t0 = time.time()
        try:
            enriched = extractor.batch_extract(articles, variant_info, clinvar_info)
        except Exception as exc:  # noqa: BLE001
            print(f"  [error] {vid} 抽取失敗: {exc}", flush=True)
            continue
        elapsed = time.time() - t0

        record = {
            "variant": vid,
            "config": args.label,
            "model": model,
            "knobs": knobs,
            "seconds": round(elapsed, 2),
            "counts": tissue_counts(enriched),
            "articles": [article_record(extractor, a, aliases) for a in enriched],
            **counter.take(),
        }
        with open(out_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(record, ensure_ascii=False, default=str) + "\n")
        print(f"  -> {elapsed:.1f}s，mention={[record['counts'][k] for k in sorted(record['counts']) if k.startswith('mention')]}",
              flush=True)

    print(f"\n完成，輸出：{out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
