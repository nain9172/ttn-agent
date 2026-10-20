#!/usr/bin/env python3
"""遺傳模式推論（utils/inheritance_inference.py）對 HGMD 的 leave-one-out 驗證。

標籤：HGMD comments 寫 homozygous / compound het → AR；否則依登記的疾病（OMIM 遺傳模式）。
HGMD 的 expected_inheritance 欄全部是 AD（基因層級預設），不拿來當答案。

每筆紀錄只給推論器「位置、變異類型」：不給它自己的疾病與紀錄（避免把答案讀回來），
類比也排除同一位置。--llm 時另外用 Gemma（vLLM）跑同一批，比較規則式與 LLM。

    python research/validate_inheritance.py                 # 只跑規則式
    python research/validate_inheritance.py --llm           # 規則式 + LLM
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from utils.inheritance_inference import AD, AR, BOTH, SHORT, UNCERTAIN, InheritanceInference  # noqa: E402


def summarize(name: str, rows, key: str) -> dict:
    n = len(rows)
    strict = sum(r[key] == r["label"] for r in rows)
    # BOTH 標籤時預測 AD 或 AR 都算對一半的情況另外列；這裡 lenient = 標籤是 BOTH 或預測是 BOTH 時算對
    lenient = sum(r[key] == r["label"] or BOTH in (r[key], r["label"]) for r in rows)
    called = [r for r in rows if r[key] != UNCERTAIN]
    per = {}
    for lab in (AD, AR, BOTH):
        sub = [r for r in rows if r["label"] == lab]
        if sub:
            per[SHORT[lab]] = {"n": len(sub), "recall": round(sum(r[key] == lab for r in sub) / len(sub), 3),
                               "pred": dict(Counter(SHORT[r[key]] for r in sub))}
    recalls = [v["recall"] for k, v in per.items() if k in ("AD", "AR")]
    out = {
        "method": name, "n": n,
        "accuracy": round(strict / n, 3),
        "accuracy_allelic_lenient": round(lenient / n, 3),
        "balanced_accuracy_AD_AR": round(sum(recalls) / len(recalls), 3) if recalls else None,
        "coverage": round(len(called) / n, 3),
        "accuracy_when_called": round(sum(r[key] == r["label"] for r in called) / len(called), 3) if called else None,
        "per_label": per,
    }
    print(json.dumps(out, indent=2, ensure_ascii=False))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--llm", action="store_true", help="also run the LLM (needs ~104 GiB free for vLLM)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", type=Path, default=ROOT / "research" / "experiments" / "inheritance_validation.json")
    args = ap.parse_args()

    import main as pipeline

    mi = pipeline.build_mechanistic_engine()
    inh = InheritanceInference(mi, pipeline.INFERENCE_HGMD_CSV_PATH)
    records = inh.labelled_hgmd()
    # 同一顆變異（同位置同 HGVS）在 HGMD 有多筆時只算一次，避免重複計分
    seen, uniq = set(), []
    for r in records:
        k = (r["pos"], r["dna"])
        if k not in seen:
            seen.add(k)
            uniq.append(r)
    if args.limit:
        uniq = uniq[: args.limit]
    print(f"[inheritance] {len(records)} labelled HGMD records → {len(uniq)} unique variants; "
          f"labels {dict(Counter(SHORT[r['label']] for r in uniq))}", flush=True)

    rows, evs = [], []
    for r in uniq:
        feats = {"pos": r["pos"], "exon": r["exon"], "band": r["band"], "vgroup": r["vgroup"], "vclass": r["vclass"]}
        ev = inh.collect(feats, exclude_acc=[r["acc"]])
        rule = inh.rule_predict(ev)
        rows.append({"acc": r["acc"], "pos": r["pos"], "variant": r["prot"] or r["dna"], "disease": r["disease"],
                     "label": r["label"], "rule": rule["inheritance"], "rule_conf": rule["confidence"]})
        evs.append((ev, rule))
    results = {"rule": summarize("rule-based", rows, "rule")}

    if args.llm:
        extractor = pipeline.get_shared_extractor(pipeline.LOCAL_LLM_BACKEND, pipeline.LOCAL_LLM_MODEL)
        is_medgemma = getattr(extractor, "is_medgemma", True)
        prompts = [inh.build_prompt(ev, rule, is_medgemma) for ev, rule in evs]
        t0 = time.time()
        raws = inh.generate(extractor, prompts)
        print(f"[inheritance] LLM: {len(prompts)} prompts in {time.time() - t0:.0f}s", flush=True)
        for row, raw, (_, rule) in zip(rows, raws, evs):
            parsed = inh.parse(raw)
            row["llm"] = parsed["inheritance"] if parsed else rule["inheritance"]
            row["llm_conf"] = parsed["confidence"] if parsed else ""
            row["llm_parsed"] = bool(parsed)
        results["llm"] = summarize("LLM", rows, "llm")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"summary": results, "rows": rows}, indent=1, ensure_ascii=False))
    print(f"[inheritance] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
