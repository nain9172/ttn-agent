#!/usr/bin/env python3
"""
用「HGMD 沒有、但 PubTator 有」的變異驗證 Track B 機制推理（leave-out validation）。

想法：
  - PubTator 提供文獻 ground truth（disease / tissue_affected）。
  - 只挑「不在 HGMD 類比 KB」的變異（連同位置都沒有）→ 對機制推理是公平的
    held-out 測試：推理拿不到 exact / same-position 捷徑，純靠 band + 類型 +
    isoform + 距離加權的 HGMD 類比。
  - 比對推理輸出（disease / tissue）與 PubTator 標註，估 tissue 準確率。

流程：
  1. 由 HGMD_TTN_processed 建立「已知位置」集合（同位置=非純推理，排除）。
  2. PubTator 逐筆用 ClinVar_Variation_ID → efetch VCV XML 取 GRCh38 座標
     （快取到 outputs/pubtator_genomic_cache.json，可續跑）。
  3. 篩出 SNV 且位置不在 HGMD → 純推理測試集。
  4. 對每筆跑 MechanisticInference（預設 rule-based；--llm 走真正模型）。
  5. 輸出比對表 CSV + tissue 準確率摘要。

用法：
    conda run -n gm4 python test_inference_pubtator.py                 # 預設 rule-based，取 20 筆
    conda run -n gm4 python test_inference_pubtator.py --target 30
    conda run -n gm4 python test_inference_pubtator.py --llm           # 用實際 Gemma 推理（需 GPU）
"""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

import config
import utils.clinvar_parser  # noqa: F401  匯入即設定 Entrez.email / api_key
from utils.mechanistic_inference import MechanisticInference

try:
    from Bio import Entrez
except ImportError as e:  # pragma: no cover
    raise SystemExit("需要 Biopython（gm4 環境）") from e

PROJECT_ROOT = Path(__file__).parent
PUBTATOR_CSV = PROJECT_ROOT / "pubtator_ttn_variants_all_with_clinvar.csv"
CACHE_PATH = PROJECT_ROOT / "outputs" / "pubtator_genomic_cache.json"
OUT_CSV = PROJECT_ROOT / "outputs" / "inference_pubtator_validation.csv"
EUTILS_SLEEP = 0.2


# ── ClinVar Variation ID → GRCh38 (chr2) genomic 座標 ────────────────────────

def _load_cache() -> Dict[str, Optional[list]]:
    if CACHE_PATH.exists():
        try:
            return json.loads(CACHE_PATH.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def _save_cache(cache: Dict) -> None:
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp = CACHE_PATH.with_suffix(".tmp")
    tmp.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")
    tmp.replace(CACHE_PATH)


def resolve_genomic(cid: str, cache: Dict) -> Optional[Tuple[int, str, str]]:
    """回傳 (pos, ref, alt)（GRCh38, chr2 正鏈）；查不到回 None。"""
    if cid in cache:
        v = cache[cid]
        return tuple(v) if v else None  # type: ignore[return-value]
    result: Optional[Tuple[int, str, str]] = None
    try:
        time.sleep(EUTILS_SLEEP)
        acc = f"VCV{int(cid):09d}"
        handle = Entrez.efetch(db="clinvar", id=acc, rettype="vcv", retmode="xml")
        xml = handle.read()
        handle.close()
        if isinstance(xml, bytes):
            xml = xml.decode("utf-8", errors="replace")
        for m in re.finditer(r"<SequenceLocation\s+([^/>]+)/?>", xml):
            attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
            if attrs.get("Assembly") == "GRCh38" and attrs.get("Chr") == "2":
                pos = attrs.get("positionVCF") or attrs.get("start")
                ref = (attrs.get("referenceAlleleVCF") or "").upper()
                alt = (attrs.get("alternateAlleleVCF") or "").upper()
                if pos and ref and alt:
                    result = (int(pos), ref, alt)
                    break
    except Exception:
        result = None
    cache[cid] = list(result) if result else None
    _save_cache(cache)
    return result


# ── tissue 正規化 / 比對 ─────────────────────────────────────────────────────

def norm_tissue(t: Optional[str]) -> str:
    if not t:
        return "Uncertain"
    s = str(t).strip().lower()
    if "both" in s:
        return "Both"
    if "cardiac" in s or "heart" in s or "dcm" in s or "cardiomyopath" in s:
        return "Cardiac"
    if "skeletal" in s or "muscle" in s or "myopath" in s or "dystrophy" in s:
        return "Skeletal"
    return "Uncertain"


def tissue_hit(pred: str, truth: str) -> str:
    """Both 視為與任一單一組織相容（保守計分）。"""
    if truth == "Uncertain" or pred == "Uncertain":
        return "n/a"
    if pred == truth:
        return "exact"
    if "Both" in (pred, truth):
        return "partial"
    return "miss"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=int, default=20, help="要收集的純推理測試變異數")
    ap.add_argument("--max-fetch", type=int, default=250, help="最多向 ClinVar 解析幾筆座標")
    ap.add_argument("--llm", action="store_true", help="用實際 Gemma 推理（需 GPU/vLLM）")
    args = ap.parse_args()

    # HGMD KB 位置集合（同位置=非純推理）
    mi = MechanisticInference(
        kb_path=config.HGMD_KB_PATH,
        band_edges=config.TTN_BAND_EDGES,
        band_labels=config.TTN_BAND_LABELS,
        top_k=config.INFERENCE_TOP_K_ANALOGS,
        transcript_interval_path=config.TRANSCRIPT_INTERVAL_PATH,
    )
    kb = mi._load_kb()
    hgmd_pos = {r["pos"] for r in kb}

    p = pd.read_csv(PUBTATOR_CSV)

    def clean(v):
        return None if pd.isna(v) else str(v).strip()

    # 依 ClinVar ID 去重，優先取有明確 disease/tissue 者
    seen_cid: set = set()
    candidates: List[Dict] = []
    for r in p.itertuples():
        if pd.isna(r.ClinVar_Variation_ID):
            continue
        cid = str(int(r.ClinVar_Variation_ID))
        if cid in seen_cid:
            continue
        seen_cid.add(cid)
        candidates.append({
            "cid": cid,
            "hgvs": clean(r.hgvs),
            "disease": clean(r.disease),
            "tissue": clean(r.tissue_affected),
            "path": clean(r.pathogenicity),
            "pmid": clean(r.pmid),
        })

    cache = _load_cache()
    extractor = None
    if args.llm:
        from utils.local_clinical_extractor import LocalClinicalExtractor
        # 預設值是 Llama-3.2-3B / TP=2，要照 config 傳入；模型在 __init__ 就載入
        extractor = LocalClinicalExtractor(
            backend=config.LOCAL_LLM_BACKEND,
            model_name=config.LOCAL_LLM_MODEL,
            tensor_parallel_size=config.LOCAL_LLM_TENSOR_PARALLEL,
            max_model_len=config.LOCAL_LLM_MAX_MODEL_LEN,
            max_context_length=config.LOCAL_LLM_MAX_CONTEXT_LENGTH,
        )

    rows: List[Dict] = []
    n_fetch = 0
    seen_var: set = set()
    for c in candidates:
        if len(rows) >= args.target or n_fetch >= args.max_fetch:
            break
        g = resolve_genomic(c["cid"], cache)
        n_fetch += 1
        if not g:
            continue
        pos, ref, alt = g
        if len(ref) != 1 or len(alt) != 1:  # 只測 SNV
            continue
        if pos in hgmd_pos:  # 同位置 → 非純推理，排除
            continue
        vkey = (pos, ref, alt)
        if vkey in seen_var:
            continue
        seen_var.add(vkey)

        variant_id = f"2-{pos}-{ref}-{alt}"
        vinfo = {"variant_id": variant_id, "pos": pos, "ref": ref, "alt": alt}
        # 用 PubTator 的 c.HGVS 幫助判類型（splice/frameshift 可判；missense/nonsense 可能仍 unknown）
        clinvar_info = {"hgvs_nucleotide": [c["hgvs"]]} if c["hgvs"] else None

        res = mi.infer(extractor, vinfo, clinvar_info, evo2_result=None)
        inf = res.get("inference", {})
        feat = res.get("features", {})

        pred_tissue = norm_tissue(inf.get("inferred_tissue"))
        truth_tissue = norm_tissue(c["tissue"])
        rows.append({
            "variant": variant_id,
            "clinvar_id": c["cid"],
            "pmid": c["pmid"],
            "band": feat.get("band"),
            "vclass": feat.get("vclass"),
            "n_analogs": len(feat.get("analogs", [])),
            "pubtator_disease": c["disease"],
            "pubtator_tissue_raw": c["tissue"],
            "pubtator_tissue": truth_tissue,
            "inferred_disease": inf.get("inferred_disease"),
            "inferred_tissue": pred_tissue,
            "confidence": inf.get("confidence"),
            "source": inf.get("source"),
            "tissue_eval": tissue_hit(pred_tissue, truth_tissue),
        })
        print(
            f"[{len(rows):>2}] {variant_id} | {feat.get('band'):<6} {feat.get('vclass'):<9} "
            f"| PT: {truth_tissue:<9}/{(c['disease'] or '')[:34]:<34} "
            f"-> INF: {pred_tissue:<9}/{(str(inf.get('inferred_disease')) or '')[:28]:<28} "
            f"[{tissue_hit(pred_tissue, truth_tissue)}]"
        )

    if not rows:
        print("沒有找到可用的純推理測試變異。")
        return

    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    evals = df["tissue_eval"]
    n = len(df)
    exact = (evals == "exact").sum()
    partial = (evals == "partial").sum()
    miss = (evals == "miss").sum()
    print("\n" + "=" * 64)
    print(f"測試變異數（HGMD 完全沒有 · SNV）：{n}")
    print(f"tissue 完全一致 exact：{exact}/{n} ({exact / n:.0%})")
    print(f"tissue 相容 partial（含 Both）：{partial}/{n} ({partial / n:.0%})")
    print(f"tissue 相容合計（exact+partial）：{(exact + partial)}/{n} ({(exact + partial) / n:.0%})")
    print(f"tissue 錯誤 miss：{miss}/{n} ({miss / n:.0%})")
    print(f"\n明細已寫入：{OUT_CSV}")


if __name__ == "__main__":
    main()
