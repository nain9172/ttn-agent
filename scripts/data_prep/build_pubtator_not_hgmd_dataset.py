#!/usr/bin/env python3
"""
建立「HGMD 沒有、但 PubTator 有」的變異資料集（加強版）。

用途：抽出可用來測試 Track B 機制推理的變異——這些變異在 HGMD 類比 KB 中不存在
（HGMD 沒有），但 PubTator 有文獻標註（可當 ground truth）。

相較初版，這一版在同一次 ClinVar VCV XML 抓取中，一併取得：
  - GRCh38(chr2) 座標（positionVCF / ref / alt）
  - MolecularConsequence（missense / nonsense / frameshift / splice ...）→ 直接判類型，
    不再因為只有 c.HGVS 而落在 unknown
  - 代表性蛋白質 HGVS（優先 NP_596869.4）

兩條變異來源路徑：
  1. 有 ClinVar_Variation_ID：直接用 ID 抓 VCV。
  2. 無 ClinVar_Variation_ID：用 variant_text / hgvs / rsid 反查 ClinVar 取得 ID 再抓 VCV
     （可用 --no-hgvs-lookup 關閉）。

快取：
  - outputs/pubtator_vcv_cache.json：cid -> {pos,ref,alt,consequence,prot}
  - outputs/pubtator_noid_clinvar_cache.json：hgvs 反查 -> cid
皆支援中斷續跑。

用法：
    conda run -n gm4 python build_pubtator_not_hgmd_dataset.py
    conda run -n gm4 python build_pubtator_not_hgmd_dataset.py --limit 50
    conda run -n gm4 python build_pubtator_not_hgmd_dataset.py --no-hgvs-lookup
"""

from __future__ import annotations

import argparse
import json
import re
import time
from bisect import bisect_left
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

import sys as _sys
from pathlib import Path as _Path

# 本檔不在 repo 根目錄：把根目錄（config / utils）與本檔所在資料夾放進 import 路徑
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
_sys.path.insert(1, str(_Path(__file__).resolve().parent))

import config
import utils.clinvar_parser  # noqa: F401  匯入即設定 Entrez.email / api_key
from clinvar_variation_id import ClinVarVariationLookup
from utils.mechanistic_inference import MechanisticInference, _classify_variant

try:
    from Bio import Entrez
except ImportError as e:  # pragma: no cover
    raise SystemExit("需要 Biopython（gm4 環境）") from e

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PUBTATOR_CSV = PROJECT_ROOT / "data" / "pubtator" / "pubtator_ttn_variants_all_with_clinvar.csv"
VCV_CACHE = PROJECT_ROOT / "outputs" / "pubtator_vcv_cache.json"
NOID_CACHE = PROJECT_ROOT / "outputs" / "pubtator_noid_clinvar_cache.json"
OUT_CSV = PROJECT_ROOT / "outputs" / "pubtator_not_in_hgmd_dataset.csv"
EUTILS_SLEEP = 0.15


# ── 快取 ─────────────────────────────────────────────────────────────────────

def _load_json(path: Path) -> Dict:
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def _save_json(path: Path, data: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)


# ── VCV XML 解析：座標 + consequence + 蛋白質 HGVS ──────────────────────────

# MolecularConsequence Type → (vclass, vgroup)，取最嚴重者
_CONSEQUENCE_RANK = [
    ("nonsense", ("nonsense", "truncating")),
    ("stop gain", ("nonsense", "truncating")),
    ("frameshift", ("frameshift", "truncating")),
    ("splice donor", ("splice", "truncating")),
    ("splice acceptor", ("splice", "truncating")),
    ("missense", ("missense", "missense")),
    ("inframe", ("inframe", "other")),
    ("synonymous", ("synonymous", "other")),
    ("splice region", ("splice_region", "other")),
    ("intron", ("intron", "other")),
    ("5 prime utr", ("utr", "other")),
    ("3 prime utr", ("utr", "other")),
]


def _map_consequence(mc_types: List[str]) -> Dict[str, str]:
    low = [t.lower() for t in mc_types]
    for key, (vclass, vgroup) in _CONSEQUENCE_RANK:
        if any(key in t for t in low):
            return {"vclass": vclass, "vgroup": vgroup}
    return {"vclass": "unknown", "vgroup": "unknown"}


def parse_vcv(cid: str, cache: Dict) -> Optional[Dict]:
    """cid -> {pos,ref,alt,consequence,prot}；查不到回 None。結果快取。"""
    if cid in cache:
        return cache[cid]
    out: Optional[Dict] = None
    try:
        time.sleep(EUTILS_SLEEP)
        acc = f"VCV{int(cid):09d}"
        handle = Entrez.efetch(db="clinvar", id=acc, rettype="vcv", retmode="xml")
        xml = handle.read()
        handle.close()
        if isinstance(xml, bytes):
            xml = xml.decode("utf-8", errors="replace")

        pos = ref = alt = None
        for m in re.finditer(r"<SequenceLocation\s+([^/>]+)/?>", xml):
            attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
            if attrs.get("Assembly") == "GRCh38" and attrs.get("Chr") == "2":
                pv = attrs.get("positionVCF") or attrs.get("start")
                rf = (attrs.get("referenceAlleleVCF") or "").upper()
                al = (attrs.get("alternateAlleleVCF") or "").upper()
                if pv and rf and al:
                    pos, ref, alt = int(pv), rf, al
                    break
        if pos is not None:
            mc = list(dict.fromkeys(re.findall(r'<MolecularConsequence[^>]*\bType="([^"]+)"', xml)))
            # 蛋白質 HGVS：優先 NP_596869.4（HGMD 用的 meta transcript）
            prot = None
            m2 = re.search(r"NP_596869\.\d+:(p\.[A-Za-z0-9*=_]+)", xml)
            if m2:
                prot = m2.group(1)
            else:
                m3 = re.search(r":(p\.[A-Za-z0-9*=_]+)", xml)
                prot = m3.group(1) if m3 else None
            out = {"pos": pos, "ref": ref, "alt": alt, "consequence": mc, "prot": prot}
    except Exception:
        out = None
    cache[cid] = out
    _save_json(VCV_CACHE, cache)
    return out


# ── 聚合工具 ─────────────────────────────────────────────────────────────────

def _mode(values: List) -> Optional[str]:
    clean = [str(v).strip() for v in values
             if pd.notna(v) and str(v).strip() and str(v).strip().lower() != "nan"]
    return Counter(clean).most_common(1)[0][0] if clean else None


def _uniq_join(values: List, sep: str = "; ") -> str:
    out: List[str] = []
    for v in values:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if s and s.lower() != "nan" and s not in out:
            out.append(s)
    return sep.join(out)


def _decide_type(vcv: Dict, hgvs_rep: Optional[str]) -> Dict[str, str]:
    """優先用 ClinVar MolecularConsequence；沒有時退回 HGVS 解析。"""
    mc = vcv.get("consequence") or []
    typ = _map_consequence(mc)
    if typ["vgroup"] != "unknown":
        return typ
    return _classify_variant(vcv.get("prot") or "", hgvs_rep or "")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="只處理前 N 個 ClinVar ID（0=全部）")
    ap.add_argument("--no-hgvs-lookup", action="store_true", help="不處理無 ClinVar ID 的變異")
    ap.add_argument("--max-noid", type=int, default=0, help="無 ID 路徑最多反查幾個 HGVS（0=全部）")
    args = ap.parse_args()

    mi = MechanisticInference(
        kb_path=config.HGMD_KB_PATH,
        band_edges=config.TTN_BAND_EDGES,
        band_labels=config.TTN_BAND_LABELS,
        top_k=config.INFERENCE_TOP_K_ANALOGS,
        transcript_interval_path=config.TRANSCRIPT_INTERVAL_PATH,
    )
    kb = mi._load_kb()
    hgmd_pos_sorted = sorted({r["pos"] for r in kb})
    hgmd_exact = {(r["pos"], str(r["ref"]).upper(), str(r["alt"]).upper()) for r in kb}
    pos_set = set(hgmd_pos_sorted)

    def nearest_hgmd_dist(pos: int) -> int:
        i = bisect_left(hgmd_pos_sorted, pos)
        best = 10 ** 12
        for j in (i - 1, i):
            if 0 <= j < len(hgmd_pos_sorted):
                best = min(best, abs(hgmd_pos_sorted[j] - pos))
        return best

    p = pd.read_csv(PUBTATOR_CSV)
    vcv_cache = _load_json(VCV_CACHE)

    rows: List[Dict] = []
    seen_var: set = set()
    n_unresolved = 0
    n_exact_in_hgmd = 0

    def add_variant(vcv: Dict, g: pd.DataFrame, source: str) -> None:
        nonlocal n_exact_in_hgmd
        pos, ref, alt = vcv["pos"], vcv["ref"], vcv["alt"]
        key = (pos, ref, alt)
        if key in seen_var:
            return
        if key in hgmd_exact:
            n_exact_in_hgmd += 1
            seen_var.add(key)
            return
        seen_var.add(key)
        hgvs_rep = _mode(g["hgvs"].tolist())
        typ = _decide_type(vcv, hgvs_rep)
        rows.append({
            "variant_id": f"2-{pos}-{ref}-{alt}",
            "chrom": "2",
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "is_snv": len(ref) == 1 and len(alt) == 1,
            "band": mi.assign_band(pos),
            "vclass": typ["vclass"],
            "vgroup": typ["vgroup"],
            "molecular_consequence": "; ".join(vcv.get("consequence") or []),
            "protein_hgvs": vcv.get("prot"),
            "clinvar_id": _mode(g["ClinVar_Variation_ID"].tolist()) if "ClinVar_Variation_ID" in g else None,
            "rsid": _mode(g["rsid"].tolist()),
            "hgvs": hgvs_rep,
            "pubtator_pathogenicity": _mode(g["pathogenicity"].tolist()),
            "pubtator_disease": _mode(g["disease"].tolist()),
            "pubtator_tissue": _mode(g["tissue_affected"].tolist()),
            "pubtator_inheritance": _mode(g["inheritance"].tolist()),
            "n_pubtator_rows": len(g),
            "pmids": _uniq_join(g["pmid"].tolist()),
            "years": _uniq_join(g["year"].tolist()),
            "source": source,
            "hgmd_exact": False,
            "hgmd_same_position": pos in pos_set,
            "nearest_hgmd_dist_bp": nearest_hgmd_dist(pos),
            "pure_inference": pos not in pos_set,
        })

    # ── Path 1：有 ClinVar ID ────────────────────────────────────────────────
    has_cid = p[p["ClinVar_Variation_ID"].notna()].copy()
    has_cid["cid"] = has_cid["ClinVar_Variation_ID"].astype(int).astype(str)
    grouped = has_cid.groupby("cid")
    cids = list(grouped.groups.keys())
    if args.limit:
        cids = cids[: args.limit]
    print(f"[Path1] unique ClinVar IDs: {len(cids)}")
    for k, cid in enumerate(cids, 1):
        vcv = parse_vcv(cid, vcv_cache)
        if not vcv:
            n_unresolved += 1
        else:
            add_variant(vcv, grouped.get_group(cid), source="clinvar_id")
        if k % 100 == 0:
            print(f"  [Path1] {k}/{len(cids)} | collected {len(rows)}")

    # ── Path 2：無 ClinVar ID，用 HGVS/rsid 反查 ClinVar ─────────────────────
    if not args.no_hgvs_lookup:
        no_cid = p[p["ClinVar_Variation_ID"].isna()].copy()
        no_cid["hgvs_key"] = no_cid["hgvs"].astype(str).str.replace(r"\s+", "", regex=True)
        noid_groups = no_cid.groupby("hgvs_key")
        keys = [k for k in noid_groups.groups.keys() if k and k.lower() != "nan"]
        if args.max_noid:
            keys = keys[: args.max_noid]
        print(f"[Path2] unique HGVS without ClinVar ID: {len(keys)}")
        lookup = ClinVarVariationLookup()
        noid_cache = _load_json(NOID_CACHE)
        for k, hk in enumerate(keys, 1):
            g = noid_groups.get_group(hk)
            if hk in noid_cache:
                cid = noid_cache[hk]
            else:
                r0 = g.iloc[0]
                cid = lookup.lookup_pubtator(
                    r0.get("variant_text"), r0.get("hgvs"), r0.get("rsid")
                )
                noid_cache[hk] = cid
                _save_json(NOID_CACHE, noid_cache)
            if not cid:
                continue
            vcv = parse_vcv(str(cid), vcv_cache)
            if not vcv:
                n_unresolved += 1
                continue
            add_variant(vcv, g, source="hgvs_lookup")
            if k % 100 == 0:
                print(f"  [Path2] {k}/{len(keys)} | collected {len(rows)}")

    if not rows:
        print("沒有可輸出的變異。")
        return

    df = pd.DataFrame(rows).sort_values(["pos"]).reset_index(drop=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    print("\n" + "=" * 64)
    print(f"HGMD 沒有、PubTator 有的變異：{len(df)}")
    print(f"  來源：clinvar_id={int((df['source']=='clinvar_id').sum())} | "
          f"hgvs_lookup={int((df['source']=='hgvs_lookup').sum())}")
    print(f"  SNV：{int(df['is_snv'].sum())} | indel/其他：{int((~df['is_snv']).sum())}")
    print(f"  pure_inference（連同位置都沒有）：{int(df['pure_inference'].sum())} | "
          f"同位置有 HGMD：{int(df['hgmd_same_position'].sum())}")
    print(f"  已知類型（非 unknown）：{int((df['vgroup'] != 'unknown').sum())}/{len(df)}")
    print(f"  vgroup 分佈：{df['vgroup'].value_counts().to_dict()}")
    print(f"  band 分佈：{df['band'].value_counts().to_dict()}")
    print(f"\n略過：exact 命中 HGMD {n_exact_in_hgmd}；座標查不到 {n_unresolved}")
    print(f"輸出：{OUT_CSV}")


if __name__ == "__main__":
    main()
