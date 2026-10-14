#!/usr/bin/env python3
"""用 data/ 的新資料重跑 Track B 推論評估。

與先前評估的差異只在「資料」，推論邏輯完全沿用 MechanisticInference：

  類比庫   data/HGMD_latest.csv（HGMD，2,193 顆，原本 1,094 顆）
           變異類型直接取 HGMD 的 mutype 欄，不再用 HGVS 正則猜
  外顯子   cardiodb.xlsx（364 個外顯子的 Hg38 座標、band、心肌/骨骼肌 PSI）
           取代原本「轉錄本模型 + 偏移校正」與「補充表 PDF 解析」兩套來源
  獨立驗證 pubtator_ttn_variants_clinvar_annotated.csv（只用有座標的變異）

用法：conda run -n gm4 python eval_inference_data.py
"""
from __future__ import annotations

import itertools
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pandas as pd

import config
from utils.disease_normalizer import (
    canonical_labels,
    compare_disease,
    is_ttn_relevant,
    is_unspecified,
    is_vague_disease,
    normalize_tissue,
    tissue_match,
    tissue_of_labels,
)
from utils.mechanistic_inference import MechanisticInference

DATA = ROOT / "data"
HGMD_CSV = DATA / "HGMD_latest.csv"
CARDIODB = DATA / "cardiodb.xlsx"
PUBTATOR = DATA / "pubtator" / "pubtator_ttn_variants_clinvar_annotated.csv"

# HGMD mutype → (vclass, vgroup)。直接用 HGMD 自己的分類，比正則可靠。
MUTYPE_MAP = {
    "frameshift": ("frameshift", "truncating"),
    "nonsense": ("nonsense", "truncating"),
    "canonical-splice": ("splice", "truncating"),
    "exonic-splice": ("splice", "truncating"),
    "gross": ("gross", "truncating"),
    "missense": ("missense", "missense"),
    "inframe": ("inframe", "missense"),
    "noncoding": ("noncoding", "unknown"),
}


class CardioDBInference(MechanisticInference):
    """把資料來源換成 data/ 的新資料，推論邏輯不動。"""

    _cardio: Optional[pd.DataFrame] = None
    _exon_rows: Optional[List] = None

    # ── cardiodb：外顯子座標 / band / PSI ─────────────────────────────────
    def _load_cardio(self) -> List[Dict[str, Any]]:
        if self._exon_rows is not None:
            return self._exon_rows
        # 欄名含空白，itertuples 會全部改成 _1.._16，容易對錯位置；改用 dict 存取。
        df = pd.read_excel(CARDIODB).iloc[1:]  # 第 0 列是 isoform 名稱標頭
        rows = []
        for _, r in df.iterrows():
            try:
                ex = int(r["Exon number"])
                a, b = int(r["Hg38 start"]), int(r["Hg38 end"])
            except (TypeError, ValueError):
                continue
            # cardiodb 對 36 個外顯子沒有骨骼肌 PSI（exon 148/150/172-205，全在
            # I-band 且心肌 PSI 也都 <90%）。缺值視為 0，它們因此落入 both-low，
            # 語意上正確：這些是兩個組織都幾乎不收錄的替代性剪接外顯子。
            ph, ps = r["PSI in heart"], r["PSI in skeletal"]
            rows.append({
                "exon": ex,
                "lo": min(a, b), "hi": max(a, b),
                "band": str(r["Band"] or "").strip(),
                "psi_heart": 0.0 if pd.isna(ph) else float(ph),
                "psi_skeletal": 0.0 if pd.isna(ps) else float(ps),
            })
        rows.sort(key=lambda x: x["lo"])
        self._exon_rows = rows
        return rows

    def _find_exon_row(self, pos: int) -> Optional[Dict[str, Any]]:
        from bisect import bisect_right
        rows = self._load_cardio()
        starts = [r["lo"] for r in rows]
        i = bisect_right(starts, pos) - 1
        if i < 0:
            return None
        r = rows[i]
        return r if r["lo"] <= pos <= r["hi"] else None

    def assign_exon(self, pos):
        r = self._find_exon_row(int(pos))
        return r["exon"] if r else None

    def assign_band(self, pos):
        r = self._find_exon_row(int(pos))
        if r and r["band"]:
            return r["band"]
        return super().assign_band(pos)          # 內含子 → 退回座標切點

    def _load_exon_inclusion(self):
        if self._exon_inclusion is not None:
            return self._exon_inclusion
        out = {}
        for r in self._load_cardio():
            h, s = r["psi_heart"], r["psi_skeletal"]
            bucket = "psi:" + self._psi_bucket(h, s, 0)
            out[r["exon"]] = {"psi_heart": h, "psi_skeletal": s, "band": r["band"],
                              "bucket": bucket, "source": "cardiodb", "n_inc": None}
        self._exon_inclusion = out
        return out

    def hgmd_lookup(self, pos, ref=None, alt=None):
        """這份 HGMD 匯出沒有 REF/ALT 欄，無法做「同位置同鹼基」比對。

        父類別用字串相等判斷 exact，而查詢與知識庫的 ref/alt 都是空字串，
        會讓同位置的任何變異都被誤判為完全命中（實測 PubTator 有 49.2% 中招，
        等於直接把答案讀出來）。這裡明確關閉 exact，只保留同位置資訊。
        """
        hit = super().hgmd_lookup(pos, ref, alt)
        hit["exact"] = None
        return hit

    # ── HGMD 類比庫 ──────────────────────────────────────────────────────
    def _load_kb(self):
        if self._kb is not None:
            return self._kb
        df = pd.read_csv(HGMD_CSV, low_memory=False)
        rows = []
        for r in df.itertuples():
            try:
                pos = int(r.startCoord)
            except (TypeError, ValueError):
                continue
            disease = str(r.disease or "").strip()
            vclass, vgroup = MUTYPE_MAP.get(str(r.mutype or "").strip(), ("unknown", "unknown"))
            rows.append({
                "pos": pos,
                "ref": "", "alt": "",
                "dna": str(getattr(r, "hgvs", "") or ""),
                "prot": str(getattr(r, "amino", "") or ""),
                "phen": disease,
                # 引擎內部用小寫組織標籤；tissue_of_labels 回傳 "Cardiac" 這種首字大寫
                "phenotype": tissue_of_labels(canonical_labels(disease)).lower(),
                "rankscore": None if pd.isna(r.rankscore) else float(r.rankscore),
                "rsid": str(getattr(r, "dbsnp", "") or ""),
                "clinvar_id": None if pd.isna(r.clinvarID) else str(r.clinvarID),
                "phen_vague": is_vague_disease(disease),
                "band": self.assign_band(pos),
                "exon": self.assign_exon(pos),
                "vclass": vclass, "vgroup": vgroup,
            })
        self._kb = rows
        return rows


def build_engine(mode="band", hotspot=True, rank=0.0, vague=1.0, incl=2.0,
                 predict_mode="analogy", m_exon=4.0, m_bucket=4.0, m_band=8.0):
    return CardioDBInference(
        kb_path=str(HGMD_CSV), band_edges=config.TTN_BAND_EDGES,
        band_labels=config.TTN_BAND_LABELS, top_k=config.INFERENCE_TOP_K_ANALOGS,
        transcript_model_path=config.TTN_TRANSCRIPT_MODEL_PATH,
        retrieval_mode=mode, exon_hotspot_rules=hotspot,
        rankscore_weighting=rank, vague_penalty=vague,
        inclusion_tissue_weighting=incl > 0, inclusion_lr_strength=incl,
        predict_mode=predict_mode, m_exon=m_exon, m_bucket=m_bucket, m_band=m_band,
    )


# ── 評估集 ─────────────────────────────────────────────────────────────────
def hgmd_cases(mi) -> List[Dict[str, Any]]:
    out = []
    for r in mi._load_kb():
        d = r["phen"]
        if not d or is_unspecified(d):
            continue
        labels = canonical_labels(d)
        if not is_ttn_relevant(labels):
            continue
        out.append({
            "variant_id": f"{r['pos']}-{r['vclass']}",
            "pos": r["pos"], "truth_disease": d,
            "truth_tissue": normalize_tissue(r["phenotype"]),
            "variant_info": {"pos": r["pos"], "ref": "", "alt": "",
                             "hgvs_vclass": r["vclass"], "hgvs_vgroup": r["vgroup"]},
        })
    return out


def pubtator_cases(pathogenic_only=True) -> List[Dict[str, Any]]:
    """只用有座標的變異；同一位置的多筆文獻標註取眾數。"""
    d = pd.read_csv(PUBTATOR, low_memory=False)
    d = d[d.position.notna()].copy()
    if pathogenic_only:
        d = d[d.pathogenicity.astype(str).str.lower().str.startswith(("pathogenic", "likely patho"))]
    d["pos"] = d.position.astype(float).astype(int)
    out = []
    for pos, g in d.groupby("pos"):
        dis = g.disease.dropna()
        dis = dis[~dis.astype(str).str.lower().isin(["not specified", "nan", ""])]
        if dis.empty:
            continue
        disease = dis.mode().iloc[0]
        labels = canonical_labels(disease)
        if not is_ttn_relevant(labels) or is_unspecified(disease):
            continue
        tis = g.tissue_affected.dropna()
        tis = tis[~tis.astype(str).str.lower().isin(["not specified", "others", "nan"])]
        truth_tissue = normalize_tissue(tis.mode().iloc[0]) if not tis.empty else tissue_of_labels(labels)
        out.append({
            "variant_id": str(pos), "pos": pos,
            "truth_disease": disease, "truth_tissue": truth_tissue,
            "variant_info": {"pos": pos, "ref": "", "alt": ""},
        })
    return out


HIT = ("exact", "compatible")


def run(mi, cases, loo: bool) -> pd.DataFrame:
    full = list(mi._load_kb())
    rows = []
    for c in cases:
        if loo:
            mi._kb = [r for r in full if r["pos"] != c["pos"]]
        try:
            f = mi.build_features(c["variant_info"])
            inf = mi._rule_based_fallback(f)
        except Exception:
            continue
        cmp_ = compare_disease(inf.get("inferred_disease"), c["truth_disease"],
                               inf.get("differential_diagnoses") or [])
        rows.append({
            "disease_eval": cmp_["level"], "pred_vague": cmp_["pred_vague"],
            "truth_labels": tuple(cmp_["truth_labels"]),
            "tissue_eval": tissue_match(normalize_tissue(inf.get("inferred_tissue")),
                                        c["truth_tissue"]),
            "band": f.get("band"),
        })
    mi._kb = full
    return pd.DataFrame(rows)


def metrics(df: pd.DataFrame) -> Dict[str, float]:
    if df.empty:
        return {}
    hit = df.disease_eval.isin(HIT)
    vague = df.pred_vague.fillna(False).astype(bool)
    per = [g.disease_eval.isin(HIT).mean() for _, g in df.groupby("truth_labels") if len(g) >= 3]
    t = df.tissue_eval
    scored = (t != "n/a").sum()
    return {"n": len(df), "hit": hit.mean(), "specific": (hit & ~vague).mean(),
            "macro": sum(per) / len(per) if per else float("nan"),
            "tissue": (t.isin(["exact", "partial"]).sum() / scored) if scored else float("nan")}


if __name__ == "__main__":
    probe = build_engine()
    kb = probe._load_kb()
    incl = probe._load_exon_inclusion()
    print(f"類比庫：{len(kb)} 顆（原 1,094）")
    print("  變異類型：", pd.Series([r["vgroup"] for r in kb]).value_counts().to_dict())
    print("  組織標註：", pd.Series([r["phenotype"] for r in kb]).value_counts().to_dict())
    print("  band    ：", pd.Series([r["band"] for r in kb]).value_counts().to_dict())
    print(f"外顯子表：{len(incl)} 個")
    print("  PSI 分桶：", pd.Series([v["bucket"] for v in incl.values()]).value_counts().to_dict())

    hg = hgmd_cases(probe)
    pt = pubtator_cases()
    print(f"\n評估集：HGMD leave-one-out n={len(hg)}   PubTator（有座標，P/LP）n={len(pt)}\n")

    CONFIGS = [
        ("baseline",  dict(mode="band", hotspot=False, rank=0.0, vague=1.0, incl=0.0)),
        ("+hotspot",  dict(mode="band", hotspot=True,  rank=0.0, vague=1.0, incl=0.0)),
        ("+exon",     dict(mode="exon", hotspot=True,  rank=0.0, vague=1.0, incl=0.0)),
        ("+psi",      dict(mode="exon", hotspot=True,  rank=0.5, vague=0.6, incl=2.0)),
        ("band-hit",  dict(mode="band", hotspot=True,  rank=0.0, vague=1.0, incl=2.0)),
    ]
    def naive(cases):
        """多數類別基準：疾病一律答 DCM、組織一律答 Cardiac。"""
        rows = []
        for c in cases:
            cm = compare_disease("Dilated cardiomyopathy", c["truth_disease"], [])
            rows.append({"disease_eval": cm["level"], "pred_vague": cm["pred_vague"],
                         "truth_labels": tuple(cm["truth_labels"]),
                         "tissue_eval": tissue_match("Cardiac", c["truth_tissue"]), "band": None})
        return metrics(pd.DataFrame(rows))

    print(f"{'組態':<12}{'HGMD hit':>10}{'specific':>10}{'macro':>9}{'tissue':>9}"
          f"{'PT hit':>9}{'specific':>10}{'tissue':>9}")
    na, nb = naive(hg), naive(pt)
    print(f"{'多數類別':<10}{na['hit']:>10.1%}{na['specific']:>10.1%}{na['macro']:>9.1%}"
          f"{na['tissue']:>9.1%}{nb['hit']:>9.1%}{nb['specific']:>10.1%}{nb['tissue']:>9.1%}"
          "   ← 一律答 DCM/Cardiac")
    for name, kw in CONFIGS:
        e = build_engine(**kw)
        a, b = metrics(run(e, hg, True)), metrics(run(e, pt, False))
        print(f"{name:<12}{a['hit']:>10.1%}{a['specific']:>10.1%}{a['macro']:>9.1%}"
              f"{a['tissue']:>9.1%}{b['hit']:>9.1%}{b['specific']:>10.1%}{b['tissue']:>9.1%}")

    out = ROOT / "outputs" / "eval_inference_data.csv"
    out.parent.mkdir(exist_ok=True)
    recs = [{"config": "majority-class", **{f"hgmd_{k}": v for k, v in na.items()},
             **{f"pt_{k}": v for k, v in nb.items()}}]
    for name, kw in CONFIGS:
        e = build_engine(**kw)
        a, b = metrics(run(e, hg, True)), metrics(run(e, pt, False))
        recs.append({"config": name, **{f"hgmd_{k}": v for k, v in a.items()},
                     **{f"pt_{k}": v for k, v in b.items()}})
    pd.DataFrame(recs).to_csv(out, index=False)
    print(f"\n明細寫入 {out}")

    # PSI 組織 LR 的一致性檢查
    print("\nPSI 組織 LR 消融（HGMD hit）：")
    win = 0
    for mode, rank, vague in itertools.product(("band", "exon"), (0.0, 0.5), (1.0, 0.6)):
        hs = []
        for i in (0.0, 2.0):
            e = build_engine(mode=mode, rank=rank, vague=vague, incl=i)
            hs.append(metrics(run(e, hg, True))["hit"])
        win += hs[1] > hs[0]
        print(f"  {mode}/rank={rank:g}/vague={vague:g}: {hs[0]:.1%} → {hs[1]:.1%}  ({hs[1]-hs[0]:+.1%})")
    print(f"  改善格數：{win}/8")
