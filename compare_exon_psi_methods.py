#!/usr/bin/env python3
"""比較「主要依賴 exon / PSI」的幾種疾病預測算法，全部以 PubTator 為 ground truth。

A. analogy       目前的距離加權類比投票（對照組）
B. exon-prior    純外顯子經驗先驗：P(disease | exon)，無資料時退回 band
C. hierarchical  階層式收縮 exon → PSI 桶 → band → 全庫（經驗貝氏 / m-estimate）
D. logistic      多項式邏輯斯迴歸，特徵 = PSI_heart, PSI_skeletal, log 比值,
                 band one-hot, vgroup one-hot, 外顯子相對位置

所有方法都在同一組 PubTator case 上評分，且都移除同位置的 HGMD 變異
（PubTator 有 49.2% 與 HGMD 同位置，不移除等於讀答案）。

用法：conda run -n gm4 python compare_exon_psi_methods.py
"""
from __future__ import annotations

import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")
ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

import eval_inference_data as E
from utils.disease_normalizer import compare_disease, normalize_tissue, tissue_match

ENGINE = E.build_engine()
KB = ENGINE._load_kb()
INCL = ENGINE._load_exon_inclusion()
CASES = E.pubtator_cases()

# 只用有外顯子、有疾病、有組織的 KB 列
KBR = [r for r in KB if r.get("exon") and r["phen"] and
       r["phenotype"] in ("cardiac", "skeletal", "both")]
print(f"類比庫可用列 {len(KBR)}   PubTator case {len(CASES)}")


def score(preds):
    """preds: list of (disease, tissue)，與 CASES 同序。"""
    rows = []
    for c, (d, t) in zip(CASES, preds):
        cm = compare_disease(d, c["truth_disease"], [])
        rows.append({"disease_eval": cm["level"], "pred_vague": cm["pred_vague"],
                     "truth_labels": tuple(cm["truth_labels"]),
                     "tissue_eval": tissue_match(normalize_tissue(t), c["truth_tissue"]),
                     "band": None})
    return E.metrics(pd.DataFrame(rows))


def tally(rows, key):
    """key(row) -> 分群鍵；回傳 {鍵: Counter(disease)} 與 {鍵: Counter(tissue)}."""
    d, t = defaultdict(Counter), defaultdict(Counter)
    for r in rows:
        k = key(r)
        if k is None:
            continue
        d[k][r["phen"]] += 1
        t[k][r["phenotype"]] += 1
    return d, t


# ── B. 純外顯子經驗先驗 ─────────────────────────────────────────────────────
def method_exon_prior():
    preds = []
    for c in CASES:
        rows = [r for r in KBR if r["pos"] != c["pos"]]        # 移除同位置
        ex = ENGINE.assign_exon(c["pos"])
        bd = ENGINE.assign_band(c["pos"])
        de, te = tally(rows, lambda r: r["exon"])
        db, tb = tally(rows, lambda r: r["band"])
        if ex in de and de[ex]:
            preds.append((de[ex].most_common(1)[0][0], te[ex].most_common(1)[0][0]))
        elif bd in db and db[bd]:
            preds.append((db[bd].most_common(1)[0][0], tb[bd].most_common(1)[0][0]))
        else:
            preds.append(("Cardiomyopathy, dilated", "cardiac"))
    return preds


# ── C. 階層式收縮 exon → PSI 桶 → band → 全庫 ───────────────────────────────
def method_hierarchical(m_exon=4.0, m_bucket=8.0, m_band=16.0):
    """每一層都往上一層收縮（m-estimate）。m 越大代表越不信任該層的計數。"""
    preds = []
    for c in CASES:
        rows = [r for r in KBR if r["pos"] != c["pos"]]
        ex = ENGINE.assign_exon(c["pos"])
        bd = ENGINE.assign_band(c["pos"])
        bk = (INCL.get(ex) or {}).get("bucket")

        de, te = tally(rows, lambda r: r["exon"])
        dk, tk = tally(rows, lambda r: (INCL.get(r["exon"]) or {}).get("bucket"))
        db, tb = tally(rows, lambda r: r["band"])
        g_d, g_t = Counter(), Counter()
        for r in rows:
            g_d[r["phen"]] += 1
            g_t[r["phenotype"]] += 1

        def hier(cnt_ex, cnt_bk, cnt_bd, cnt_g, labels):
            n_g = sum(cnt_g.values()) or 1
            out = {}
            for lab in labels:
                p = cnt_g[lab] / n_g                                     # 全庫
                n = sum(cnt_bd.values())
                p = (cnt_bd[lab] + m_band * p) / (n + m_band)            # band
                n = sum(cnt_bk.values())
                p = (cnt_bk[lab] + m_bucket * p) / (n + m_bucket)        # PSI 桶
                n = sum(cnt_ex.values())
                p = (cnt_ex[lab] + m_exon * p) / (n + m_exon)            # 外顯子
                out[lab] = p
            return max(out, key=out.get)

        dis = hier(de.get(ex, Counter()), dk.get(bk, Counter()),
                   db.get(bd, Counter()), g_d, set(g_d))
        tis = hier(te.get(ex, Counter()), tk.get(bk, Counter()),
                   tb.get(bd, Counter()), g_t, set(g_t))
        preds.append((dis, tis))
    return preds


# ── D. 多項式邏輯斯迴歸 ─────────────────────────────────────────────────────
def method_logistic():
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    BANDS = ["Z-disk", "I-band", "A-band", "M-band"]
    VG = ["truncating", "missense", "unknown"]

    def feat(pos, vgroup):
        ex = ENGINE.assign_exon(pos)
        info = INCL.get(ex) or {}
        h = info.get("psi_heart", 0.0) or 0.0
        s = info.get("psi_skeletal", 0.0) or 0.0
        bd = ENGINE.assign_band(pos)
        return ([h / 100, s / 100, np.log((s + 1) / (h + 1)), (ex or 0) / 364]
                + [1.0 if bd == b else 0.0 for b in BANDS]
                + [1.0 if vgroup == v else 0.0 for v in VG])

    Xtr = np.array([feat(r["pos"], r["vgroup"]) for r in KBR])
    yd = np.array([r["phen"] for r in KBR])
    yt = np.array([r["phenotype"] for r in KBR])
    pos_tr = np.array([r["pos"] for r in KBR])

    # 疾病類別太細會學不動，只保留出現 >=10 次的
    keep = {d for d, n in Counter(yd).items() if n >= 10}
    mask = np.array([d in keep for d in yd])

    sc = StandardScaler().fit(Xtr)
    preds = []
    for c in CASES:
        keep_row = pos_tr != c["pos"]                       # 移除同位置
        md = LogisticRegression(max_iter=800, C=1.0)
        mt = LogisticRegression(max_iter=800, C=1.0)
        md.fit(sc.transform(Xtr[mask & keep_row]), yd[mask & keep_row])
        mt.fit(sc.transform(Xtr[keep_row]), yt[keep_row])
        x = sc.transform([feat(c["pos"], "truncating")])
        preds.append((md.predict(x)[0], mt.predict(x)[0]))
    return preds


if __name__ == "__main__":
    results = {}

    # A. 現行類比投票
    df = E.run(E.build_engine(mode="band", hotspot=True, rank=0.0, vague=1.0, incl=2.0),
               CASES, True)
    results["A. analogy（現行）"] = E.metrics(df)

    # 多數類別基準
    results["0. 多數類別基準"] = score([("Dilated cardiomyopathy", "Cardiac")] * len(CASES))

    print("\n計算 B. exon-prior ...");    results["B. exon-prior"] = score(method_exon_prior())
    print("計算 C. hierarchical ...");    results["C. hierarchical"] = score(method_hierarchical())
    print("計算 D. logistic ...（較慢）"); results["D. logistic"] = score(method_logistic())

    print(f"\n{'方法':<22}{'hit':>8}{'specific':>10}{'macro':>9}{'tissue':>9}")
    for k in sorted(results):
        m = results[k]
        print(f"{k:<22}{m['hit']:>8.1%}{m['specific']:>10.1%}{m['macro']:>9.1%}{m['tissue']:>9.1%}")

    out = ROOT / "outputs" / "exon_psi_methods.csv"
    pd.DataFrame(results).T.to_csv(out)
    print(f"\n寫入 {out}")
