#!/usr/bin/env python3
"""把 pubtator-only 驗證結果的統計數字整理成 JSON（供報表 / canvas 使用）。

用法：
    conda run -n gm4 python summarize_validation_stats.py
    conda run -n gm4 python summarize_validation_stats.py --json outputs/validation_stats.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

HIT = ("exact", "compatible")
DISEASE_ORDER = ["exact", "compatible", "differential", "same_tissue", "miss", "n/a"]


def _vague(df: pd.DataFrame) -> pd.Series:
    if "pred_vague" not in df.columns:
        return pd.Series(False, index=df.index)
    return df["pred_vague"].fillna(False).astype(bool)


def block(df: pd.DataFrame) -> Dict:
    """一個子集的核心指標。"""
    n = len(df)
    if not n:
        return {}
    d = df["disease_eval"].fillna("n/a")
    hit = int(d.isin(HIT).sum())
    specific = int((d.isin(HIT) & ~_vague(df)).sum())
    with_diff = hit + int((d == "differential").sum())
    t = df["tissue_eval"].fillna("n/a")
    t_scored = int((t != "n/a").sum())
    t_ok = int(t.isin(("exact", "partial")).sum())
    return {
        "n": n,
        "exact": int((d == "exact").sum()),
        "hit": hit,
        "hit_pct": round(100 * hit / n, 1),
        "specific": specific,
        "specific_pct": round(100 * specific / n, 1),
        "with_diff": with_diff,
        "with_diff_pct": round(100 * with_diff / n, 1),
        "tissue_ok": t_ok,
        "tissue_scored": t_scored,
        "tissue_pct": round(100 * t_ok / t_scored, 1) if t_scored else None,
        "dist": {k: int((d == k).sum()) for k in DISEASE_ORDER if int((d == k).sum())},
    }


def by(df: pd.DataFrame, col: str, min_n: int = 1) -> List[Dict]:
    """依欄位分組的命中率（含具體命中，用來看籠統答案的比重）。"""
    out = []
    for key, g in df.groupby(df[col].fillna("(未知)")):
        if len(g) < min_n:
            continue
        d = g["disease_eval"].fillna("n/a")
        hit = int(d.isin(HIT).sum())
        specific = int((d.isin(HIT) & ~_vague(g)).sum())
        out.append({
            "key": str(key),
            "n": len(g),
            "hit": hit,
            "hit_pct": round(100 * hit / len(g), 1),
            "specific": specific,
            "specific_pct": round(100 * specific / len(g), 1),
        })
    return sorted(out, key=lambda r: -r["n"])


def confusion(df: pd.DataFrame, top: int = 12) -> List[Dict]:
    g = (df.assign(t=df["truth_labels"].fillna("?"), p=df["pred_labels"].fillna("?"))
           .groupby(["t", "p", "disease_eval"]).size()
           .sort_values(ascending=False).head(top))
    return [{"truth": t, "pred": p, "level": lvl, "n": int(n)}
            for (t, p, lvl), n in g.items()]


def confidence_calibration(df: pd.DataFrame) -> List[Dict]:
    """LLM 自報信心 vs 實際命中率——看它知不知道自己什麼時候不確定。"""
    if "pred_confidence" not in df.columns:
        return []
    out = []
    for key in ("High", "Moderate", "Low"):
        g = df[df["pred_confidence"] == key]
        if not len(g):
            continue
        d = g["disease_eval"].fillna("n/a")
        hit = int(d.isin(HIT).sum())
        out.append({"key": key, "n": len(g), "hit": hit,
                    "hit_pct": round(100 * hit / len(g), 1)})
    return out


def summarize(path: Path, label: str) -> Optional[Dict]:
    if not path.exists():
        return None
    df = pd.read_csv(path)
    pure = df[df["pure_inference"].fillna(False).astype(bool)]
    return {
        "label": label,
        "file": path.name,
        "overall": block(df),
        "pure": block(pure),
        "hgmd_same_pos": block(df[~df["pure_inference"].fillna(False).astype(bool)]),
        "by_band": by(df, "band"),
        "by_vclass": by(df, "vclass"),
        "by_vgroup": by(df, "vgroup"),
        "by_resolver": by(df.assign(src=df["resolved_by"].str.split(":").str[0]), "src"),
        "by_analog_distance": by(df.assign(bucket=pd.cut(
            df["nearest_analog_bp"].fillna(10 ** 9),
            bins=[-1, 0, 100, 1000, 10000, 10 ** 12],
            labels=["同位置", "≤100bp", "≤1kb", "≤10kb", ">10kb"])), "bucket"),
        "confidence": confidence_calibration(df),
        "confusion": confusion(df),
        "truth_top": [{"key": str(k), "n": int(v)} for k, v in
                      df["truth_labels"].fillna("?").value_counts().head(8).items()],
        "pred_top": [{"key": str(k), "n": int(v)} for k, v in
                     df["pred_labels"].fillna("?").value_counts().head(8).items()],
    }


def show(s: Dict) -> None:
    o = s["overall"]
    print(f"\n{'=' * 70}\n{s['label']}（{s['file']}）  n={o['n']}")
    print(f"  疾病 exact           {o['exact']:>4}/{o['n']} ({100 * o['exact'] / o['n']:.1f}%)")
    print(f"  疾病命中 exact+compat {o['hit']:>4}/{o['n']} ({o['hit_pct']}%)")
    print(f"  其中具體（非籠統詞）  {o['specific']:>4}/{o['n']} ({o['specific_pct']}%)")
    print(f"  含 differential       {o['with_diff']:>4}/{o['n']} ({o['with_diff_pct']}%)")
    print(f"  組織相容              {o['tissue_ok']:>4}/{o['tissue_scored']} ({o['tissue_pct']}%)")
    print(f"  純推理子集 n={s['pure']['n']}：命中 {s['pure']['hit_pct']}%、"
          f"具體 {s['pure']['specific_pct']}%、組織 {s['pure']['tissue_pct']}%")
    for name, key in (("band", "by_band"), ("變異類型", "by_vclass"),
                      ("類比距離", "by_analog_distance"), ("信心", "confidence")):
        if not s.get(key):
            continue
        print(f"  分{name}：")
        for r in s[key]:
            extra = f"  具體 {r['specific_pct']}%" if "specific_pct" in r else ""
            print(f"    {r['key']:<10} {r['hit']:>4}/{r['n']:<5} {r['hit_pct']:>5}%{extra}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rule", type=Path,
                    default=Path("outputs/inference_pubtator_only_validation.csv"))
    ap.add_argument("--llm", type=Path,
                    default=Path("outputs/inference_pubtator_only_validation_llm.csv"))
    ap.add_argument("--json", type=Path, default=Path("outputs/validation_stats.json"))
    args = ap.parse_args()

    stats = [s for s in (summarize(args.llm, "LLM 類比推理"),
                         summarize(args.rule, "規則式基線")) if s]
    for s in stats:
        show(s)

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(stats, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON 已寫入：{args.json}")


if __name__ == "__main__":
    main()
