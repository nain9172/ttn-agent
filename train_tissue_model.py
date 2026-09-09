#!/usr/bin/env python3
"""
訓練並驗證組織（心肌 / 骨骼肌）判定器。

背景：benchmark_inference.py 顯示組織判定是 Track B 的閘門——組織判對時疾病命中率
69.8%，判錯時只剩 10.0%。現行邏輯讓疾病投票順帶決定組織，這裡改成獨立預測。

要打敗的基準線是現行規則式組織判定：
    HGMD leave-one-out  73.1%
    PubTator 跨來源      59.1%

驗證設計（重點是不能自我欺騙）：
  - HGMD 的類比特徵一律 leave-one-out（連同位置的變異一起移除），否則答案會直接
    從自己身上洩漏進 analog_cardiac_frac。
  - 除了同源的交叉驗證，一定要看跨來源結果（HGMD 訓練 → PubTator 測試，以及反向）。
    同源分數再高，換一個標註來源就垮掉的話沒有意義。
  - 兩個模型（logistic regression / gradient boosting）一起比，避免把單一模型的
    配適雜訊當成結論。

用法：
    conda run -n gm4 python train_tissue_model.py                 # 只評估，不存檔
    conda run -n gm4 python train_tissue_model.py --save          # 評估後訓練最終模型並存檔
    conda run -n gm4 python train_tissue_model.py --model logreg
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

import config
from benchmark_inference import build_engine, build_hgmd_cases, build_pubtator_cases
from utils.mechanistic_inference import MechanisticInference
from utils.tissue_classifier import FEATURE_NAMES, TissueClassifier, tissue_features

PROJECT_ROOT = Path(__file__).parent
DEFAULT_MODEL_OUT = PROJECT_ROOT / "data" / "tissue_classifier.joblib"
DEFAULT_DATA_OUT = PROJECT_ROOT / "outputs" / "tissue_training_data.csv"


def build_dataset(
    mi: MechanisticInference, cases: List[Dict[str, Any]], loo: bool, source: str
) -> pd.DataFrame:
    """對每個 case 抽組織特徵，同時記下現行規則式判定當基準線。"""
    full_kb = list(mi._load_kb())
    rows = []
    for case in cases:
        if loo:
            # 同位置的變異一起移除：留著等於把答案放在 0bp 的類比裡
            mi._kb = [r for r in full_kb if r["pos"] != case["pos"]]
        try:
            features = mi.build_features(case["variant_info"])
            heuristic = mi._rule_based_fallback(features)
        except Exception:  # noqa: BLE001
            continue

        row = dict(tissue_features(features))
        row.update({
            "source": source,
            "variant_id": case["variant_id"],
            "pos": case["pos"],
            "truth_tissue": case["truth_tissue"],
            "heuristic_tissue": heuristic.get("inferred_tissue"),
        })
        rows.append(row)

    mi._kb = full_kb
    df = pd.DataFrame(rows)
    # 只留二元標籤明確的案例：Both / Uncertain 無法用來訓練二元分類器
    df = df[df["truth_tissue"].isin(["Cardiac", "Skeletal"])].copy()
    df["y"] = (df["truth_tissue"] == "Cardiac").astype(int)
    return df


def make_model(kind: str):
    """
    kind 後綴 "-bal" 表示啟用 class_weight="balanced"。

    平衡與否是這個問題的核心取捨，不是實作細節：HGMD 的組織分佈是 70/30、
    PubTator 是 82/18。開了 balanced 之後模型會大量改判骨骼肌，同源 balanced
    accuracy 漂亮，但換到基準率更偏心肌的 PubTator 就整片誤判。
    """
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    balanced = "balanced" if kind.endswith("-bal") else None
    if kind.startswith("logreg"):
        return make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, class_weight=balanced, C=1.0),
        )
    return GradientBoostingClassifier(
        n_estimators=200, max_depth=3, learning_rate=0.05, random_state=0
    )


def accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float((y_true == y_pred).mean()) if len(y_true) else float("nan")


def balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """兩類各自的正確率再平均。偏斜資料下這個才看得出模型有沒有真的在分類。"""
    parts = [
        (y_pred[y_true == c] == c).mean() for c in (0, 1) if (y_true == c).any()
    ]
    return float(np.mean(parts)) if parts else float("nan")


def heuristic_scores(df: pd.DataFrame) -> Dict[str, float]:
    pred = (df["heuristic_tissue"] == "Cardiac").astype(int).to_numpy()
    y = df["y"].to_numpy()
    return {"acc": accuracy(y, pred), "bal": balanced_accuracy(y, pred)}


def evaluate(train: pd.DataFrame, test: pd.DataFrame, kind: str) -> Dict[str, float]:
    model = make_model(kind)
    model.fit(train[FEATURE_NAMES].to_numpy(), train["y"].to_numpy())
    y = test["y"].to_numpy()
    pred = model.predict(test[FEATURE_NAMES].to_numpy())
    return {"acc": accuracy(y, pred), "bal": balanced_accuracy(y, pred), "n": len(test)}


def cross_val(df: pd.DataFrame, kind: str, folds: int = 5) -> Dict[str, float]:
    """同源交叉驗證。依 genomic 位置分組切 fold，避免同一個群聚同時出現在訓練和測試。"""
    from sklearn.model_selection import GroupKFold

    # 以 10kb 為一個 group：TTN 的致病變異高度群聚，隨機切 fold 會讓同一群聚
    # 的鄰居橫跨訓練/測試集，分數會虛高。
    groups = (df["pos"] // 10000).to_numpy()
    X, y = df[FEATURE_NAMES].to_numpy(), df["y"].to_numpy()
    accs, bals = [], []
    for tr, te in GroupKFold(n_splits=folds).split(X, y, groups):
        model = make_model(kind)
        model.fit(X[tr], y[tr])
        pred = model.predict(X[te])
        accs.append(accuracy(y[te], pred))
        bals.append(balanced_accuracy(y[te], pred))
    return {"acc": float(np.mean(accs)), "bal": float(np.mean(bals)), "n": len(df)}


def combined_cv(df: pd.DataFrame, kind: str, folds: int = 5) -> Dict[str, Dict[str, float]]:
    """合併兩個來源訓練，但預測結果分來源評分。"""
    from sklearn.model_selection import GroupKFold

    groups = (df["pos"] // 10000).to_numpy()
    X, y = df[FEATURE_NAMES].to_numpy(), df["y"].to_numpy()
    pred = np.zeros_like(y)
    for tr, te in GroupKFold(n_splits=folds).split(X, y, groups):
        model = make_model(kind)
        model.fit(X[tr], y[tr])
        pred[te] = model.predict(X[te])

    out = {}
    for src in ("hgmd", "pubtator"):
        m = (df["source"] == src).to_numpy()
        out[src] = {
            "acc": accuracy(y[m], pred[m]),
            "bal": balanced_accuracy(y[m], pred[m]),
            "n": int(m.sum()),
        }
    return out


def _row(label: str, s: Dict[str, float], base: Dict[str, float] = None) -> str:
    acc, bal = s["acc"], s["bal"]
    cells = f"{acc:>10.1%}{bal:>14.1%}"
    if base:
        cells = f"{acc:>7.1%} {acc - base['acc']:+6.1%}{bal:>8.1%} {bal - base['bal']:+6.1%}"
    return f"  {label:<34}{s.get('n', ''):>6}{cells}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="logreg,logreg-bal,gbdt",
                    help="逗號分隔；後綴 -bal 代表 class_weight=balanced")
    ap.add_argument("--save", action="store_true", help="評估後訓練最終模型並存檔")
    ap.add_argument("--model-out", type=Path, default=DEFAULT_MODEL_OUT)
    ap.add_argument("--data-out", type=Path, default=DEFAULT_DATA_OUT)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument(
        "--all-evidence", dest="pathogenic_only", action="store_false", default=True,
        help="連 VUS / 良性證據的 PubTator 標註也納入。預設只採信 P/LP 證據——"
             "不濾的話約三分之一的「錯誤」其實是標註雜訊（收案世代的疾病），"
             "而且會讓結構特徵看起來完全沒有資訊量",
    )
    args = ap.parse_args()

    mi = build_engine("band", hotspot=True)

    print("建立 HGMD 資料集（leave-one-out）…", flush=True)
    hgmd = build_dataset(mi, build_hgmd_cases(mi, args.limit), loo=True, source="hgmd")
    print(f"建立 PubTator 資料集（{'僅 P/LP 證據' if args.pathogenic_only else '全部證據'}）…",
          flush=True)
    pt = build_dataset(
        mi,
        build_pubtator_cases(
            mi, PROJECT_ROOT / "pubtator_only_variants.csv", args.limit,
            pathogenic_only=args.pathogenic_only,
        ),
        loo=False,
        source="pubtator",
    )

    args.data_out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat([hgmd, pt]).to_csv(args.data_out, index=False)

    for name, df in (("HGMD", hgmd), ("PubTator", pt)):
        share = df["y"].mean()
        print(f"{name}: {len(df)} 顆（心肌 {share:.0%} / 骨骼肌 {1 - share:.0%}）")

    kinds = [k.strip() for k in args.model.split(",") if k.strip()]
    h_hgmd, h_pt = heuristic_scores(hgmd), heuristic_scores(pt)

    print(f"\n{'=' * 78}\n組織判定準確率\n{'=' * 78}")
    print(f"  {'':<34}{'n':>6}{'accuracy':>10}{'balanced acc':>14}")
    print("\n── HGMD（同源，10kb 分組交叉驗證）")
    print(_row("現行規則式（基準線）", {**h_hgmd, "n": len(hgmd)}))
    for kind in kinds:
        print(_row(kind, cross_val(hgmd, kind), h_hgmd))

    print("\n── PubTator（跨來源：HGMD 訓練 → PubTator 測試）")
    print(_row("現行規則式（基準線）", {**h_pt, "n": len(pt)}))
    for kind in kinds:
        print(_row(kind, evaluate(hgmd, pt, kind), h_pt))

    print("\n── 反向（PubTator 訓練 → HGMD 測試）")
    print(_row("現行規則式（基準線）", {**h_hgmd, "n": len(hgmd)}))
    for kind in kinds:
        print(_row(kind, evaluate(pt, hgmd, kind), h_hgmd))

    # 兩個來源幾乎不重疊（只有 32 顆共同變異），合併訓練等於同時看到兩種基準率，
    # 應該比只吃單一來源穩健。分 fold 時仍照位置分組，並且分來源回報，才能和
    # 各自的基準線直接比。
    both_df = pd.concat([hgmd, pt], ignore_index=True)
    print("\n── 合併訓練（HGMD+PubTator，10kb 分組交叉驗證，分來源回報）")
    print(_row("現行規則式（HGMD）", {**h_hgmd, "n": len(hgmd)}))
    print(_row("現行規則式（PubTator）", {**h_pt, "n": len(pt)}))
    for kind in kinds:
        per_source = combined_cv(both_df, kind)
        for src, base in (("hgmd", h_hgmd), ("pubtator", h_pt)):
            print(_row(f"{kind}（{src}）", per_source[src], base))

    # 特徵重要性：確認模型學到的是機制上講得通的東西，不是雜訊
    print("\n── 特徵權重（logreg，HGMD+PubTator 全量；正=偏心肌）")
    lr = make_model("logreg")
    lr.fit(both_df[FEATURE_NAMES].to_numpy(), both_df["y"].to_numpy())
    coefs = sorted(
        zip(FEATURE_NAMES, lr.named_steps["logisticregression"].coef_[0]),
        key=lambda kv: -abs(kv[1]),
    )
    for name, c in coefs[:10]:
        print(f"    {name:<32}{c:+.3f}")

    if args.save:
        kind = kinds[0] if len(kinds) == 1 else "logreg"
        model = make_model(kind)
        model.fit(both_df[FEATURE_NAMES].to_numpy(), both_df["y"].to_numpy())
        clf = TissueClassifier(model, FEATURE_NAMES, {
            "kind": kind,
            "n_train": int(len(both_df)),
            "sources": ["hgmd_loo", "pubtator"],
            "cardiac_share": float(both_df["y"].mean()),
            "baseline_hgmd_acc": h_hgmd["acc"],
            "baseline_pubtator_acc": h_pt["acc"],
        })
        clf.save(args.model_out)
        print(f"\n模型已存到 {args.model_out}")

    print(f"訓練資料已存到 {args.data_out}")


if __name__ == "__main__":
    main()
