#!/usr/bin/env python3
"""實驗二／實驗三共用的評估工具（無 GPU、無網路）。

兩個實驗問的是不同問題，但用同一套 ground truth 與同一組指標：

  實驗二  換模型（MedGemma-27B vs Gemma 4 31B）對抽取品質的影響
  實驗三  開關四個防幻覺機制（alias gate / guided JSON / evidence grounding /
          self-consistency）對抽取品質的影響

Ground truth 沿用論文 3.4 節 Strategy 1 的做法：

  HGMD_TTN_with_clinvar_pubcount.xlsx 的 `phenotype` 欄（cardiac / skeletal /
  both / unsure）是真值；預測值是「有提到輸入變異」的文獻的 tissue 票數
  （mention_Cardiac / mention_Skeletal / mention_Both / mention_Not_Specified）
  取 argmax。mention 與 tissue 的判定直接呼叫 utils/html_report.py 裡產生報告
  用的同兩個函式，避免評估與報告漂移。

與論文 3.4 節的差別只有一處：那裡的 122 顆變異是「mention 票數 > 0」，也就是
用當時那個模型的輸出決定評估集，換模型或關掉 gate 後評估集會跟著變，數字不可
比。這裡改成固定評估集（曾經處理過 ≥1 篇文獻的 582 顆變異，或其分層子集），
把「一票都沒有」記成 abstain 而不是排除，因此每個設定都報三個數字：

  coverage   有給出答案的比例（1 - abstain rate）
  accuracy   在自己有答案的子集上的正確率（可與論文 78.7% 對照）
  overall    在固定評估集上的正確率（abstain 一律算錯）
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "research" / "experiments"
CORPUS_DIR = EXP_DIR / "corpus"          # build_llm_corpus.py 的輸出（兩輪共用，以變異為鍵）

# 兩輪評估集並存，用環境變數 TTN_EXP_SET 切換，免得第二輪覆蓋第一輪的結果：
#   未設定  第一輪：從 582 顆池子分層抽樣的 60 顆（結果已封存在 round1_random60/）
#   "hub"   第二輪：提及最多變異的 5 篇 hub 文獻所提及的 57 顆，且只餵那幾篇文獻
#           （評估集由 find_hub_articles.py 產生）
SET_NAME = os.environ.get("TTN_EXP_SET", "").strip()
SUFFIX = f"_{SET_NAME}" if SET_NAME else ""
ROUND1_DIR = EXP_DIR / "round1_random60"
RUNS_DIR = EXP_DIR / f"runs_{SET_NAME}" if SET_NAME else ROUND1_DIR / "runs"
EVAL_SET_CSV = (EXP_DIR / f"llm_eval_set{SUFFIX}.csv" if SET_NAME
                else ROUND1_DIR / "llm_eval_set.csv")

GT_XLSX = PROJECT_ROOT / "HGMD_TTN_with_clinvar_pubcount.xlsx"

TISSUES = ["Cardiac", "Skeletal", "Both", "Not_Specified"]
MENTION_COLS = [f"mention_{t}" for t in TISSUES]
NO_MENTION_COLS = [f"no_mention_{t}" for t in TISSUES]
# mention_* 四欄的 argmax → ground truth phenotype 的標籤（順序必須與 TISSUES 對齊）
LABELS4 = ["cardiac", "skeletal", "both", "others"]
# 真值裡 unsure 與預測的 others 是同一類（都代表「文獻沒指明組織」）
GT_LABEL_MAP = {"cardiac": "cardiac", "skeletal": "skeletal", "both": "both",
                "unsure": "others", "others": "others"}


# ── Ground truth ─────────────────────────────────────────────────────────────

def load_gt() -> pd.DataFrame:
    """讀 HGMD 真值表，附上「當初有幾篇文獻被處理過」。"""
    df = pd.read_excel(GT_XLSX)
    df["articles_seen"] = df[MENTION_COLS + NO_MENTION_COLS].sum(axis=1)
    df["gt_label"] = df["phenotype"].map(GT_LABEL_MAP)
    return df


def build_eval_set(
    n: int = 60,
    seed: int = 0,
    min_per_class: int = 8,
    out_csv: Optional[Path] = EVAL_SET_CSV,
) -> pd.DataFrame:
    """從「曾經處理過 ≥1 篇文獻」的變異中抽分層子集，當成兩個實驗的固定評估集。

    both / unsure 兩類本來就少（14 / 4 顆），全數納入；cardiac / skeletal 依
    比例抽，但 skeletal 至少 `min_per_class` 顆，否則 confusion matrix 會沒東西
    可看。n 給 None 或 0 代表整池（582 顆）都要。
    """
    gt = load_gt()
    pool = gt[(gt["articles_seen"] > 0) & gt["input"].notna()].copy()

    if not n:
        picked = pool
    else:
        rng = np.random.default_rng(seed)
        chunks: List[pd.DataFrame] = []
        rare = pool[pool["gt_label"].isin(["both", "others"])]
        chunks.append(rare)
        budget = max(0, n - len(rare))
        common = pool[pool["gt_label"].isin(["cardiac", "skeletal"])]
        share = common["gt_label"].value_counts(normalize=True)
        quota = {
            "skeletal": max(min_per_class, int(round(budget * share.get("skeletal", 0)))),
        }
        quota["cardiac"] = max(0, budget - quota["skeletal"])
        for label, k in quota.items():
            sub = common[common["gt_label"] == label]
            k = min(k, len(sub))
            if k:
                idx = rng.choice(sub.index.to_numpy(), size=k, replace=False)
                chunks.append(sub.loc[idx])
        picked = pd.concat(chunks).sort_values("articles_seen", ascending=False)

    cols = ["input", "#CHROM", "POS", "REF", "ALT", "ID", "DNA", "PROT", "PHEN",
            "phenotype", "gt_label", "articles_seen", "ClinVar_Variation_ID"]
    picked = picked[[c for c in cols if c in picked.columns]].reset_index(drop=True)
    if out_csv is not None:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        picked.to_csv(out_csv, index=False)
    return picked


def load_eval_set() -> pd.DataFrame:
    """載入目前這一輪的評估集（由 TTN_EXP_SET 決定）。

    第一輪的 CSV 若還沒建過就即時抽樣建出來；第二輪（hub）的 CSV 由
    find_hub_articles.py 產生，找不到就直接報錯，免得默默換成別的集合。
    """
    if EVAL_SET_CSV.exists():
        return pd.read_csv(EVAL_SET_CSV)
    if SET_NAME:
        raise FileNotFoundError(
            f"找不到 {EVAL_SET_CSV}；先跑 "
            f"python research/experiments/find_hub_articles.py --top 5 --out {EVAL_SET_CSV}"
        )
    return build_eval_set()


# ── 把一次抽取的結果換算成 mention 票數 ───────────────────────────────────────

def tissue_counts(articles: Iterable[Dict]) -> Dict[str, int]:
    """把 batch_extract 後的 articles 換算成 8 個計數欄。

    判定邏輯直接借用報告產生器的兩個 staticmethod，確保與 outputs/*.html 裡
    Tissue Involvement Statistics 的數字一致。
    """
    from utils.html_report import HTMLReportGenerator as R

    counts = {c: 0 for c in MENTION_COLS + NO_MENTION_COLS}
    for article in articles:
        mention = R._article_mentions_variant(article)
        tissue = R._normalize_tissue((article.get("clinical_info") or {}).get("tissue_affected"))
        key = tissue.replace(" ", "_")
        prefix = "mention" if mention else "no_mention"
        col = f"{prefix}_{key}"
        if col in counts:
            counts[col] += 1
    return counts


def strategy1_predict(counts: Dict[str, int]) -> Optional[str]:
    """論文 Strategy 1：mention_* 四欄 argmax。四欄全 0 → None（abstain）。"""
    values = [int(counts.get(c, 0) or 0) for c in MENTION_COLS]
    if sum(values) == 0:
        return None
    return LABELS4[int(np.argmax(values))]


# ── 讀 run_extraction_config.py 的輸出 ────────────────────────────────────────

def load_run(path: Path) -> pd.DataFrame:
    """讀一個設定的 JSONL 輸出，回傳「一列一顆變異」的表。"""
    rows: List[Dict] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            counts = rec.get("counts", {})
            articles = rec.get("articles", [])
            row = {
                "input": rec["variant"],
                "config": rec.get("config", path.stem),
                "n_articles": len(articles),
                "pred": strategy1_predict(counts),
                "seconds": rec.get("seconds"),
                # 防幻覺機制的動作次數
                "n_gate_blocked": sum(
                    1 for a in articles if a.get("extraction_source") == "deterministic_no_alias"
                ),
                "n_ungrounded_downgraded": sum(
                    1 for a in articles if "ungrounded" in str(a.get("extraction_source", ""))
                ),
                "n_json_failures": rec.get("json_parse_failures", 0),
                "n_samples_parsed": rec.get("samples_parsed", 0),
                # 事後探針（與機制開關無關，永遠算）
                "n_alias_absent": sum(1 for a in articles if a.get("alias_present") is False),
                "n_claims": sum(
                    1 for a in articles if a.get("disease") not in (None, "", "Not specified",
                                                                    "JSON Parsing Failed")
                ),
                "n_claims_without_alias": sum(
                    1 for a in articles
                    if a.get("alias_present") is False
                    and a.get("disease") not in (None, "", "Not specified", "JSON Parsing Failed")
                ),
                "n_claims_ungrounded_probe": sum(
                    1 for a in articles
                    if a.get("evidence_grounded_probe") is False
                    and a.get("disease") not in (None, "", "Not specified", "JSON Parsing Failed")
                ),
                "n_llm_calls": sum(1 for a in articles if a.get("prompted")),
            }
            row.update(counts)
            rows.append(row)
    return pd.DataFrame(rows)


def load_articles(path: Path) -> pd.DataFrame:
    """同一份 JSONL，但攤平成「一列一篇 (變異, 文獻)」，用來看逐篇差異。"""
    rows: List[Dict] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            for a in rec.get("articles", []):
                rows.append({"input": rec["variant"], "config": rec.get("config", path.stem), **a})
    return pd.DataFrame(rows)


# ── 指標 ─────────────────────────────────────────────────────────────────────

def score_run(run: pd.DataFrame, eval_set: pd.DataFrame) -> Dict[str, float]:
    """一個設定的摘要指標。abstain（pred is None）在 overall 裡算錯。"""
    merged = eval_set[["input", "gt_label"]].merge(run, on="input", how="inner")
    n = len(merged)
    answered = merged[merged["pred"].notna()]
    correct_answered = (answered["pred"] == answered["gt_label"]).sum()

    per_class_f1 = []
    for label in LABELS4:
        tp = int(((merged["pred"] == label) & (merged["gt_label"] == label)).sum())
        fp = int(((merged["pred"] == label) & (merged["gt_label"] != label)).sum())
        fn = int(((merged["pred"] != label) & (merged["gt_label"] == label)).sum())
        if tp + fp + fn == 0:
            continue
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        per_class_f1.append(2 * prec * rec / (prec + rec) if prec + rec else 0.0)

    totals = {c: int(merged[c].sum()) for c in ["n_articles", "n_llm_calls", "n_gate_blocked",
                                                "n_claims", "n_claims_without_alias",
                                                "n_claims_ungrounded_probe",
                                                "n_ungrounded_downgraded", "n_json_failures",
                                                "n_samples_parsed"] if c in merged}
    claims = max(1, totals.get("n_claims", 0))
    samples = max(1, totals.get("n_samples_parsed", 0))
    return {
        "n_variants": n,
        "coverage": len(answered) / n if n else float("nan"),
        "accuracy_answered": correct_answered / len(answered) if len(answered) else float("nan"),
        "accuracy_overall": correct_answered / n if n else float("nan"),
        "macro_f1": float(np.mean(per_class_f1)) if per_class_f1 else float("nan"),
        # 幻覺相關：分母是「這個設定做出的疾病主張數」
        "claims": totals.get("n_claims", 0),
        "claim_rate_no_alias": totals.get("n_claims_without_alias", 0) / claims,
        "claim_rate_ungrounded": totals.get("n_claims_ungrounded_probe", 0) / claims,
        "json_failure_rate": totals.get("n_json_failures", 0) / samples,
        "llm_calls": totals.get("n_llm_calls", 0),
        "minutes": float(merged["seconds"].sum() / 60) if "seconds" in merged else float("nan"),
        **totals,
    }


def confusion(run: pd.DataFrame, eval_set: pd.DataFrame,
              labels: Optional[List[str]] = None) -> pd.DataFrame:
    """真值 × 預測的混淆矩陣；abstain 自成一欄（pred = abstain）。"""
    labels = labels or LABELS4
    merged = eval_set[["input", "gt_label"]].merge(run, on="input", how="inner")
    pred = merged["pred"].fillna("abstain")
    mat = pd.crosstab(merged["gt_label"], pred)
    return mat.reindex(index=labels, columns=labels + ["abstain"], fill_value=0)


def summary_table(runs: Dict[str, pd.DataFrame], eval_set: pd.DataFrame) -> pd.DataFrame:
    """多個設定的指標併成一張表（列 = 設定）。"""
    return pd.DataFrame({name: score_run(run, eval_set) for name, run in runs.items()}).T


def pairwise_agreement(a: pd.DataFrame, b: pd.DataFrame) -> Dict[str, float]:
    """兩個設定之間的一致性（變異層級與文獻層級都看）。"""
    m = a[["input", "pred"]].merge(b[["input", "pred"]], on="input", suffixes=("_a", "_b"))
    same = (m["pred_a"].fillna("abstain") == m["pred_b"].fillna("abstain")).mean()
    return {"n": len(m), "variant_level_agreement": float(same)}


def discover_runs(pattern: str = "*.jsonl", runs_dir: Path = RUNS_DIR) -> Dict[str, Path]:
    return {p.stem: p for p in sorted(runs_dir.glob(pattern))}
