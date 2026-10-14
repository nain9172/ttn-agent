#!/usr/bin/env python3
"""
比較 Track B 機制推理的檢索策略：band（現行）vs exon（外顯子層級分層檢索）。

為什麼需要這支腳本：band 只有 4 格、跨度 28 萬 bp，但文獻談 TTN 的表型專一性
一律以外顯子為單位（RELATED_WORK.md：Savarese 2020、Roberts 2015）。改成外顯子
層級檢索之後，必須有辦法證明準確率真的提升，而不是換個寫法而已。

兩套獨立的 ground truth：

  hgmd      對 HGMD 每顆變異做 leave-one-out：先把它自己（預設連同位置的其他
            變異一起）從類比庫移除，再預測，然後和它自己的 PHEN 比對。
            這是唯一涵蓋整條基因的標註，但和類比庫同源，屬於同源評估。

  pubtator  用 pubtator_only_variants.csv 的文獻標註當答案。資料量小、偏向被
            研究得多的變異，但和 HGMD 類比庫完全獨立，是跨來源的類推驗證。

            預設只採信原文判為 Pathogenic / Likely pathogenic 的標註。這不是可有
            可無的調校：PubTator 的 disease 欄是 LLM 從全文抽的，大量標註來自基因
            套組的疾病清單或世代篩檢表格，抽到的是「收案世代的疾病」而非變異表型
            （最明顯的是 HCM——TTN 並非公認的 HCM 致病基因，卻佔了未過濾集的 17%，
            其中 57% 的原文分級是 VUS 或良性）。不過濾的後果不只是分數低估 10 個
            百分點，還會讓結論整個反過來：組織基準率校正在未過濾集 +1.6%、乾淨集
            -0.4%（它其實在配適雜訊），而結構性特徵在未過濾集看起來毫無資訊量
            （balanced accuracy 50%，等同亂猜）、在乾淨集則有 62%。

組態逐項累加，才知道每個槓桿各貢獻多少：

  baseline    這次改動前的行為（band 檢索、無熱點規則、類比等權投票）
  +hotspot    加上 exon 344/364 熱點硬規則
  +exon       改用外顯子分層檢索
  +rank       類比票數依 HGMD RANKSCORE 打折（不可靠的標註不該等值投票）
  production  再加上籠統疾病類比的折扣，即目前的生產預設

只跑 rule-based（不載 LLM）：這次改動的是檢索與投票邏輯，rule-based 是可重現、
可歸因的那一層；LLM 會把 prompt 變動的效果混進來。

用法：
    conda run -n gm4 python benchmark_inference.py
    conda run -n gm4 python benchmark_inference.py --truth hgmd --limit 300
    conda run -n gm4 python benchmark_inference.py --loo-scope exact
    conda run -n gm4 python benchmark_inference.py --out outputs/bench.csv
"""

from __future__ import annotations

import argparse
import re
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

import sys as _sys
from pathlib import Path as _Path

# 本檔不在 repo 根目錄：把根目錄（config / utils）與本檔所在資料夾放進 import 路徑
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
_sys.path.insert(1, str(_Path(__file__).resolve().parent))

import config
from utils.disease_normalizer import (
    canonical_labels,
    compare_disease,
    is_ttn_relevant,
    is_unspecified,
    normalize_tissue,
    tissue_match,
)
from utils.mechanistic_inference import MechanisticInference

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PUBTATOR_CSV = PROJECT_ROOT / "pubtator_only_variants.csv"
DEFAULT_OUT = PROJECT_ROOT / "outputs" / "benchmark_inference.csv"

# (組態名稱, retrieval_mode, 是否啟用外顯子熱點硬規則)
# (組態名稱, retrieval_mode, 外顯子熱點硬規則, RANKSCORE 折扣強度, 籠統類比折扣)
# 末欄為 inclusion_lr_strength：外顯子 isoform 收錄度（PSI 代理）→ 組織 LR 的強度，
# 0 = 停用。詳見 mechanistic_inference.inclusion_tissue_lr()。
CONFIGS = [
    ("baseline", "band", False, 0.0, 1.0, 0.0),  # 這次改動前的行為
    ("+hotspot", "band", True, 0.0, 1.0, 0.0),
    ("+exon", "exon", True, 0.0, 1.0, 0.0),
    ("+rank", "exon", True, 0.5, 1.0, 0.0),
    ("+vague", "exon", True, 0.5, 0.6, 0.0),  # 前一版的生產預設
    ("+psi", "exon", True, 0.5, 0.6, 2.0),  # 加上收錄度組織 LR（生產預設）
    ("band-hit", "band", True, 0.0, 1.0, 2.0),  # hit 最高，但 specific/tissue 較低
]
PRODUCTION_CONFIG = ("exon", True, 0.5, 0.6, 2.0)

HIT_LEVELS = ("exact", "compatible")


def build_engine(
    retrieval_mode: str,
    hotspot: bool,
    tissue_prior: Optional[Dict[str, float]] = None,
    rankscore: float = 0.5,
    vague_penalty: float = 0.6,
    inclusion_lr: float = 1.0,
) -> MechanisticInference:
    return MechanisticInference(
        kb_path=config.HGMD_KB_PATH,
        band_edges=config.TTN_BAND_EDGES,
        band_labels=config.TTN_BAND_LABELS,
        top_k=config.INFERENCE_TOP_K_ANALOGS,
        transcript_interval_path=config.TRANSCRIPT_INTERVAL_PATH,
        transcript_model_path=config.TTN_TRANSCRIPT_MODEL_PATH,
        reference_fasta_path=config.REFERENCE_GENOME_PATH,
        retrieval_mode=retrieval_mode,
        exon_hotspot_rules=hotspot,
        tissue_prior=tissue_prior,
        rankscore_weighting=rankscore,
        vague_penalty=vague_penalty,
        inclusion_tissue_weighting=inclusion_lr > 0,
        inclusion_lr_strength=inclusion_lr,
    )


# ── 評估集 ───────────────────────────────────────────────────────────────────


def build_hgmd_cases(mi: MechanisticInference, limit: int = 0) -> List[Dict[str, Any]]:
    """HGMD 每顆變異一個 case，ground truth 是它自己的 PHEN / phenotype。"""
    cases = []
    for r in mi._load_kb():
        disease = (r.get("phen") or "").strip()
        if not disease or is_unspecified(disease):
            continue
        labels = canonical_labels(disease)
        if not is_ttn_relevant(labels):
            continue
        cases.append({
            "variant_id": f"2-{r['pos']}-{r['ref']}-{r['alt']}",
            "pos": r["pos"],
            "variant_info": {
                "pos": r["pos"], "ref": r["ref"], "alt": r["alt"],
                "hgvs_c": r.get("dna"), "hgvs_p": r.get("prot"),
            },
            "truth_disease": disease,
            "truth_tissue": normalize_tissue(r.get("phenotype") or disease),
            "tissue_is_explicit": bool(r.get("phenotype")),
        })
    cases.sort(key=lambda c: c["pos"])
    return _subsample(cases, limit)


# PubTator 的 tissue_affected 只有這幾個值算「文獻真的講了組織」；
# "Not specified" 和 "Others" 都不算。
EXPLICIT_TISSUES = {"cardiac", "skeletal", "both"}

# 只有這兩級的證據能用來定義「這顆變異造成什麼病」。
#
# 為什麼一定要濾：PubTator 的 disease 欄是 LLM 從全文抽的，很多來自基因套組的
# 疾病清單或世代篩檢表格，抽到的是「收案世代的疾病」而不是變異本身的表型。
# 最明顯的是 HCM——TTN 並非公認的 HCM 致病基因，但 376 列 HCM 標註裡有 148 列
# 的變異被原文分類為 VUS、42 列可能良性、27 列良性。拿這種標註當 ground truth，
# 量到的是標註雜訊，不是模型能力。
PATHOGENIC_EVIDENCE = {"pathogenic", "likely pathogenic"}

# pathogenicity 欄還不夠：原文寫 "VUS-LP"（VUS 偏可能致病，本質仍是 VUS）時，
# 抽取會把 "VUS-" 前綴丟掉、記成 "Likely pathogenic"。這類列佔 P/LP 的 8%，
# 而且高度集中在 HCM（89 列裡有 55 列），得回頭看原始證據句才擋得掉。
_VUS_IN_EVIDENCE = re.compile(r"\bVUS\b", re.IGNORECASE)


def build_pubtator_cases(
    mi: MechanisticInference,
    csv_path: Path,
    limit: int = 0,
    explicit_tissue_only: bool = False,
    pathogenic_only: bool = False,
) -> List[Dict[str, Any]]:
    """
    重用 test_inference_pubtator_only 的 HGVS 解析 / 多列合併邏輯。

    組織 ground truth 分兩種來源，必須分開看：文獻明確標了 tissue_affected 的，
    以及沒標、只能從疾病名反推的。後者會讓「組織準確率」和「疾病準確率」變成同
    一件事的兩種寫法（反推本身就是用疾病名做的）。explicit_tissue_only=True 只留前者。

    pathogenic_only=True 時，只用原文把該變異判為 Pathogenic / Likely pathogenic
    的那些列來決定疾病與組織（見 PATHOGENIC_EVIDENCE 的說明）。
    """
    from test_inference_pubtator_only import _best_disease, _mode, collect_variants

    groups, _ = collect_variants(csv_path, mi)
    cases = []
    for vid, g in groups.items():
        rows = g["rows"]
        if pathogenic_only:
            rows = [
                r for r in rows
                if str(r.get("pathogenicity") or "").strip().lower() in PATHOGENIC_EVIDENCE
                and not _VUS_IN_EVIDENCE.search(str(r.get("evidence") or ""))
            ]
            if not rows:
                continue

        disease = _best_disease([r["disease"] for r in rows])
        if not disease:
            continue
        labels = canonical_labels(disease)
        if not is_ttn_relevant(labels):
            continue

        explicit = [
            r["tissue"] for r in rows
            if r["tissue"] and str(r["tissue"]).strip().lower() in EXPLICIT_TISSUES
        ]
        if explicit_tissue_only and not explicit:
            continue
        # 有明確標註就用明確標註的多數，沒有才退回從疾病名反推
        tissue_is_explicit = bool(explicit)
        tissue_src = _mode(explicit) if explicit else disease

        rv = g["resolved"]
        cases.append({
            "variant_id": vid,
            "pos": rv.pos,
            "variant_info": rv.to_variant_info(),
            "truth_disease": disease,
            "truth_tissue": normalize_tissue(tissue_src),
            "tissue_is_explicit": tissue_is_explicit,
        })
    cases.sort(key=lambda c: c["pos"])
    return _subsample(cases, limit)


def _subsample(cases: List[Dict], limit: int) -> List[Dict]:
    """有 --limit 時等距抽樣，避免只取到基因某一端的 band。"""
    if not limit or limit >= len(cases):
        return cases
    stride = len(cases) / limit
    return [cases[int(i * stride)] for i in range(limit)]


# ── 執行單一組態 ─────────────────────────────────────────────────────────────


def run_config(
    mi: MechanisticInference,
    cases: List[Dict[str, Any]],
    loo: bool,
    loo_scope: str,
    config_name: str,
    truth_name: str,
) -> List[Dict[str, Any]]:
    """跑完一個組態的所有 case，回傳逐筆評分結果。"""
    full_kb = list(mi._load_kb())
    rows = []
    t0 = time.time()

    for i, case in enumerate(cases):
        if loo:
            # leave-one-out：把答案（以及預設連同一位置的其他變異）抽掉再預測，
            # 否則 hgmd_exact_match 會直接把答案讀出來，量到的是查表不是推論。
            pos, ref, alt = case["pos"], case["variant_info"]["ref"], case["variant_info"]["alt"]
            if loo_scope == "position":
                mi._kb = [r for r in full_kb if r["pos"] != pos]
            else:
                mi._kb = [
                    r for r in full_kb
                    if not (r["pos"] == pos and r["ref"] == ref and r["alt"] == alt)
                ]

        try:
            features = mi.build_features(case["variant_info"])
            inf = mi._rule_based_fallback(features)
        except Exception as e:  # noqa: BLE001  單顆失敗不該中斷整輪
            print(f"  ! {case['variant_id']} 推理失敗：{type(e).__name__}: {e}", flush=True)
            continue

        pred_disease = inf.get("inferred_disease")
        differentials = inf.get("differential_diagnoses") or []
        cmp_res = compare_disease(pred_disease, case["truth_disease"], differentials)
        pred_tissue = normalize_tissue(inf.get("inferred_tissue"))
        analogs = features.get("analogs", [])

        rows.append({
            "truth_set": truth_name,
            "config": config_name,
            "variant_id": case["variant_id"],
            "pos": case["pos"],
            "band": features.get("band"),
            "exon": features.get("exon"),
            "vclass": features.get("vclass"),
            "vgroup": features.get("vgroup"),
            "n_analogs": len(analogs),
            "nearest_analog_bp": min([a.get("distance", 10**9) for a in analogs], default=None),
            "top_tier": analogs[0].get("match_tier") if analogs else None,
            "n_exon_analogs": sum(
                1 for a in analogs if str(a.get("match_tier", "")).startswith("exon")
            ),
            "truth_disease": case["truth_disease"],
            "truth_labels": ",".join(cmp_res["truth_labels"]),
            "truth_tissue": case["truth_tissue"],
            "tissue_is_explicit": case.get("tissue_is_explicit"),
            "pred_disease": pred_disease,
            "pred_labels": ",".join(cmp_res["pred_labels"]),
            "pred_tissue": pred_tissue,
            "pred_confidence": inf.get("confidence"),
            "pred_source": inf.get("source", "analog_vote"),
            "disease_eval": cmp_res["level"],
            "tissue_eval": tissue_match(pred_tissue, case["truth_tissue"]),
            "pred_vague": cmp_res["pred_vague"],
        })

        if (i + 1) % 250 == 0:
            print(f"    {i + 1}/{len(cases)}  ({time.time() - t0:.0f}s)", flush=True)

    mi._kb = full_kb
    return rows


# ── 指標 ─────────────────────────────────────────────────────────────────────


def metrics(df: pd.DataFrame) -> Dict[str, float]:
    """
    hit       疾病命中率：預測等於答案，或和答案相容（一般詞 ↔ 特定病名）。
    specific  hit 中再排除只答 Titinopathy / Cardiomyopathy 這類籠統詞的。
    macro     先算每個答案類別各自的 hit 率再平均。TTN 的答案分佈極度偏斜
              （DCM 佔大宗），hit 會被多數類別灌水，macro 才看得出少數表型。
    tissue    組織相容率（exact + partial），只計有標到組織的 case。
    """
    n = len(df)
    if not n:
        return {}
    d = df["disease_eval"]
    hit_mask = d.isin(HIT_LEVELS)
    vague = df["pred_vague"].fillna(False).astype(bool)

    per_class = [
        g["disease_eval"].isin(HIT_LEVELS).mean()
        for _, g in df.groupby("truth_labels") if len(g) >= 3
    ]
    t = df["tissue_eval"]
    t_scored = (t != "n/a").sum()

    return {
        "n": n,
        "hit": hit_mask.mean(),
        "specific": (hit_mask & ~vague).mean(),
        "macro": sum(per_class) / len(per_class) if per_class else float("nan"),
        "tissue": (t.isin(["exact", "partial"]).sum() / t_scored) if t_scored else float("nan"),
    }


def print_comparison(df: pd.DataFrame, truth_name: str, baseline: str = "band") -> None:
    print(f"\n{'=' * 78}\n{truth_name} ground truth\n{'=' * 78}")
    order = list(dict.fromkeys(df["config"]))  # 保留執行順序
    base = metrics(df[df["config"] == baseline]) if baseline in order else {}

    keys = ("hit", "specific", "macro", "tissue")
    print(f"{'組態':<14}{'n':>6}" + "".join(f"{k:>16}" for k in keys))
    for cfg in order:
        m = metrics(df[df["config"] == cfg])
        if not m:
            continue
        cells = ""
        for key in keys:
            val = m[key]
            if val != val:  # NaN
                cells += f"{'—':>16}"
            elif cfg == baseline or not base or base.get(key) != base.get(key):
                cells += f"{val:>16.1%}"
            else:
                cells += f"{val:>9.1%} {val - base[key]:+6.1%}"
        print(f"{cfg:<14}{m['n']:>6}{cells}")

    # 只看兩個組態意見不同的 case，才知道改動實際動到了什麼
    if baseline in order:
        for cfg in order:
            if cfg == baseline:
                continue
            _print_disagreement(df, baseline, cfg)


def _print_disagreement(df: pd.DataFrame, base_cfg: str, new_cfg: str) -> None:
    a = df[df["config"] == base_cfg].set_index("variant_id")
    b = df[df["config"] == new_cfg].set_index("variant_id")
    shared = a.index.intersection(b.index)
    if not len(shared):
        return
    a_hit = a.loc[shared, "disease_eval"].isin(HIT_LEVELS)
    b_hit = b.loc[shared, "disease_eval"].isin(HIT_LEVELS)
    fixed = int((~a_hit & b_hit).sum())
    broke = int((a_hit & ~b_hit).sum())
    changed = int((a.loc[shared, "pred_labels"].fillna("") != b.loc[shared, "pred_labels"].fillna("")).sum())
    print(f"  {base_cfg} → {new_cfg}：預測改變 {changed} 顆，其中修好 {fixed}、弄壞 {broke}"
          f"（淨 {fixed - broke:+d}）")


def print_breakdown(df: pd.DataFrame, cfg: str) -> None:
    """新組態的錯誤分佈：看提升／退步集中在哪裡。"""
    sub = df[df["config"] == cfg]
    if sub.empty:
        return
    print(f"\n[{cfg}] 依最近類比的定位層級")
    for tier, g in sub.groupby("top_tier", dropna=False):
        hit = g["disease_eval"].isin(HIT_LEVELS).mean()
        print(f"  {str(tier):<12} {len(g):>5} 顆   hit {hit:.1%}")

    # 明確標註 vs 從疾病名反推：如果兩邊分數差很多，代表現在的總分有一部分
    # 是在量標註品質，而不是在量模型
    if sub["tissue_is_explicit"].notna().any():
        print(f"\n[{cfg}] 依組織標註來源")
        for explicit, g in sub.groupby(sub["tissue_is_explicit"].astype(bool)):
            m = metrics(g)
            label = "文獻明確標註" if explicit else "由疾病名反推"
            print(f"  {label:<14} {m['n']:>5} 顆   hit {m['hit']:.1%}   "
                  f"specific {m['specific']:.1%}   tissue {m['tissue']:.1%}")

    print(f"\n[{cfg}] 依 band")
    for band, g in sub.groupby("band", dropna=False):
        hit = g["disease_eval"].isin(HIT_LEVELS).mean()
        print(f"  {str(band):<12} {len(g):>5} 顆   hit {hit:.1%}")

    print(f"\n[{cfg}] 最常見的錯誤（前 8 組）")
    miss = sub[~sub["disease_eval"].isin(HIT_LEVELS)]
    pairs = Counter(
        (r.truth_labels or "?", r.pred_labels or "?") for r in miss.itertuples()
    )
    for (t, p), n in pairs.most_common(8):
        print(f"  {n:>4}  truth={str(t)[:28]:<30} pred={str(p)[:28]:<30}")


# ── 主流程 ───────────────────────────────────────────────────────────────────


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", choices=["hgmd", "pubtator", "both"], default="both")
    ap.add_argument("--configs", default=",".join(c[0] for c in CONFIGS),
                    help="要比較的組態，逗號分隔")
    ap.add_argument("--limit", type=int, default=0, help="每套 ground truth 最多測幾顆（0=全部）")
    ap.add_argument("--loo-scope", choices=["position", "exact"], default="position",
                    help="HGMD leave-one-out 移除範圍：position=同位置全移（較嚴格）")
    ap.add_argument("--pubtator-csv", type=Path, default=DEFAULT_PUBTATOR_CSV)
    ap.add_argument(
        "--pubtator-tissue", choices=["all", "explicit"], default="all",
        help="explicit=只留文獻真的標了 tissue_affected 的變異；all 會把沒標的用疾病名反推",
    )
    # 預設就過濾：不濾的話 1028 顆裡有 465 顆的「答案」來自 VUS / 良性變異的
    # 收案世代疾病，量出來的命中率會低估約 10 個百分點，而且會讓「把預測往心臟推」
    # 這種迎合雜訊的改動看起來有效（基準率校正在髒資料 +1.6%、乾淨資料 -0.4%）。
    ap.add_argument(
        "--pubtator-all-evidence", dest="pubtator_pathogenic",
        action="store_false", default=True,
        help="連 VUS / 良性證據的標註也納入（預設只採信 P/LP 證據）",
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument(
        "--prior-sweep", default="",
        help="逗號分隔的心臟基準率，例 0.68,0.75,0.85。會在 band 組態上加測這些校正強度；"
             "同一組設定要在 HGMD 和 PubTator 兩邊都變好，才不是配適雜訊",
    )
    ap.add_argument(
        "--vague-sweep", default="",
        help="逗號分隔的籠統疾病類比折扣，例 0.5,0.25。在生產預設組態上加測；"
             "HGMD 約三成 PHEN 是 Cardiomyopathy / Titinopathy 這類上位詞",
    )
    ap.add_argument("--incl-sweep", default="",
                    help="掃描 inclusion_lr_strength，例如 0,0.5,1,1.5,2")
    args = ap.parse_args()

    wanted = {c.strip() for c in args.configs.split(",") if c.strip()}
    configs = [c for c in CONFIGS if c[0] in wanted]

    engines = {name: build_engine(mode, hot, rankscore=rank, vague_penalty=vag,
                                  inclusion_lr=incl)
               for name, mode, hot, rank, vag, incl in configs}
    p_mode, p_hot, p_rank, p_vag, p_incl = PRODUCTION_CONFIG
    for raw in (p.strip() for p in args.prior_sweep.split(",") if p.strip()):
        card = float(raw)
        name = f"prior{card:.2f}"
        # 只校正 cardiac / skeletal 兩類，其餘（both / unsure）維持原權重
        engines[name] = build_engine(
            p_mode, p_hot, tissue_prior={"cardiac": card, "skeletal": 1.0 - card},
            rankscore=p_rank, vague_penalty=p_vag, inclusion_lr=p_incl,
        )
        configs.append((name, p_mode, p_hot, p_rank, p_vag, p_incl))
    for raw in (v.strip() for v in args.vague_sweep.split(",") if v.strip()):
        pen = float(raw)
        name = f"vague{pen:.2f}"
        engines[name] = build_engine(p_mode, p_hot, rankscore=p_rank, vague_penalty=pen,
                                     inclusion_lr=p_incl)
        configs.append((name, p_mode, p_hot, p_rank, pen, p_incl))
    for raw in (v.strip() for v in args.incl_sweep.split(",") if v.strip()):
        st = float(raw)
        name = f"incl{st:.2f}"
        engines[name] = build_engine(p_mode, p_hot, rankscore=p_rank,
                                     vague_penalty=p_vag, inclusion_lr=st)
        configs.append((name, p_mode, p_hot, p_rank, p_vag, st))
    if not configs:
        raise SystemExit(f"沒有符合的組態：{args.configs}")
    probe = next(iter(engines.values()))

    truth_sets: Dict[str, List[Dict]] = {}
    if args.truth in ("hgmd", "both"):
        truth_sets["HGMD (leave-one-out)"] = build_hgmd_cases(probe, args.limit)
    if args.truth in ("pubtator", "both"):
        if args.pubtator_csv.exists():
            explicit_only = args.pubtator_tissue == "explicit"
            quals = [q for q, on in (("明確組織標註", explicit_only),
                                     ("僅 P/LP 證據", args.pubtator_pathogenic)) if on]
            label = "PubTator (跨來源" + ("，" + "，".join(quals) if quals else "") + ")"
            truth_sets[label] = build_pubtator_cases(
                probe, args.pubtator_csv, args.limit,
                explicit_tissue_only=explicit_only,
                pathogenic_only=args.pubtator_pathogenic,
            )
        else:
            print(f"找不到 {args.pubtator_csv}，跳過 PubTator")

    all_rows: List[Dict] = []
    for truth_name, cases in truth_sets.items():
        loo = truth_name.startswith("HGMD")
        print(f"\n### {truth_name}：{len(cases)} 顆變異"
              f"{f'（LOO 移除範圍={args.loo_scope}）' if loo else ''}", flush=True)
        for name in engines:
            print(f"  跑組態 {name} …", flush=True)
            all_rows += run_config(
                engines[name], cases, loo, args.loo_scope, name, truth_name
            )

    if not all_rows:
        print("沒有可評分的變異。")
        return

    df = pd.DataFrame(all_rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    for truth_name in truth_sets:
        sub = df[df["truth_set"] == truth_name]
        if sub.empty:
            continue
        print_comparison(sub, truth_name)
        best = configs[-1][0]
        if best in set(sub["config"]):
            print_breakdown(sub, best)

    print(f"\n明細已寫入：{args.out}（{len(df)} 列）")


if __name__ == "__main__":
    main()
