#!/usr/bin/env python3
"""
用 pubtator_only_variants.csv 批量測試 Track B 機制推理，並和檔案裡的 disease 比對。

pubtator_only_variants.csv 是「只有 PubTator 文獻標註」的 TTN 變異（status=pubtator_only），
每列有 hgvs / rsid 與 LLM 從全文抽出的 disease / tissue_affected，可當 ground truth。

流程：
  1. 用 utils/hgvs_resolver 把 hgvs（c./p./g.）離線換算成 GRCh38 座標
     （不需連網、不需 ClinVar；同一篇論文若同時有 c. 與 p.，會互相驗證轉錄本編號）。
  2. 同一顆變異（解析後的 variant_id）跨列合併，disease / tissue 取多數。
  3. 對每顆變異跑 MechanisticInference（預設 rule-based，--llm 走實際模型）。
  4. 用 utils/disease_normalizer 把疾病名映射成 canonical label 後比對：
       exact / compatible（一般↔特定）/ differential / same_tissue / miss
  5. 輸出明細 CSV + 準確率摘要（另外分出「HGMD 完全沒有」的純推理子集）。

每跑完一顆變異就立刻 append 到輸出 CSV，並且預設 resume：再次執行時會跳過
CSV 裡已經有的變異，所以 LLM 全量測試中途斷線 / 中斷都能接續（--no-resume 可重跑）。
只想重新產生摘要時用 --report-only。

用法：
    conda run -n gm4 python test_inference_pubtator_only.py                  # 全部，rule-based
    conda run -n gm4 python test_inference_pubtator_only.py --limit 100
    conda run -n gm4 python test_inference_pubtator_only.py --exclude-hgmd position
    conda run -n gm4 python test_inference_pubtator_only.py --llm            # 用 Gemma 推理（需 GPU）
    conda run -n gm4 python test_inference_pubtator_only.py --llm --report-only
"""

from __future__ import annotations

import argparse
import csv
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd

import config
from utils.disease_normalizer import (
    canonical_labels,
    compare_disease,
    is_ttn_relevant,
    is_unspecified,
    normalize_tissue,
    tissue_match,
)
from utils.hgvs_resolver import HgvsResolveError
from utils.mechanistic_inference import MechanisticInference

PROJECT_ROOT = Path(__file__).parent
DEFAULT_CSV = PROJECT_ROOT / "pubtator_only_variants.csv"
DEFAULT_OUT = PROJECT_ROOT / "outputs" / "inference_pubtator_only_validation.csv"
DEFAULT_OUT_LLM = PROJECT_ROOT / "outputs" / "inference_pubtator_only_validation_llm.csv"

# 輸出 CSV 的欄位（固定順序，才能安全 append / resume）
FIELDNAMES = [
    "variant_id", "pos", "input_hgvs", "transcript", "hgvs_c", "hgvs_p",
    "resolved_by", "resolve_confidence", "band", "vclass", "vgroup",
    "n_analogs", "nearest_analog_bp", "hgmd_exact", "hgmd_same_position",
    "pure_inference", "n_pubtator_rows", "pmids", "truth_disease", "truth_labels",
    "truth_tissue", "pred_disease", "pred_labels", "pred_tissue",
    "pred_confidence", "pred_source", "differentials", "disease_eval",
    "tissue_eval", "pred_vague", "rationale",
]


def _clean(v) -> Optional[str]:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).strip()
    return s if s and s.lower() != "nan" else None


def _mode(values: List[Optional[str]]) -> Optional[str]:
    vals = [v for v in values if v]
    return Counter(vals).most_common(1)[0][0] if vals else None


def _best_disease(values: List[Optional[str]]) -> Optional[str]:
    """多列合併時，優先取有明確病名者（忽略 Not specified），再取最常見。"""
    named = [v for v in values if v and not is_unspecified(v) and canonical_labels(v)]
    return _mode(named) if named else None


# ── 從 CSV 建立「解析後的變異 → 文獻標註」清單 ────────────────────────────────


def collect_variants(
    csv_path: Path, mi: MechanisticInference, use_rsid: bool = False
) -> tuple:
    df = pd.read_csv(csv_path)
    print(f"讀入 {csv_path.name}：{len(df)} 列")

    # 同一篇論文裡的 p.HGVS：用來交叉驗證 c.HGVS 該用哪個轉錄本編號
    p_hint: Dict[tuple, str] = {}
    for r in df.itertuples():
        h = _clean(getattr(r, "norm_hgvs", None)) or _clean(getattr(r, "hgvs", None))
        rsid = _clean(getattr(r, "rsid", None))
        if h and h.lower().startswith("p.") and rsid:
            p_hint.setdefault((getattr(r, "pmid", None), rsid.split("(")[0]), h)

    resolve_cache: Dict[tuple, object] = {}
    groups: Dict[str, Dict] = {}
    fail_reasons: Counter = Counter()
    fail_examples: List[str] = []
    n_no_id = 0

    for r in df.itertuples():
        hgvs = _clean(getattr(r, "norm_hgvs", None)) or _clean(getattr(r, "hgvs", None))
        rsid = _clean(getattr(r, "rsid", None))
        rsid = rsid.split("(")[0] if rsid else None
        queries = [q for q in (hgvs, _clean(getattr(r, "variant_text", None))) if q]
        if use_rsid and rsid:
            queries.append(rsid)
        if not queries:
            n_no_id += 1
            continue

        hint = p_hint.get((getattr(r, "pmid", None), rsid)) if rsid else None
        if hint and hgvs and hgvs.lower().startswith("p."):
            hint = None

        key = (tuple(queries), hint)
        if key not in resolve_cache:
            resolved = None
            first_err = ""
            for q in queries:
                try:
                    resolved = mi.resolver.resolve(q, hgvs_p=hint)
                    break
                except HgvsResolveError as e:
                    first_err = first_err or str(e)
                except Exception as e:  # noqa: BLE001
                    first_err = first_err or f"{type(e).__name__}: {e}"
            resolve_cache[key] = resolved or first_err
        got = resolve_cache[key]
        if isinstance(got, str):
            fail_reasons[_reason_bucket(got)] += 1
            if len(fail_examples) < 12:
                fail_examples.append(f"{queries[0]} → {got}")
            continue

        g = groups.setdefault(got.variant_id, {"resolved": got, "rows": []})
        g["rows"].append({
            "pmid": _clean(getattr(r, "pmid", None)),
            "year": _clean(getattr(r, "year", None)),
            "hgvs": hgvs,
            "disease": _clean(getattr(r, "disease", None)),
            "tissue": _clean(getattr(r, "tissue_affected", None)),
            "pathogenicity": _clean(getattr(r, "pathogenicity", None)),
            # 原文的證據句。pathogenicity 欄是 LLM 抽的，會把 "VUS-LP" 記成
            # "Likely pathogenic"，得靠原句才看得出那其實還是 VUS。
            "evidence": _clean(getattr(r, "evidence_sentence", None)),
        })

    n_merged = _merge_protein_only_groups(groups)
    print(
        f"解析成功：{len(groups)} 顆不同變異"
        f"（無法解析 {sum(fail_reasons.values())} 列，無識別碼 {n_no_id} 列，"
        f"p.HGVS 併入同 codon 的 c.HGVS 變異 {n_merged} 組）"
    )
    if fail_reasons:
        print("  無法解析的原因分佈：")
        for reason, n in fail_reasons.most_common():
            print(f"    {n:>5}  {reason}")
        print("  範例：")
        for e in fail_examples[:6]:
            print(f"    - {e}")
    return groups, fail_reasons


def _merge_protein_only_groups(groups: Dict[str, Dict]) -> int:
    """
    只給 p.HGVS 的變異只能定位到 codon 起點（沒有 ref/alt），會和同一顆變異的
    c.HGVS 版本各自成組。若同一 codon（pos, pos-1, pos-2；TTN 負鏈）已有帶
    ref/alt 且類型相同的變異，就併進去，避免同一顆變異被算成兩筆。
    """
    with_alleles: Dict[int, List[str]] = defaultdict(list)
    for vid, g in groups.items():
        rv = g["resolved"]
        if rv.ref and rv.alt:
            with_alleles[rv.pos].append(vid)

    merged = 0
    for vid in list(groups):
        rv = groups[vid]["resolved"]
        if rv.ref and rv.alt:
            continue
        target = None
        for delta in (0, -1, -2):
            for cand in with_alleles.get(rv.pos + delta, []):
                if groups[cand]["resolved"].vclass == rv.vclass:
                    target = cand
                    break
            if target:
                break
        if target:
            groups[target]["rows"].extend(groups.pop(vid)["rows"])
            merged += 1
    return merged


def _iteration_order(groups: Dict[str, Dict], limit: int) -> List[tuple]:
    """依 genomic 位置排序；有 --limit 時改用等距抽樣，避免只測到基因一端的 band。"""
    ordered = sorted(groups.items(), key=lambda kv: kv[1]["resolved"].pos)
    if not limit or limit >= len(ordered):
        return ordered
    stride = max(1, len(ordered) // (limit * 3))
    return [ordered[i] for offset in range(stride) for i in range(offset, len(ordered), stride)]


def _reason_bucket(msg: str) -> str:
    m = msg.lower()
    if "grch37" in m:
        return "g.HGVS 是 GRCh37 座標"
    if "不在 grch38" in m:
        return "座標不在 TTN 區間"
    if "沒有轉錄本在 aa" in m:
        return "p.HGVS 的胺基酸與所有轉錄本都不符"
    if "參考鹼基相符" in m:
        return "c.HGVS 參考鹼基與所有轉錄本都不符"
    if "無法定位" in m:
        return "c.HGVS 位置無法定位（intronic 邊界不符）"
    if "不認得" in m or "無法解析" in m:
        return "無法解析的寫法（rsID / 自由文字）"
    if "rsid" in m:
        return "rsID（需 --rsid-lookup）"
    return "其他"


# ── 主流程 ───────────────────────────────────────────────────────────────────


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--out", type=Path, default=None,
                    help=f"預設 {DEFAULT_OUT.name}（--llm 時為 {DEFAULT_OUT_LLM.name}）")
    ap.add_argument("--limit", type=int, default=0, help="最多測幾顆變異（0=全部）")
    ap.add_argument("--llm", action="store_true", help="用實際 LLM 推理（需 GPU/vLLM）")
    ap.add_argument(
        "--exclude-hgmd",
        choices=["none", "exact", "position"],
        default="none",
        help="排除 HGMD 已收錄的變異：exact=同 pos+ref+alt；position=同位置（最嚴格的純推理）",
    )
    ap.add_argument("--include-nonttn", action="store_true",
                    help="連癌症等與 TTN 機制無關的表型也納入評分")
    ap.add_argument("--rsid-lookup", action="store_true",
                    help="沒有 HGVS 時用 rsID 查 dbSNP（需連網，結果會快取）")
    ap.add_argument("--batch", type=int, default=0,
                    help="一次送幾個 prompt 給 LLM（0=自動：--llm 用 16，rule-based 用 64）")
    ap.add_argument("--verbose", action="store_true", help="逐筆列印")
    ap.add_argument("--no-resume", dest="resume", action="store_false",
                    help="不續跑，直接覆寫輸出 CSV")
    ap.add_argument("--report-only", action="store_true",
                    help="不跑推理，只用現有 CSV 重新產生摘要")
    ap.add_argument("--rescore", action="store_true",
                    help="用目前的疾病正規化規則重評現有 CSV（會覆寫評分欄位）並重出摘要")
    args = ap.parse_args()
    if args.out is None:
        args.out = DEFAULT_OUT_LLM if args.llm else DEFAULT_OUT
    if not args.batch:
        args.batch = 16 if args.llm else 64

    if args.report_only or args.rescore:
        if not args.out.exists():
            raise SystemExit(f"找不到 {args.out}")
        df = pd.read_csv(args.out)
        if args.rescore:
            before = df["disease_eval"].copy()
            df = rescore(df)
            changed = int((before.fillna("") != df["disease_eval"].fillna("")).sum())
            df.to_csv(args.out, index=False)
            print(f"已用現行規則重評 {len(df)} 筆，{changed} 筆評分改變 -> {args.out}")
        report(df, Counter(), args)
        return

    mi = MechanisticInference(
        kb_path=config.HGMD_KB_PATH,
        band_edges=config.TTN_BAND_EDGES,
        band_labels=config.TTN_BAND_LABELS,
        top_k=config.INFERENCE_TOP_K_ANALOGS,
        transcript_interval_path=config.TRANSCRIPT_INTERVAL_PATH,
        transcript_model_path=config.TTN_TRANSCRIPT_MODEL_PATH,
        reference_fasta_path=config.REFERENCE_GENOME_PATH,
        allow_network=args.rsid_lookup,
    )

    groups, _ = collect_variants(args.csv, mi, use_rsid=args.rsid_lookup)

    done: Set[str] = set()
    if args.resume and args.out.exists():
        try:
            done = set(pd.read_csv(args.out)["variant_id"].astype(str))
            print(f"續跑：{args.out.name} 已有 {len(done)} 顆變異，將跳過")
        except Exception as e:  # noqa: BLE001
            print(f"讀取既有輸出失敗（將重新開始）：{e}")
    elif args.out.exists():
        args.out.unlink()

    extractor = None
    if args.llm:
        from utils.local_clinical_extractor import LocalClinicalExtractor

        # LocalClinicalExtractor 的預設值是 Llama-3.2-3B / TP=2，必須照 config 傳入。
        # 模型在 __init__ 就載入（沒有 load_model()）。
        extractor = LocalClinicalExtractor(
            backend=config.LOCAL_LLM_BACKEND,
            model_name=config.LOCAL_LLM_MODEL,
            tensor_parallel_size=config.LOCAL_LLM_TENSOR_PARALLEL,
            max_model_len=config.LOCAL_LLM_MAX_MODEL_LEN,
            max_context_length=config.LOCAL_LLM_MAX_CONTEXT_LENGTH,
        )
        print(
            f"LLM 已就緒：{config.LOCAL_LLM_MODEL}（TP={config.LOCAL_LLM_TENSOR_PARALLEL}）",
            flush=True,
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    # 續跑時沿用檔案裡既有的表頭：若程式後來多加了欄位，用新表頭續寫會讓
    # 新舊列的欄位對不上，整個 CSV 就毀了。
    fieldnames = FIELDNAMES
    if done:
        try:
            with open(args.out, newline="", encoding="utf-8") as fh:
                existing = next(csv.reader(fh))
            if existing and existing != FIELDNAMES:
                fieldnames = existing
                print(f"沿用既有表頭（{len(existing)} 欄）以免續寫錯位")
        except Exception as e:  # noqa: BLE001
            print(f"讀取既有表頭失敗（用預設欄位）：{e}")

    writer_fh = open(args.out, "a", newline="", encoding="utf-8")
    writer = csv.DictWriter(writer_fh, fieldnames=fieldnames, extrasaction="ignore")
    if not done:
        writer.writeheader()
        writer_fh.flush()

    # 先把要評分的變異挑出來（座標解析 / ground truth 過濾都很便宜），
    # 再分批送進推理，這樣 LLM 模式可以一次 generate 多個 prompt。
    skipped = Counter()
    todo: List[Dict] = []
    for vid, g in _iteration_order(groups, args.limit):
        rv = g["resolved"]
        if vid in done:
            skipped["已在輸出 CSV（續跑跳過）"] += 1
            continue
        truth_disease = _best_disease([r["disease"] for r in g["rows"]])
        if not truth_disease:
            skipped["ground truth 沒有明確疾病"] += 1
            continue
        truth_labels = canonical_labels(truth_disease)
        if not args.include_nonttn and not is_ttn_relevant(truth_labels):
            skipped["非 TTN 機制表型（癌症等）"] += 1
            continue

        hgmd_hit = mi.hgmd_lookup(rv.pos, rv.ref, rv.alt)
        if args.exclude_hgmd == "exact" and hgmd_hit["exact"]:
            skipped["HGMD 已收錄同一變異"] += 1
            continue
        if args.exclude_hgmd == "position" and hgmd_hit["same_position"]:
            skipped["HGMD 有同位置變異"] += 1
            continue

        todo.append({
            "vid": vid, "resolved": rv, "group": g, "hgmd_hit": hgmd_hit,
            "truth_disease": truth_disease, "truth_labels": truth_labels,
            "truth_tissue": normalize_tissue(
                _mode([r["tissue"] for r in g["rows"]]) or truth_disease
            ),
        })
        if args.limit and len(todo) >= args.limit:
            break

    print(f"待推理變異：{len(todo)}（batch={args.batch}）", flush=True)

    rows: List[Dict] = []
    t0 = time.time()
    for start in range(0, len(todo), args.batch):
        chunk = todo[start: start + args.batch]
        try:
            results = mi.infer_batch(extractor, [c["resolved"].to_variant_info() for c in chunk])
        except Exception as e:  # noqa: BLE001  一批失敗不該讓整輪中斷
            print(f"  ! batch @{start} 推理失敗，跳過：{type(e).__name__}: {e}", flush=True)
            skipped["推理失敗"] += len(chunk)
            continue

        for item, res in zip(chunk, results):
            feat, inf = res["features"], res["inference"]
            rv, hgmd_hit = item["resolved"], item["hgmd_hit"]
            pred_disease = inf.get("inferred_disease")
            differentials = inf.get("differential_diagnoses") or []
            cmp_res = compare_disease(pred_disease, item["truth_disease"], differentials)
            pred_tissue = normalize_tissue(inf.get("inferred_tissue"))
            analogs = feat.get("analogs", [])

            row = {
                "variant_id": item["vid"],
                "pos": rv.pos,
                "input_hgvs": rv.input,
                "transcript": rv.transcript,
                "hgvs_c": rv.hgvs_c,
                "hgvs_p": rv.hgvs_p or feat.get("hgvs_p"),
                "resolved_by": rv.resolved_by,
                "resolve_confidence": rv.confidence,
                "band": feat.get("band"),
                "vclass": feat.get("vclass"),
                "vgroup": feat.get("vgroup"),
                "n_analogs": len(analogs),
                "nearest_analog_bp": min(
                    [a.get("distance", 10**9) for a in analogs], default=None
                ),
                "hgmd_exact": bool(hgmd_hit["exact"]),
                "hgmd_same_position": bool(hgmd_hit["same_position"]),
                "pure_inference": not hgmd_hit["same_position"],
                "n_pubtator_rows": len(item["group"]["rows"]),
                "pmids": "; ".join(
                    sorted({r["pmid"] for r in item["group"]["rows"] if r["pmid"]})
                ),
                "truth_disease": item["truth_disease"],
                "truth_labels": ",".join(item["truth_labels"]),
                "truth_tissue": item["truth_tissue"],
                "pred_disease": pred_disease,
                "pred_labels": ",".join(cmp_res["pred_labels"]),
                "pred_tissue": pred_tissue,
                "pred_confidence": inf.get("confidence"),
                "pred_source": inf.get("source"),
                "differentials": "; ".join(str(d) for d in differentials),
                "disease_eval": cmp_res["level"],
                "tissue_eval": tissue_match(pred_tissue, item["truth_tissue"]),
                "pred_vague": cmp_res["pred_vague"],
                "rationale": inf.get("rationale"),
            }
            rows.append(row)
            # 逐筆落地 + flush：中斷或斷線都不會丟掉已算好的結果
            writer.writerow(row)
            writer_fh.flush()

            if args.verbose:
                print(
                    f"  {item['vid']:<24} {str(feat.get('band')):<7}"
                    f"{str(feat.get('vclass')):<11}"
                    f"PT:{item['truth_disease'][:30]:<32}-> "
                    f"{str(pred_disease)[:28]:<30}[{cmp_res['level']}]",
                    flush=True,
                )

        elapsed = time.time() - t0
        rate = elapsed / max(len(rows), 1)
        eta = rate * (len(todo) - len(rows)) / 60
        hit = sum(1 for r in rows if r["disease_eval"] in ("exact", "compatible"))
        print(
            f"[{len(rows):>4}/{len(todo)}] {rate:.1f}s/變異 | "
            f"目前疾病命中 {hit}/{len(rows)} ({hit / max(len(rows), 1):.0%}) | "
            f"剩餘約 {eta:.0f} 分鐘",
            flush=True,
        )

    writer_fh.close()
    if not args.out.exists():
        print("沒有可評分的變異。")
        return

    df = pd.read_csv(args.out)  # 含續跑前已完成的部分
    report(df, skipped, args)
    print(f"\n明細已寫入：{args.out}（本次新增 {len(rows)} 筆，總計 {len(df)} 筆）")


# ── 重新評分 ─────────────────────────────────────────────────────────────────


def rescore(df: pd.DataFrame) -> pd.DataFrame:
    """用目前的 disease_normalizer 重算評分欄位。

    明細 CSV 存的是 truth_disease / pred_disease / differentials 原文，所以
    正規化規則修好之後可以直接重評，不必再跑一次（LLM 推理很貴）。
    """
    out = df.copy()
    cols = {"truth_labels": [], "pred_labels": [], "truth_tissue": [],
            "pred_tissue": [], "disease_eval": [], "tissue_eval": [], "pred_vague": []}
    for r in df.itertuples():
        diffs = [d.strip() for d in str(getattr(r, "differentials", "") or "").split(";") if d.strip()]
        cmp_res = compare_disease(_clean(r.pred_disease), _clean(r.truth_disease), diffs)
        # 組織沿用推理當下的輸出（LLM 有自己的 inferred_tissue），
        # 只有欄位缺值時才退回從疾病名推。
        pred_tissue = normalize_tissue(_clean(getattr(r, "pred_tissue", None))
                                       or _clean(r.pred_disease))
        truth_tissue = normalize_tissue(_clean(getattr(r, "truth_tissue", None))
                                        or _clean(r.truth_disease))
        cols["truth_labels"].append(",".join(cmp_res["truth_labels"]))
        cols["pred_labels"].append(",".join(cmp_res["pred_labels"]))
        cols["truth_tissue"].append(truth_tissue)
        cols["pred_tissue"].append(pred_tissue)
        cols["disease_eval"].append(cmp_res["level"])
        cols["tissue_eval"].append(tissue_match(pred_tissue, truth_tissue))
        cols["pred_vague"].append(cmp_res["pred_vague"])
    for k, v in cols.items():
        out[k] = v
    return out[[c for c in FIELDNAMES if c in out.columns]]


# ── 摘要 ─────────────────────────────────────────────────────────────────────

_DISEASE_ORDER = ["exact", "compatible", "differential", "same_tissue", "miss", "n/a"]
_TISSUE_ORDER = ["exact", "partial", "miss", "n/a"]


def _dist(series: pd.Series, order: List[str]) -> str:
    n = len(series)
    counts = series.value_counts()
    parts = [
        f"{k}={int(counts.get(k, 0))} ({counts.get(k, 0) / n:.0%})"
        for k in order if counts.get(k, 0)
    ]
    return "  ".join(parts)


def _summary_block(df: pd.DataFrame, title: str) -> None:
    n = len(df)
    if not n:
        return
    d = df["disease_eval"]
    hit = int(((d == "exact") | (d == "compatible")).sum())
    with_diff = hit + int((d == "differential").sum())
    t = df["tissue_eval"]
    t_ok = int(((t == "exact") | (t == "partial")).sum())
    t_scored = int((t != "n/a").sum())
    print(f"\n{title}（n={n}）")
    print(f"  disease：{_dist(d, _DISEASE_ORDER)}")
    print(f"    疾病命中（exact+compatible）：{hit}/{n} ({hit / n:.0%})")
    if "pred_vague" in df.columns:
        vague = df["pred_vague"].fillna(False).astype(bool)
        specific = int((((d == "exact") | (d == "compatible")) & ~vague).sum())
        print(f"    其中具體命中（排除只答 Titinopathy/Cardiomyopathy 這類籠統詞）："
              f"{specific}/{n} ({specific / n:.0%})")
    print(f"    含 differential 命中：{with_diff}/{n} ({with_diff / n:.0%})")
    print(f"  tissue：{_dist(t, _TISSUE_ORDER)}")
    if t_scored:
        print(f"    組織相容（exact+partial）：{t_ok}/{t_scored} ({t_ok / t_scored:.0%})")


def report(df: pd.DataFrame, skipped: Counter, args) -> None:
    print("\n" + "=" * 78)
    print(f"評分變異數：{len(df)}   （推理模式：{'LLM' if args.llm else 'rule-based'}）")
    if skipped:
        print("略過：" + "，".join(f"{k} {v}" for k, v in skipped.most_common()))

    _summary_block(df, "整體")
    _summary_block(df[df["pure_inference"]], "純推理子集（HGMD 連同位置都沒有）")
    _summary_block(df[~df["pure_inference"]], "HGMD 有同位置變異的子集")

    print("\n分 band（疾病 exact+compatible 命中率）")
    for band, g in df.groupby("band"):
        hit = int(g["disease_eval"].isin(["exact", "compatible"]).sum())
        print(f"  {str(band):<8} {hit:>4}/{len(g):<4} ({hit / len(g):.0%})")

    print("\n分變異類型")
    for vg, g in df.groupby("vclass"):
        hit = int(g["disease_eval"].isin(["exact", "compatible"]).sum())
        print(f"  {str(vg):<12} {hit:>4}/{len(g):<4} ({hit / len(g):.0%})")

    print("\n分 HGVS 解析來源")
    for src, g in df.groupby(df["resolved_by"].str.split(":").str[0]):
        hit = int(g["disease_eval"].isin(["exact", "compatible"]).sum())
        print(f"  {str(src):<16} {hit:>4}/{len(g):<4} ({hit / len(g):.0%})")

    print("\n最常見的 ground truth vs 推論（前 12 組）")
    pairs = Counter(
        (r.truth_labels or "?", r.pred_labels or "?", r.disease_eval) for r in df.itertuples()
    )
    for (t, p, lvl), n in pairs.most_common(12):
        print(f"  {n:>4}  truth={t[:26]:<28} pred={p[:26]:<28} [{lvl}]")


if __name__ == "__main__":
    main()
