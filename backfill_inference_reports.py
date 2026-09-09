#!/usr/bin/env python3
"""為既有的 *_no-llm.html 報告補上 Track B 類比推理區塊。

這些報告是「ClinVar/PubMed 都找不到任何文獻」的變異（Literature Review 區塊為空，
因此檔名後綴是 no-llm）。它們正是 Track B 該上場的案例，但產生當時管線還沒有
mechanistic inference，所以報告裡缺這一段。

為什麼不直接重跑 main.py：
  main.py 每個 process 都要重新載入 vLLM（3–5 分鐘）與 Evo2，512 顆光是模型載入
  就要 40 小時以上。這支腳本只載入一次 LLM，其餘可重用的結果（Evo2 分數、蛋白
  示意圖、ClinVar 標頭）直接沿用報告裡既有的內容，總時間約 1 小時。

作法是「注入」而非整份重建：保留原報告的標頭、Evo2、示意圖與 ClinVar 連結
（有 11 份含 ClinVar 資訊，重建會因為無法還原爬取結果而掉資料），
只補上推理區塊並把 <style> 換成目前版本（舊 CSS 沒有 Track B 的樣式）。

用法：
    conda run -n gm4 python backfill_inference_reports.py --limit 3 --dry-run
    conda run -n gm4 python backfill_inference_reports.py            # 全部
    conda run -n gm4 python backfill_inference_reports.py --force    # 連已有的也重算
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import config
from utils.html_report import HTMLReportGenerator
from utils.mechanistic_inference import MechanisticInference

RE_STYLE = re.compile(r"<style>.*?</style>", re.S)
RE_PREDICTION = re.compile(
    r'prediction-box prediction-(\w+)"><div class="prediction-result">([A-Z]+)</div>'
)
RE_SCORE = re.compile(r'<div class="score-value">(-?[\d.]+)</div>')
# Literature Review 區塊的開頭：推理區塊要插在它前面（與 main.py 的排序一致）
LIT_MARKER = '<section class="section"><h2>Literature Review'
INFERENCE_MARKER = "inference-section"


def parse_variant_id(stem: str) -> Optional[Dict]:
    """從檔名 stem（如 2-178527533-G-T）還原 variant_info。"""
    parts = stem.split("-")
    if len(parts) != 4:
        return None
    chrom, pos, ref, alt = parts
    if not pos.isdigit():
        return None
    return {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "variant_id": stem,
    }


def parse_evo2(text: str) -> Optional[Dict]:
    """從既有報告還原 evo2_result（Evo2 對同一顆變異是決定性的，可直接沿用）。"""
    pred = RE_PREDICTION.search(text)
    scores = RE_SCORE.findall(text)
    if not pred or len(scores) < 3:
        return None
    delta, ref_score, var_score = (float(s) for s in scores[:3])
    return {
        "success": True,
        "prediction": pred.group(1),
        "delta_score": delta,
        "ref_score": ref_score,
        "var_score": var_score,
        "confidence": abs(delta),
        "threshold": getattr(config, "PATHOGENIC_THRESHOLD", None),
    }


def inject(text: str, section_html: str, css: str) -> str:
    """把推理區塊插進 Literature Review 前面，並換成目前版本的 CSS。"""
    if LIT_MARKER not in text:
        raise ValueError("找不到 Literature Review 區塊，報告格式與預期不符")
    text = RE_STYLE.sub(lambda _: f"<style>\n{css}\n    </style>", text, count=1)
    return text.replace(LIT_MARKER, section_html.strip() + "\n        " + LIT_MARKER, 1)


def collect(outputs: Path, force: bool) -> Tuple[List[Dict], Dict[str, int]]:
    skipped: Dict[str, int] = {}

    def bump(k: str) -> None:
        skipped[k] = skipped.get(k, 0) + 1

    todo: List[Dict] = []
    for path in sorted(outputs.glob("*_no-llm.html")):
        text = path.read_text(encoding="utf-8", errors="ignore")
        if INFERENCE_MARKER in text and not force:
            bump("已有推理區塊")
            continue
        vi = parse_variant_id(path.name[: -len("_no-llm.html")])
        if vi is None:
            bump("檔名無法解析成變異")
            continue
        if LIT_MARKER not in text:
            bump("報告格式不符（找不到 Literature Review）")
            continue
        todo.append({"path": path, "text": text, "vi": vi, "evo2": parse_evo2(text)})
    return todo, skipped


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", type=Path, default=config.OUTPUT_DIR)
    ap.add_argument("--limit", type=int, default=0, help="最多處理幾份（0=全部）")
    ap.add_argument("--batch", type=int, default=16, help="一次送幾個 prompt 給 LLM")
    ap.add_argument("--force", action="store_true", help="連已有推理區塊的也重算")
    ap.add_argument("--dry-run", action="store_true",
                    help="不覆寫原檔，改寫到 <name>.preview.html 供檢查")
    ap.add_argument("--no-llm", action="store_true",
                    help="不載入 LLM，用規則式 fallback（快速驗證流程用）")
    args = ap.parse_args()

    todo, skipped = collect(args.outputs, args.force)
    if args.limit:
        todo = todo[: args.limit]
    print(f"待補推理區塊的報告：{len(todo)}")
    if skipped:
        print("略過：" + "，".join(f"{k} {v}" for k, v in skipped.items()))
    if not todo:
        return

    mi = MechanisticInference(
        kb_path=config.HGMD_KB_PATH,
        band_edges=config.TTN_BAND_EDGES,
        band_labels=config.TTN_BAND_LABELS,
        top_k=config.INFERENCE_TOP_K_ANALOGS,
        transcript_interval_path=config.TRANSCRIPT_INTERVAL_PATH,
        transcript_model_path=config.TTN_TRANSCRIPT_MODEL_PATH,
        reference_fasta_path=config.REFERENCE_GENOME_PATH,
        allow_network=False,
    )

    extractor = None
    if not args.no_llm:
        from utils.local_clinical_extractor import LocalClinicalExtractor

        extractor = LocalClinicalExtractor(
            backend=config.LOCAL_LLM_BACKEND,
            model_name=config.LOCAL_LLM_MODEL,
            tensor_parallel_size=config.LOCAL_LLM_TENSOR_PARALLEL,
            max_model_len=config.LOCAL_LLM_MAX_MODEL_LEN,
            max_context_length=config.LOCAL_LLM_MAX_CONTEXT_LENGTH,
        )
        print(f"LLM 已就緒：{config.LOCAL_LLM_MODEL}", flush=True)

    gen = HTMLReportGenerator()
    css = gen._get_css()

    done = 0
    failed = 0
    t0 = time.time()
    for start in range(0, len(todo), args.batch):
        chunk = todo[start: start + args.batch]
        try:
            results = mi.infer_batch(
                extractor,
                [c["vi"] for c in chunk],
                evo2_results=[c["evo2"] for c in chunk],
            )
        except Exception as e:  # noqa: BLE001  一批失敗不該讓整輪中斷
            print(f"  ! batch @{start} 推理失敗，跳過：{type(e).__name__}: {e}", flush=True)
            failed += len(chunk)
            continue

        for item, res in zip(chunk, results):
            path: Path = item["path"]
            try:
                section = gen._generate_inference_section(res)
                if not section.strip():
                    raise ValueError("推理區塊為空")
                out_text = inject(item["text"], section, css)
                target = (path.with_suffix(".preview.html") if args.dry_run else path)
                # 先寫暫存檔再 rename：中途失敗不會毀掉原本的報告
                tmp = target.with_suffix(".tmp")
                tmp.write_text(out_text, encoding="utf-8")
                os.replace(tmp, target)
                done += 1
            except Exception as e:  # noqa: BLE001
                print(f"  ! {path.name} 寫入失敗：{type(e).__name__}: {e}", flush=True)
                failed += 1

        elapsed = time.time() - t0
        rate = elapsed / max(done, 1)
        eta = rate * (len(todo) - done) / 60
        print(
            f"[{done:>4}/{len(todo)}] {rate:.1f}s/份 | 失敗 {failed} | 剩餘約 {eta:.0f} 分鐘",
            flush=True,
        )

    print(f"\n完成：{done} 份已補上推理區塊，{failed} 份失敗。")
    if args.dry_run:
        print("（--dry-run：結果寫在 *.preview.html，原報告未變更）")


if __name__ == "__main__":
    sys.exit(main())
