#!/usr/bin/env python3
"""把 NTUH 回饋後的報告改版套用到 outputs/ 裡既有的報告，不重跑 pipeline。

每份報告會：
  - 換成新版樣式（字型統一、資料庫卡片對齊、推論標示、側邊目錄）並加上目錄程式
  - 標題、header、Variant Summary 改用基因組正股 ID，加上轉錄本 c. 寫法
  - Evo2 區塊加上分數說明與 cutoff
  - 重畫 isoform 示意圖（心肌／骨骼肌 isoform 框、exon 49）
  - Track B 推論區塊與文獻 LLM 欄位標成推論、沒直接提到變異的文獻改灰
  - 加上遺傳模式（AD / AR）推論區塊（--llm 時用 Gemma 判斷，否則規則式）

可以重複執行：每一步都是「找到舊版就換成新版」，已經是新版的會原樣保留。

    python backfill_report_ui.py                    # 全部報告，遺傳模式用規則式
    python backfill_report_ui.py --llm              # 遺傳模式改用 LLM（需要 vLLM 的 GPU 記憶體）
    python backfill_report_ui.py --limit 20 --dry-run
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import os
import re
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

OUTPUT_DIR = ROOT / "outputs"
STEM_RE = re.compile(r"^(2)-(\d+)-([ACGT])-([ACGT])_(.+)$")

RE_STYLE = re.compile(r"<style>.*?</style>", re.S)
RE_TITLE = re.compile(r"<title>.*?</title>", re.S)
RE_HEADER_VARIANT = re.compile(
    r'Variant: (?:<a href="([^"]*)"[^>]*>)?[^<]*?(?:</a>)?<br>(?:<span class="coding-strand">.*?</span><br>)?', re.S
)
RE_SUMMARY = re.compile(r'<section class="section"><h2>Variant Summary</h2>.*?</section>', re.S)
RE_EVO2 = re.compile(r'<section class="section"><h2>Evo2 Prediction.*?</section>', re.S)
RE_PRED = re.compile(r'<div class="prediction-result">([A-Z]+)</div>')
RE_SCORE = re.compile(r'<div class="score-value">(-?[\d.]+)</div>')
RE_IMAGE = re.compile(r'<section class="section"><h2>Protein Domain Localization</h2>.*?</section>', re.S)
RE_INHERIT = re.compile(r'\s*<section class="section inference-section inherit-section">.*?</section>', re.S)
RE_OLD_INF_NOTE = re.compile(r'<div class="inference-note">\s*No published article was found.*?</div>', re.S)
RE_OLD_TOC = re.compile(r"\s*<script>\(function \(\) \{\s*var secs = document\.querySelectorAll.*?\}\)\(\);</script>", re.S)
RE_ARTICLE = re.compile(r'<div class="pubmed-article(?: indirect)?" style="')
RE_STATS_INDIRECT = re.compile(
    r'<div style="margin-bottom: (\d+)px;">(\s*<h4[^>]*>Articles that do NOT directly mention)')
RE_EXTRACT_SLOT = re.compile(
    r'(<div class="article-info-label">Inheritance</div><div>[^<]*</div></div>\s*</div>)(?!\s*<div class="extract-note">)')
RE_TRACKB_DISEASE = re.compile(r'Predicted Disease \(inference\)</div><div[^>]*>(.*?)</div>', re.S)
RE_TRACKB_TYPE = re.compile(r'Variant Type</div><div[^>]*>(.*?)</div>', re.S)
RE_ART_FIELD = re.compile(
    r'<div class="article-meta">[^<]*\| PMID: (\d+)</div>.*?Disease</div><div>(.*?)</div></div>.*?'
    r'Inheritance</div><div>(.*?)</div>', re.S)

RE_AI_LABEL = re.compile(
    r'<div style="font-weight: bold; color: #[0-9A-Fa-f]{6}; margin-bottom: 5px;">AI Analysis</div>'
    r'|<div class="llm-extract-label">AI Analysis <span class="ai-tag"[^>]*>[^<]*</span></div>')
RE_REASONING = re.compile(r'AI Analysis.*?<div style="font-style: italic;[^"]*">(.*?)</div>', re.S)
RE_ANALOG_NAN = re.compile(r'(<tr>\s*<td>(?:<a [^>]*>)?)nan((?:</a>)?</td>)')
EXTRACT_NOTE = ('<div class="extract-note">Disease, tissue, age and inheritance were extracted from this article by '
                'the language model (LLM); check the article before relying on them.</div>')


def clinvar_names() -> Dict[str, Dict[str, str]]:
    """data/clinvar_ttn.csv 的 HGVS 名稱（用來判斷變異類型），以內部 ID 為鍵。"""
    from utils import variant_ids

    out: Dict[str, Dict[str, str]] = {}
    path = ROOT / "data" / "clinvar_ttn.csv"
    if not path.exists():
        return out
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            try:
                vid = variant_ids.internal_id(row["input"])
            except variant_ids.VariantIdError:
                continue
            name = row.get("Name") or ""
            m = re.search(r"(c\.[^\s(]+)", name)
            p = re.search(r"\((p\.[^)]+)\)", name)
            out[vid] = {"hgvs_nucleotide": m.group(1) if m else "", "hgvs_protein": p.group(1) if p else ""}
    return out


def mention_articles(src: str) -> List[Dict]:
    """直接提到變異的文獻（mention-true）的 PMID、疾病、inheritance，當作遺傳模式推論的文獻證據。"""
    out = []
    for block in src.split('<div class="pubmed-article')[1:]:
        if "mention-true" not in block.split('<div class="article-info">', 1)[0]:
            continue
        m = RE_ART_FIELD.search(block)
        if m:
            disease = re.sub(r"<[^>]+>", "", m.group(2)).replace("(inferred from article)", "").strip()
            out.append({"pmid": m.group(1), "clinical_info": {"disease": html.unescape(disease),
                                                             "inheritance": html.unescape(m.group(3).strip())}})
    return out


# ── 證據收集（主程序）──────────────────────────────────────────────────────

def collect_inheritance(paths: List[Path], use_llm: bool) -> Dict[str, Dict]:
    import main as pipeline
    from utils.inheritance_inference import InheritanceInference
    from utils.variant_db_lookup import get_shared_lookup

    mi = pipeline.build_mechanistic_engine()
    inh = InheritanceInference(mi, pipeline.INFERENCE_HGMD_CSV_PATH)
    db = get_shared_lookup(
        hgmd_csv=pipeline.HGMD_FULL_CSV_PATH, hgmd_kb=pipeline.HGMD_KB_PATH,
        pubtator_csv=pipeline.PUBTATOR_VARIANTS_CSV_PATH, cache_path=pipeline.PUBTATOR_COORD_CACHE_PATH,
        pubmind_json=pipeline.PUBMIND_JSON_PATH,
    )
    names = clinvar_names()
    items = []
    t0 = time.time()
    for i, path in enumerate(paths, 1):
        m = STEM_RE.match(path.stem)
        if not m:
            continue
        chrom, pos, ref, alt, _ = m.groups()
        vid = f"{chrom}-{pos}-{ref}-{alt}"
        src = path.read_text(encoding="utf-8", errors="replace")
        variant_info = {"variant_id": vid, "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt}
        clinvar_info = names.get(vid) or {}
        tb_type = RE_TRACKB_TYPE.search(src)
        if not clinvar_info and tb_type:  # Track B 區塊已經判過類型（例：「nonsense (truncating)」）
            t = re.match(r"\s*(\w+)\s*\((\w+)\)", html.unescape(tb_type.group(1)))
            if t and t.group(2) != "unknown":
                variant_info.update(hgvs_vclass=t.group(1), hgvs_vgroup=t.group(2))
        features = mi.build_features(variant_info, clinvar_info or None, None)
        articles = mention_articles(src)
        tb = RE_TRACKB_DISEASE.search(src)
        disease = html.unescape(re.sub(r"<[^>]+>", "", tb.group(1))).strip() if tb else ""
        if not disease or disease.lower() == "uncertain":
            disease = pipeline._main_disease(None, articles) or ""
        dbann = db.lookup(chrom=chrom, pos=int(pos), ref=ref, alt=alt)
        ev = inh.collect(features, dbann, articles, disease or None)
        items.append((path.name, ev, inh.rule_predict(ev)))
        if i % 500 == 0:
            print(f"  evidence {i}/{len(paths)} ({time.time() - t0:.0f}s)", flush=True)

    results: Dict[str, Dict] = {}
    for name, ev, rule in items:
        from utils.inheritance_inference import _public_evidence

        results[name] = {"inheritance": rule["inheritance"], "confidence": rule["confidence"],
                         "rationale": "Rule-based: " + "; ".join(rule["reasons"]) + ".",
                         "source": "rule-based", "rule": rule, "evidence": _public_evidence(ev)}
    if use_llm and items:
        extractor = pipeline.get_shared_extractor(pipeline.LOCAL_LLM_BACKEND, pipeline.LOCAL_LLM_MODEL)
        is_medgemma = getattr(extractor, "is_medgemma", True)
        prompts = [inh.build_prompt(ev, rule, is_medgemma) for _, ev, rule in items]
        t0 = time.time()
        raws = inh.generate(extractor, prompts)
        print(f"  LLM: {len(prompts)} prompts in {time.time() - t0:.0f}s", flush=True)
        label = "LLM (" + str(getattr(extractor, "model_name", "") or pipeline.LOCAL_LLM_MODEL).split("/")[-1] + ")"
        for (name, _, _), raw in zip(items, raws):
            parsed = inh.parse(raw)
            if parsed:
                results[name].update(parsed, source=label)
            else:
                results[name]["source"] = "rule-based (LLM output unparsable)"
    return results


# ── 單份報告改版（子程序）────────────────────────────────────────────────

def transform(path_str: str, inheritance: Optional[Dict], rerender: bool, dry_run: bool) -> str:
    import matplotlib

    matplotlib.use("Agg")
    from config import REPORT_TITLE
    from utils.html_report import HTMLReportGenerator, INFERENCE_WARNING_HTML, REPORT_TOC_JS

    path = Path(path_str)
    m = STEM_RE.match(path.stem)
    if not m:
        return "skip"
    chrom, pos, ref, alt, _ = m.groups()
    vid = f"{chrom}-{pos}-{ref}-{alt}"
    gen = HTMLReportGenerator()
    vi = gen._display({"variant_id": vid, "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt})
    src = path.read_text(encoding="utf-8", errors="replace")
    out = src

    out = RE_STYLE.sub(lambda _: "<style>\n        " + gen._get_css() + "\n    </style>", out, count=1)
    out = RE_TITLE.sub(lambda _: f"<title>{REPORT_TITLE} - {vi['display_id']}</title>", out, count=1)

    def _header(mm: re.Match) -> str:
        url = mm.group(1)
        link = (f'<a href="{url}" target="_blank" rel="noopener" class="clinvar-link">{vi["display_id"]}</a>'
                if url else vi["display_id"])
        coding = (f'<span class="coding-strand">TTN coding strand: {ref}&gt;{alt}'
                  + (f' · {html.escape(vi["hgvs_c_display"])}' if vi.get("hgvs_c_display") else "")
                  + "</span><br>")
        return f"Variant: {link}<br>{coding}"

    out = RE_HEADER_VARIANT.sub(_header, out, count=1)

    evo2_block = RE_EVO2.search(out)
    evo2 = None
    if evo2_block:
        pred = RE_PRED.search(evo2_block.group(0))
        scores = RE_SCORE.findall(evo2_block.group(0))
        if pred and len(scores) >= 3:
            evo2 = {"success": True, "prediction": pred.group(1).lower(), "delta_score": float(scores[0]),
                    "ref_score": float(scores[1]), "var_score": float(scores[2])}
            out = out.replace(evo2_block.group(0), gen.evo2_section_html(
                evo2["prediction"], evo2["delta_score"], evo2["ref_score"], evo2["var_score"]), 1)
    out = RE_SUMMARY.sub(lambda _: gen._generate_summary_section(vi, evo2), out, count=1)

    if rerender and RE_IMAGE.search(out):
        from utils.image_generator import ImageGenerator
        import config

        ig = ImageGenerator()
        with tempfile.TemporaryDirectory() as tmp:
            ig.output_dir = Path(tmp)
            img = ig.generate_titin_schematic(dict(vi), str(config.TRANSCRIPT_INTERVAL_PATH))
            b64 = base64.b64encode(img.read_bytes()).decode("ascii")
        out = RE_IMAGE.sub(lambda _: gen._generate_image_section(b64, vi), out, count=1)

    # Track B：標成推論
    out = out.replace("INFERENCE · NO DIRECT LITERATURE", "INFERENCE · NOT OBSERVED DATA")
    out = RE_OLD_INF_NOTE.sub(lambda _: INFERENCE_WARNING_HTML, out, count=1)
    for label in ("Predicted Tissue", "Confidence", "Differential"):
        out = out.replace(f'<div class="inference-card"><div class="label">{label}</div>',
                          f'<div class="inference-card llm"><div class="label">{label}</div>')
    out = out.replace('<div class="inference-card" style="border-left-width:6px;"><div class="label">Predicted Disease',
                      '<div class="inference-card llm" style="border-left-width:6px;"><div class="label">Predicted Disease')
    out = out.replace('<div class="inference-rationale">', '<div class="inference-rationale llm">')

    # 遺傳模式：放在文獻區塊前（沒有文獻區塊就放 footer 前）
    out = RE_INHERIT.sub("", out)
    if inheritance:
        section = gen._generate_inheritance_section(inheritance)
        anchor = out.find('<section class="section"><h2>Literature Review')
        if anchor < 0:
            anchor = out.find("<footer>")
        if anchor >= 0:
            out = out[:anchor] + section.strip() + "\n        " + out[anchor:]

    # 文獻：LLM 擷取欄位標示、沒直接提到變異的文獻改灰
    parts = RE_ARTICLE.split(out)
    if len(parts) > 1:
        rebuilt = [parts[0]]
        for part in parts[1:]:
            head = part[:3000]
            indirect = "mention-false" in head.split('<div class="article-info">', 1)[0]
            rebuilt.append(f'<div class="pubmed-article{" indirect" if indirect else ""}" style="')
            if not indirect:
                part = RE_EXTRACT_SLOT.sub(lambda mm: mm.group(1) + "\n                " + EXTRACT_NOTE, part, count=1)
            reasoning = RE_REASONING.search(part)
            part = RE_AI_LABEL.sub(
                lambda _: f'<div class="llm-extract-label">AI Analysis '
                          f'{HTMLReportGenerator.ai_tag(reasoning.group(1) if reasoning else "")}</div>', part, count=1)
            rebuilt.append(part)
        out = "".join(rebuilt)
    out = RE_STATS_INDIRECT.sub(r'<div class="stats-row stats-indirect" style="margin-bottom: \1px;">\2', out)
    out = RE_ANALOG_NAN.sub(r"\1—\2", out)
    out = out.replace('<div class="dbann-disease dbann-none">—</div>', '<div class="dbann-disease dbann-none">No record</div>')

    # 側邊目錄
    out = RE_OLD_TOC.sub("", out)
    out = out.replace("</body>", f"    <script>{REPORT_TOC_JS}</script>\n</body>", 1)

    if out == src:
        return "unchanged"
    if not dry_run:
        tmp = path.with_suffix(".html.tmp")
        tmp.write_text(out, encoding="utf-8")
        os.replace(tmp, path)
    return "updated"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", type=Path, default=OUTPUT_DIR)
    ap.add_argument("--llm", action="store_true", help="infer inheritance with the LLM (vLLM) instead of rules")
    ap.add_argument("--no-inheritance", action="store_true", help="skip the inheritance section")
    ap.add_argument("--no-image", action="store_true", help="keep the existing schematics")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", nargs="*", help="only these report file names")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 4))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    paths = sorted(p for p in args.outputs.glob("*.html") if STEM_RE.match(p.stem) and "smoketest" not in p.stem)
    if args.only:
        wanted = set(args.only)
        paths = [p for p in paths if p.name in wanted]
    if args.limit:
        paths = paths[: args.limit]
    print(f"[backfill-ui] {len(paths)} report(s)", flush=True)

    inheritance: Dict[str, Dict] = {}
    if not args.no_inheritance:
        print("[backfill-ui] collecting inheritance evidence…", flush=True)
        inheritance = collect_inheritance(paths, args.llm)
        cache = ROOT / "research" / "experiments" / "inheritance_backfill.json"
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps({k: {kk: v[kk] for kk in ("inheritance", "confidence", "source")}
                                     for k, v in inheritance.items()}, ensure_ascii=False, indent=0))

    counts: Dict[str, int] = {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = [pool.submit(transform, str(p), inheritance.get(p.name), not args.no_image, args.dry_run)
                for p in paths]
        for i, f in enumerate(futs, 1):
            try:
                r = f.result()
            except Exception as e:  # noqa: BLE001
                r = "error"
                print(f"  [error] {paths[i - 1].name}: {e}", flush=True)
            counts[r] = counts.get(r, 0) + 1
            if i % 250 == 0:
                print(f"  {i}/{len(paths)} ({time.time() - t0:.0f}s) {counts}", flush=True)
    print(f"[backfill-ui] done in {time.time() - t0:.0f}s: {counts}" + (" (dry run)" if args.dry_run else ""))
    return 0 if not counts.get("error") else 1


if __name__ == "__main__":
    raise SystemExit(main())
