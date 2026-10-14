#!/usr/bin/env python3
"""Build a browsable static site from the per-variant HTML reports in outputs/.

The generated reports are self-contained single files (inline CSS + a ~470 KB
base64 schematic). This script turns them into a static site:

  site/
    index.html            searchable / sortable index of every report
    assets/               shared CSS, JS and the extracted schematic PNGs
    data/reports.json     the parsed metadata (also inlined into index.html)
    reports/*.html        the reports, rewritten to use the shared assets

index.html works when opened directly from the filesystem: the metadata is
inlined, so no fetch() (and therefore no local web server) is required.

Usage:
    python build_site.py                       # LLM reports only
    python build_site.py --include all         # every report in outputs/
    python build_site.py --site docs           # e.g. for GitHub Pages
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import html as html_mod
import json
import logging
import re
import shutil
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

from utils.pubmind_lookup import PubMindLookup
from utils.pubmind_lookup import disease_counts as pubmind_disease_counts

logger = logging.getLogger("build_site")

# TTN 結構 band 的 genomic 分箱（與 config.TTN_BAND_EDGES / TTN_BAND_LABELS 相同）。
# 報告的 Track B 區塊有 "Structural Band" 卡片時以那個為準；沒跑過推論的報告
# （約 8%）就用座標推算補上，兩者在有交集的報告上一致度 99.97%。
BAND_EDGES = [178525989, 178535982, 178618584, 178775356, 178807423]
BAND_LABELS = ["M-band", "A-band", "I-band", "Z-disk"]
BAND_CHOICES = ["Z-disk", "I-band", "A-band", "M-band"]


def band_from_pos(pos: int) -> str:
    for i, label in enumerate(BAND_LABELS):
        if BAND_EDGES[i] < pos <= BAND_EDGES[i + 1]:
            return label
    return ""


TISSUE_KEYS = ("cardiac", "skeletal", "both", "not_specified")
TISSUE_LABELS = {
    "cardiac": "Cardiac",
    "skeletal": "Skeletal",
    "both": "Both",
    "not_specified": "Not specified",
}
TISSUE_COLORS = {
    "cardiac": "#e74c3c",
    "skeletal": "#f39c12",
    "both": "#9b59b6",
    "not_specified": "#95a5a6",
}

# ── report parsing ────────────────────────────────────────────────────────────

RE_STYLE = re.compile(r"<style>(.*?)</style>", re.S)
RE_IMG = re.compile(r'<img src="data:image/png;base64,([A-Za-z0-9+/=]+)"([^>]*)>')
RE_GENERATED = re.compile(r"Generated: (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})")
RE_CLINVAR = re.compile(r"clinvar/variation/(\d+)")
RE_PREDICTION_BOX = re.compile(
    r'<div class="prediction-box prediction-(\w+)"><div class="prediction-result">([A-Z]+)</div>'
)
RE_SCORE = re.compile(r'<div class="score-value">(-?[\d.]+)</div>')
RE_N_ARTICLES = re.compile(r"Literature Review \((\d+) articles\)")
RE_MENTION = re.compile(r'mention-badge mention-(true|false)')
RE_ARTICLE_LINK = re.compile(r'<a href="(https://pubmed\.ncbi\.nlm\.nih\.gov/\d+/)"[^>]*>(.*?)</a>', re.S)
RE_ARTICLE_META = re.compile(r'<div class="article-meta">(.*?) \((.*?)\) \| PMID: (\d+)</div>', re.S)
RE_INFERENCE_SECTION = re.compile(
    r'<section class="section inference-section">(.*?)</section>', re.S
)
RE_INFERENCE_SOURCE = re.compile(r"Inference source: (.*?)</p>", re.S)
RE_TAG = re.compile(r"<[^>]+>")
# Curated Database Annotations 區塊（HGMD / PubTator 查表）
RE_DB_SECTION = re.compile(r'<section class="section dbann-section', re.S)
RE_DB_CARD = re.compile(
    r'<div class="dbann-db dbann-(hgmd|pubtator)">'
    r'.*?<div class="dbann-disease">(.*?)</div>'
    r'\s*<div class="dbann-db-meta">(.*?)</div>',
    re.S,
)
RE_DB_COUNT = re.compile(r"^(\d+) (?:report|record)")
RE_DB_MULTIPLIER = re.compile(r"\s*×\s*\d+$")


def strip_tags(text: str) -> str:
    return html_mod.unescape(RE_TAG.sub(" ", text)).replace("\xa0", " ").strip()


def collapse_ws(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def field_value(block: str, label: str) -> str:
    """Read one `<label>` / `<value>` pair out of an info-card or article-info item."""
    m = re.search(
        re.escape(label) + r"</div>\s*<div[^>]*>(.*?)</div>\s*</div>", block, re.S
    )
    if not m:
        return ""
    return collapse_ws(strip_tags(m.group(1)))


def normalize_tissue(raw: str) -> str:
    t = (raw or "").strip().lower()
    if "both" in t:
        return "both"
    if "cardiac" in t:
        return "cardiac"
    if "skeletal" in t:
        return "skeletal"
    return "not_specified"


def parse_articles(body: str) -> List[Dict]:
    """Parse the per-article blocks of the Literature Review section."""
    articles: List[Dict] = []
    blocks = body.split('<div class="pubmed-article"')[1:]
    for block in blocks:
        link = RE_ARTICLE_LINK.search(block)
        meta = RE_ARTICLE_META.search(block)
        mention = RE_MENTION.search(block)
        disease = field_value(block, "Disease")
        inferred = "(inferred from article)" in disease
        disease = collapse_ws(disease.replace("(inferred from article)", ""))
        articles.append(
            {
                "pmid": meta.group(3) if meta else "",
                "title": collapse_ws(strip_tags(link.group(2))) if link else "",
                "url": link.group(1) if link else "",
                "journal": collapse_ws(strip_tags(meta.group(1))) if meta else "",
                "year": collapse_ws(strip_tags(meta.group(2))) if meta else "",
                "mention": bool(mention and mention.group(1) == "true"),
                "disease": disease,
                "disease_inferred": inferred,
                "tissue": normalize_tissue(field_value(block, "Tissue Affected")),
            }
        )
    return articles


def parse_inference(source: str) -> Optional[Dict]:
    m = RE_INFERENCE_SECTION.search(source)
    if not m:
        return None
    block = m.group(1)
    src = RE_INFERENCE_SOURCE.search(block)
    return {
        "disease": field_value(block, "Predicted Disease (inference)"),
        "tissue": field_value(block, "Predicted Tissue"),
        "confidence": field_value(block, "Confidence"),
        "differential": field_value(block, "Differential"),
        "band": field_value(block, "Structural Band"),
        "vclass": field_value(block, "Variant Type"),
        "source": collapse_ws(strip_tags(src.group(1))) if src else "",
    }


def parse_databases(source: str) -> Dict:
    """讀出報告裡 HGMD / PubTator 各自登記的疾病。

    `checked` 表示這份報告有跑過資料庫查表；沒有這個區塊的舊報告不能當成
    「不在資料庫裡」，所以在索引上會標成未查詢而不是未命中。
    """
    result: Dict = {
        "checked": bool(RE_DB_SECTION.search(source)),
        "hgmd": [],
        "pubtator": [],
        "n_pubtator": 0,
    }
    if not result["checked"]:
        return result

    for source_name, disease_html, meta in RE_DB_CARD.findall(source):
        diseases = []
        for part in collapse_ws(strip_tags(disease_html)).split("·"):
            name = RE_DB_MULTIPLIER.sub("", part.strip()).strip()
            if name and name != "—":
                diseases.append(name)
        result[source_name] = diseases
        if source_name == "pubtator":
            count = RE_DB_COUNT.match(collapse_ws(strip_tags(meta)))
            result["n_pubtator"] = int(count.group(1)) if count else len(diseases)

    return result


# ── PubMind ──────────────────────────────────────────────────────────────────

def pubmind_for(lookup: Optional["PubMindLookup"], pos: int, ref: str, alt: str) -> Optional[Dict]:
    """索引頁用的 PubMind 摘要；沒有載入資料時回傳 None（表格顯示 not checked）。"""
    if lookup is None or not lookup.loaded:
        return None
    hits = lookup.lookup(pos, ref, alt)
    return {
        "diseases": [d for d, _ in pubmind_disease_counts(hits)],
        "hits": [{k: h[k] for k in ("pvid", "label", "level", "sig", "papers")} for h in hits],
    }


def parse_report(path: Path) -> Optional[Dict]:
    """Extract index metadata from one generated report."""
    source = path.read_text(encoding="utf-8", errors="replace")

    stem = path.stem
    if "_" not in stem:
        logger.warning("skipping %s (cannot split variant/model)", path.name)
        return None
    variant_id, model = stem.rsplit("_", 1)
    parts = variant_id.split("-")
    if len(parts) != 4 or not parts[1].isdigit():
        logger.warning("skipping %s (unexpected variant id %r)", path.name, variant_id)
        return None
    chrom, pos, ref, alt = parts

    evo2: Dict = {}
    box = RE_PREDICTION_BOX.search(source)
    if box:
        scores = RE_SCORE.findall(source)
        evo2 = {
            "prediction": box.group(2).lower(),
            "delta": float(scores[0]) if len(scores) > 0 else None,
            "ref": float(scores[1]) if len(scores) > 1 else None,
            "var": float(scores[2]) if len(scores) > 2 else None,
        }

    articles = parse_articles(source)
    n_declared = RE_N_ARTICLES.search(source)
    tissue_counts = {k: 0 for k in TISSUE_KEYS}
    mention_tissue = {k: 0 for k in TISSUE_KEYS}
    diseases: Counter = Counter()
    for art in articles:
        tissue_counts[art["tissue"]] += 1
        if art["mention"]:
            mention_tissue[art["tissue"]] += 1
        if art["disease"] and art["disease"].lower() not in ("not specified", "n/a", "none"):
            diseases[art["disease"]] += 1

    generated = RE_GENERATED.search(source)
    clinvar = RE_CLINVAR.search(source)
    inference = parse_inference(source)
    databases = parse_databases(source)

    # 檔名後綴 no-llm 指的是「沒有文獻可讓 LLM 擷取臨床資訊」，但這種變異會走
    # Track B，推理內容仍是 LLM 產生的。標成 no-llm 會誤導，改標 inference。
    if model == "no-llm" and inference:
        model = "inference"

    n_mention = sum(1 for a in articles if a["mention"])
    dominant = max(TISSUE_KEYS, key=lambda k: (tissue_counts[k], k != "not_specified"))
    if tissue_counts[dominant] == 0:
        dominant = "not_specified"

    return {
        "variant_id": variant_id,
        "model": model,
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "generated": generated.group(1) if generated else "",
        "clinvar_id": clinvar.group(1) if clinvar else "",
        "evo2": evo2,
        "n_articles": int(n_declared.group(1)) if n_declared else len(articles),
        "n_mention": n_mention,
        "tissue": tissue_counts,
        "mention_tissue": mention_tissue,
        "dominant_tissue": dominant,
        "diseases": [d for d, _ in diseases.most_common(6)],
        "inference": inference,
        "band": ((inference or {}).get("band") or "").strip() or band_from_pos(int(pos)),
        "db": databases,
        "source_file": path.name,
        "source_bytes": len(source.encode("utf-8")),
    }


# ── report rewriting ─────────────────────────────────────────────────────────


def rewrite_report(
    source: str,
    assets_dir: Path,
    css_names: Dict[str, str],
    img_names: Dict[str, str],
    nav_html: str,
) -> str:
    """Point a report at the shared stylesheet + extracted image, and add nav."""
    img_dir = assets_dir / "img"

    def _css(match: re.Match) -> str:
        body = match.group(1)
        digest = hashlib.md5(body.encode("utf-8")).hexdigest()
        name = css_names.get(digest)
        if name is None:
            name = f"report-{digest[:8]}.css"
            css_names[digest] = name
            (assets_dir / name).write_text(body, encoding="utf-8")
        return (
            f'<link rel="stylesheet" href="../assets/{name}">\n'
            f'    <link rel="stylesheet" href="../assets/report-nav.css">'
        )

    def _img(match: re.Match) -> str:
        data = base64.b64decode(match.group(1))
        digest = hashlib.sha1(data).hexdigest()[:16]
        name = img_names.get(digest)
        if name is None:
            name = f"{digest}.png"
            img_names[digest] = name
            (img_dir / name).write_bytes(data)
        attrs = match.group(2)
        if "loading=" not in attrs:
            attrs += ' loading="lazy" decoding="async"'
        return f'<img src="../assets/img/{name}"{attrs}>'

    out = RE_STYLE.sub(_css, source, count=1)
    out = RE_IMG.sub(_img, out)
    out = out.replace('<div class="container">', nav_html + '<div class="container">', 1)
    return out


def make_nav(meta: Dict, prev: Optional[Dict], nxt: Optional[Dict]) -> str:
    def link(item: Optional[Dict], label: str, cls: str) -> str:
        if not item:
            return f'<span class="nav-btn nav-disabled">{label}</span>'
        return (
            f'<a class="nav-btn {cls}" href="{report_filename(item)}" '
            f'title="{html_mod.escape(item["variant_id"])}">{label}</a>'
        )

    return (
        '<nav class="report-nav">'
        '<a class="nav-btn nav-home" href="../index.html">&#8592; All reports</a>'
        f'<span class="nav-title">{html_mod.escape(meta["variant_id"])}'
        f'<span class="nav-model">{html_mod.escape(meta["model"])}</span></span>'
        f'{link(prev, "&#8249; Prev", "nav-prev")}{link(nxt, "Next &#8250;", "nav-next")}'
        "</nav>"
    )


def report_filename(meta: Dict) -> str:
    return f"{meta['variant_id']}__{meta['model']}.html"


# ── static assets ────────────────────────────────────────────────────────────

REPORT_NAV_CSS = """/* Navigation bar injected into each report by build_site.py */
.report-nav {
    max-width: 1200px;
    margin: 0 auto 16px;
    display: flex;
    align-items: center;
    gap: 10px;
    flex-wrap: wrap;
    background: rgba(255, 255, 255, 0.92);
    border-radius: 10px;
    padding: 10px 16px;
    box-shadow: 0 4px 18px rgba(0, 0, 0, 0.15);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
}
.report-nav .nav-title {
    font-weight: 700;
    color: #2c3e50;
    margin-right: auto;
    display: flex;
    align-items: baseline;
    gap: 8px;
}
.report-nav .nav-model {
    font-weight: 400;
    font-size: 0.8em;
    color: #7f8c9a;
    background: #eef1f6;
    border-radius: 4px;
    padding: 2px 7px;
}
.report-nav .nav-btn {
    text-decoration: none;
    color: #3d4d63;
    background: #eef1f6;
    border-radius: 6px;
    padding: 6px 12px;
    font-size: 0.9em;
    transition: background 0.15s, color 0.15s;
}
.report-nav .nav-btn:hover { background: #667eea; color: #fff; }
.report-nav .nav-home { background: #667eea; color: #fff; font-weight: 600; }
.report-nav .nav-home:hover { background: #5a6fd8; }
.report-nav .nav-disabled { opacity: 0.4; cursor: default; }
@media print { .report-nav { display: none; } }
"""

SITE_CSS = """/* TTN variant report index */
:root {
    --ink: #24303f;
    --muted: #6b7a90;
    --line: #e4e9f0;
    --brand: #667eea;
    --brand-2: #764ba2;
    --cardiac: #e74c3c;
    --skeletal: #f39c12;
    --both: #9b59b6;
    --none: #b9c2cf;
}
* { margin: 0; padding: 0; box-sizing: border-box; }
body {
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color: var(--ink);
    background: #f4f6fa;
    line-height: 1.55;
}
a { color: var(--brand); }
header.site {
    background: linear-gradient(135deg, var(--brand) 0%, var(--brand-2) 100%);
    color: #fff;
    padding: 34px 24px 28px;
}
header.site .inner { max-width: 1400px; margin: 0 auto; }
header.site h1 { font-size: 1.9em; letter-spacing: -0.01em; }
header.site p { opacity: 0.88; margin-top: 6px; }
header.site .site-nav { margin-top: 12px; }
header.site .site-nav a {
    color: #fff; text-decoration: none; font-weight: 600; font-size: 0.9em;
    border: 1px solid rgba(255, 255, 255, 0.45); border-radius: 999px; padding: 5px 14px;
}
header.site .site-nav a:hover { background: rgba(255, 255, 255, 0.16); }
main { max-width: 1400px; margin: 0 auto; padding: 22px 24px 60px; }

.stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 14px; margin: -46px 0 22px; }
.stat {
    background: #fff; border-radius: 10px; padding: 16px 18px;
    box-shadow: 0 6px 20px rgba(31, 45, 68, 0.08);
}
.stat .k { font-size: 0.78em; text-transform: uppercase; letter-spacing: 0.06em; color: var(--muted); font-weight: 700; }
.stat .v { font-size: 1.65em; font-weight: 700; margin-top: 2px; }
.stat .sub { font-size: 0.82em; color: var(--muted); }

.toolbar {
    position: sticky; top: 0; z-index: 5;
    background: rgba(244, 246, 250, 0.95);
    backdrop-filter: blur(6px);
    padding: 12px 0 14px;
    border-bottom: 1px solid var(--line);
    display: flex; flex-direction: column; gap: 10px;
}
.toolbar-row { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; }
.db-row-label {
    font-size: 0.78em; text-transform: uppercase; letter-spacing: 0.06em;
    color: var(--muted); font-weight: 700; margin-right: 2px;
}
.toolbar button.db-op {
    width: 34px; height: 34px; border: 1px solid var(--line); border-radius: 999px;
    background: #fff; color: var(--brand); font-size: 1.2em; font-weight: 700;
    font-family: 'Cambria Math', 'STIX Two Math', 'Segoe UI Symbol', 'DejaVu Sans', serif;
    line-height: 1; cursor: pointer;
}
.toolbar button.db-op:hover:not(:disabled) { border-color: var(--brand); background: #f0f4ff; }
.toolbar button.db-op:disabled { color: #c3cad5; cursor: default; }
.db-expr { font-size: 0.85em; color: var(--muted); margin-left: 6px; }
.toolbar button.db-group {
    height: 34px; padding: 0 12px; border: 1px solid var(--line); border-radius: 8px;
    background: #fff; color: var(--brand); font-weight: 600; cursor: pointer; white-space: nowrap;
    font-family: 'Cambria Math', 'STIX Two Math', 'Segoe UI Symbol', 'DejaVu Sans', serif;
}
.toolbar button.db-group:hover:not(:disabled) { border-color: var(--brand); background: #f0f4ff; }
.toolbar button.db-group:disabled { color: #c3cad5; cursor: default; }
.toolbar input[type="search"] {
    flex: 1 1 300px; min-width: 220px;
    padding: 10px 14px; border: 1px solid var(--line); border-radius: 8px;
    font-size: 0.95em; background: #fff; color: var(--ink);
}
.toolbar input[type="search"]:focus { outline: 2px solid var(--brand); border-color: transparent; }
.toolbar select {
    padding: 9px 12px; border: 1px solid var(--line); border-radius: 8px;
    background: #fff; color: var(--ink); font-size: 0.9em;
}
.toolbar .toggle {
    display: inline-flex; align-items: center; gap: 6px;
    background: #fff; border: 1px solid var(--line); border-radius: 8px;
    padding: 9px 12px; font-size: 0.9em; cursor: pointer; user-select: none;
}
.toolbar .toggle input { accent-color: var(--brand); }
.toolbar button.reset {
    padding: 9px 14px; border: 1px solid var(--line); border-radius: 8px;
    background: #fff; color: var(--muted); cursor: pointer; font-size: 0.9em;
}
.toolbar button.reset:hover { color: var(--ink); border-color: var(--muted); }
.result-count { color: var(--muted); font-size: 0.9em; margin: 14px 0 8px; }

.table-wrap { background: #fff; border-radius: 12px; box-shadow: 0 6px 20px rgba(31, 45, 68, 0.07); overflow: hidden; }
table { width: 100%; border-collapse: collapse; font-size: 0.92em; }
thead th {
    background: #f7f9fc; color: var(--muted); text-align: left;
    padding: 12px 14px; font-size: 0.85em; text-transform: uppercase;
    letter-spacing: 0.04em; border-bottom: 1px solid var(--line); white-space: nowrap;
}
thead th.sortable { cursor: pointer; user-select: none; }
thead th.sortable:hover { color: var(--ink); }
thead th .arrow { opacity: 0; margin-left: 4px; }
thead th.sorted .arrow { opacity: 1; }
tbody td { padding: 11px 14px; border-bottom: 1px solid #f1f4f8; vertical-align: middle; }
tbody tr:hover { background: #f9fbff; }
tbody tr:last-child td { border-bottom: none; }
td.variant a { font-weight: 700; text-decoration: none; white-space: nowrap; font-family: ui-monospace, 'SFMono-Regular', Menlo, monospace; }
td.variant a:hover { text-decoration: underline; }
td.variant .sub { display: block; color: var(--muted); font-size: 0.85em; font-family: inherit; }
td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }

.pill { display: inline-block; padding: 2px 9px; border-radius: 999px; font-size: 0.8em; font-weight: 700; white-space: nowrap; }
.pill.pathogenic { background: #fdecea; color: #c0392b; }
.pill.benign { background: #e8f8ef; color: #1e8449; }
.pill.na { background: #eef1f6; color: var(--muted); }
.pill.model { background: #eef1f6; color: #55637a; font-weight: 600; }
.pill.db-hgmd { background: #e8f5e9; color: #1b5e20; }
.pill.db-pubtator { background: #e3f0fc; color: #10508f; }
.pill.db-pubmind { background: #f3e8fd; color: #6a1b9a; }
.pill.db-miss { background: #f2f4f8; color: #a5aeba; font-weight: 600; }
.db-cell { display: flex; flex-wrap: wrap; gap: 4px; align-items: center; max-width: 220px; }
.db-cell .db-disease {
    flex: 1 1 100%; color: var(--muted); font-size: 0.82em;
    white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
}
.pill.conf-high { background: #e8f8ef; color: #1e8449; }
.pill.conf-medium { background: #fff6e0; color: #b9770e; }
.pill.conf-low { background: #eef1f6; color: var(--muted); }

.bar { display: flex; height: 9px; border-radius: 999px; overflow: hidden; background: #f1f4f8; min-width: 96px; }
.bar span { display: block; }
.bar .cardiac { background: var(--cardiac); }
.bar .skeletal { background: var(--skeletal); }
.bar .both { background: var(--both); }
.bar .not_specified { background: var(--none); }
.bar-legend { font-size: 0.8em; color: var(--muted); margin-top: 3px; white-space: nowrap; }
.legend { display: flex; gap: 14px; flex-wrap: wrap; font-size: 0.85em; color: var(--muted); margin: 10px 0 0; }
.legend i { display: inline-block; width: 10px; height: 10px; border-radius: 3px; margin-right: 5px; }

.mention-strong { font-weight: 700; color: #1e8449; }
.mention-zero { color: var(--muted); }
.diseases {
    color: var(--muted); font-size: 0.88em; max-width: 320px;
    /* 這一欄常被擠窄，不夾行會把整列拉得很高 */
    display: -webkit-box; -webkit-line-clamp: 2; -webkit-box-orient: vertical; overflow: hidden;
}
.inference-cell { font-size: 0.88em; }
.inference-cell .dis { display: block; color: #8a5a00; font-weight: 600; }
.empty { padding: 60px 20px; text-align: center; color: var(--muted); }
footer.site { max-width: 1400px; margin: 0 auto; padding: 0 24px 40px; color: var(--muted); font-size: 0.85em; }
@media (max-width: 900px) {
    .stats { margin-top: 16px; }
    .hide-sm { display: none; }
}
/* ── Analyze new variants（TTN-Site 啟動程式提供）─────────────────────────── */
.pill.local { background: #fff4e5; color: #b26a00; }
.analyze {
    background: #fff; border-radius: 10px; padding: 18px 20px; margin-bottom: 18px;
    box-shadow: 0 6px 20px rgba(31, 45, 68, 0.08); border-left: 5px solid #b26a00;
}
.analyze h2 { font-size: 1.1em; display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }
.analyze .hint { color: var(--muted); font-size: 0.88em; margin: 4px 0 12px; }
.analyze summary { cursor: pointer; font-weight: 600; color: var(--brand); margin: 6px 0; }
.an-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 10px 14px; margin: 10px 0; }
.an-grid label { display: flex; flex-direction: column; font-size: 0.8em; color: var(--muted); font-weight: 600; gap: 3px; }
.analyze input[type=text], .analyze input[type=number], .analyze input[type=password], .analyze textarea {
    font: inherit; font-size: 0.95em; padding: 7px 10px; border: 1px solid var(--line); border-radius: 7px; color: var(--ink);
}
.analyze textarea { width: 100%; min-height: 92px; font-family: Consolas, Menlo, monospace; resize: vertical; }
.an-actions { display: flex; gap: 10px; align-items: center; flex-wrap: wrap; margin-top: 8px; }
.analyze button {
    font: inherit; font-weight: 600; padding: 7px 16px; border-radius: 7px; border: 1px solid var(--brand);
    background: var(--brand); color: #fff; cursor: pointer;
}
.analyze button.secondary { background: #fff; color: var(--brand); }
.analyze button:disabled { opacity: 0.5; cursor: default; }
.analyze label.file { font-size: 0.9em; color: var(--brand); font-weight: 600; cursor: pointer; }
.analyze label.file input { display: none; }
.analyze label.an-check { font-size: 0.9em; color: var(--ink); display: flex; align-items: center; gap: 5px; cursor: pointer; }
.an-msg { font-size: 0.88em; }
.an-msg.ok { color: #1e8449; }
.an-msg.err { color: #c0392b; }
.an-box { background: #fff8ec; border: 1px solid #f3d19c; border-radius: 7px; padding: 10px 12px; margin-top: 10px; font-size: 0.9em; }
.an-box code { font-size: 0.95em; word-break: break-all; }
.an-preview { font-size: 0.85em; margin-top: 8px; color: var(--muted); }
.an-preview a { font-weight: 600; }
.an-jobs { width: 100%; border-collapse: collapse; margin-top: 14px; font-size: 0.88em; }
.an-jobs th { text-align: left; color: var(--muted); font-size: 0.85em; text-transform: uppercase; letter-spacing: 0.04em; padding: 6px 8px; border-bottom: 1px solid var(--line); }
.an-jobs td { padding: 7px 8px; border-bottom: 1px solid var(--line); vertical-align: top; }
.an-jobs .vlist { font-family: Consolas, Menlo, monospace; font-size: 0.92em; }
.state { font-weight: 700; }
.state.queued { color: var(--muted); }
.state.waiting { color: #b26a00; }
.state.running { color: #b26a00; }
.state.done { color: #1e8449; }
.state.failed { color: #c0392b; }
.an-offline { color: var(--muted); font-size: 0.85em; margin: -8px 0 16px; }
"""

SITE_JS = r"""(function () {
  'use strict';

  var DATA = window.__TTN_REPORTS__ || [];
  var META = window.__TTN_META__ || {};
  var TISSUES = ['cardiac', 'skeletal', 'both', 'not_specified'];
  var TISSUE_LABELS = { cardiac: 'Cardiac', skeletal: 'Skeletal', both: 'Both', not_specified: 'Not specified' };
  var TISSUE_SHORT = { cardiac: 'Card', skeletal: 'Skel', both: 'Both', not_specified: 'NS' };

  function mergeReports(records) {
    records.forEach(function (r) {
      r.local = true;
      if (r.pos != null) prepare(r);
      for (var i = DATA.length - 1; i >= 0; i--) {
        if (DATA[i].variant_id === r.variant_id) DATA.splice(i, 1);
      }
      DATA.push(r);
    });
  }

  DATA.forEach(prepare);
  // 使用者自己在 Spark 上跑的報告（由 TTN-Site 啟動程式寫入 data/local-reports.js），
  // 同一顆變異以本機版本為準
  mergeReports(window.__TTN_LOCAL_REPORTS__ || []);

  function prepare(r) {
    var bits = [r.variant_id, r.model, r.ref + '>' + r.alt, String(r.pos), r.clinvar_id];
    bits = bits.concat(r.diseases || []);
    if (r.inference) bits.push(r.inference.disease, r.inference.tissue, r.inference.band, r.inference.differential);
    if (r.band) bits.push(r.band);
    if (r.evo2 && r.evo2.prediction) bits.push(r.evo2.prediction);
    if (r.db) bits = bits.concat(r.db.hgmd || [], r.db.pubtator || []);
    if (r.pubmind) {
      bits = bits.concat(r.pubmind.diseases || []);
      (r.pubmind.hits || []).forEach(function (h) { bits.push(h.label, h.pvid); });
    }
    bits.push(TISSUE_LABELS[r.dominant_tissue]);
    if (r.local) bits.push('spark local new');
    r._q = bits.filter(Boolean).join(' ').toLowerCase();
  }

  // 資料庫命中狀態：'in' 有紀錄、'out' 查過但沒有、'' 這份報告沒跑過查表
  function dbState(r, source) {
    // PubMind 是建站時直接比對 JSON，不靠報告內的查表區塊；沒載入資料時為 null
    if (source === 'pubmind') {
      if (!r.pubmind) return '';
      return r.pubmind.hits.length ? 'in' : 'out';
    }
    if (!r.db || !r.db.checked) return '';
    return (r.db[source] || []).length ? 'in' : 'out';
  }

  var DB_SOURCES = ['hgmd', 'pubtator', 'pubmind'];

  function freshState() {
    return {
      q: '', pred: '', tissue: '', band: '', mention: false, inference: false,
      hgmd: '', pubtator: '', pubmind: '', op1: 'and', op2: 'and', group: 'left', sort: 'pos', dir: 1
    };
  }
  var state = freshState();

  var els = {
    q: document.getElementById('q'),
    pred: document.getElementById('f-pred'),
    tissue: document.getElementById('f-tissue'),
    hgmd: document.getElementById('f-hgmd'),
    pubtator: document.getElementById('f-pubtator'),
    pubmind: document.getElementById('f-pubmind'),
    dbops: document.querySelectorAll('button.db-op'),
    dbexpr: document.getElementById('db-expr'),
    dbgroup: document.getElementById('db-group'),
    band: document.getElementById('f-band'),
    mention: document.getElementById('f-mention'),
    inference: document.getElementById('f-inference'),
    reset: document.getElementById('reset'),
    count: document.getElementById('count'),
    body: document.getElementById('rows'),
    stats: document.getElementById('stats'),
    head: document.querySelectorAll('th.sortable')
  };

  function esc(s) {
    return String(s == null ? '' : s).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }

  // op1 是 HGMD 與 PubTator 之間的符號，op2 是 PubTator 與 PubMind 之間的符號。
  // 每個符號歸屬於它右邊的資料庫，用來把該資料庫接到左邊已累積的結果上。
  var DB_OPS = { pubtator: 'op1', pubmind: 'op2' };
  var DB_LABELS = { hgmd: 'HGMD', pubtator: 'PubTator', pubmind: 'PubMind' };

  function activeDbs() {
    return DB_SOURCES.filter(function (k) { return state[k]; });
  }

  // 預設由左到右結合：(A op1 B) op2 C；三個都設定時可切成 A op1 (B op2 C)。
  // 中間的資料庫是 any 時直接跳過，右邊那個資料庫改用它自己左側的符號接上去。
  function groupRight(active) {
    return active.length === 3 && state.group === 'right';
  }

  function combine(op, a, b) { return op === 'or' ? (a || b) : (a && b); }

  function dbMatch(r, active) {
    var m = active.map(function (k) { return dbState(r, k) === state[k]; });
    if (groupRight(active)) return combine(state.op1, m[0], combine(state.op2, m[1], m[2]));
    var acc = m[0];
    for (var i = 1; i < active.length; i++) acc = combine(state[DB_OPS[active[i]]], acc, m[i]);
    return acc;
  }

  function dbExpression(active) {
    var t = active.map(function (k) { return (state[k] === 'out' ? 'not ' : '') + DB_LABELS[k]; });
    var sym = function (k) { return state[DB_OPS[k]] === 'or' ? ' ∪ ' : ' ∩ '; };
    if (active.length === 2) return t[0] + sym(active[1]) + t[1];
    if (groupRight(active)) return t[0] + sym(active[1]) + '(' + t[1] + sym(active[2]) + t[2] + ')';
    return '(' + t[0] + sym(active[1]) + t[1] + ')' + sym(active[2]) + t[2];
  }

  function filtered() {
    var terms = state.q.toLowerCase().split(/\s+/).filter(Boolean);
    var active = activeDbs();
    return DATA.filter(function (r) {
      if (state.pred && (!r.evo2 || r.evo2.prediction !== state.pred)) return false;
      if (state.tissue && r.dominant_tissue !== state.tissue) return false;
      if (state.band && r.band !== state.band) return false;
      if (active.length && !dbMatch(r, active)) return false;
      if (state.mention && !r.n_mention) return false;
      if (state.inference && !r.inference) return false;
      for (var i = 0; i < terms.length; i++) if (r._q.indexOf(terms[i]) === -1) return false;
      return true;
    });
  }

  var SORTERS = {
    pos: function (r) { return r.pos; },
    variant: function (r) { return r.variant_id; },
    model: function (r) { return r.model; },
    pred: function (r) { return (r.evo2 && r.evo2.prediction) || 'zz'; },
    delta: function (r) { return r.evo2 && r.evo2.delta != null ? Math.abs(r.evo2.delta) : -1; },
    articles: function (r) { return r.n_articles; },
    mention: function (r) { return r.n_mention; },
    tissue: function (r) { return TISSUES.indexOf(r.dominant_tissue); },
    // HGMD 權重最高，其次 PubTator、PubMind；三個都沒查過排最後
    db: function (r) {
      var score = 0, checked = false;
      [['hgmd', 4], ['pubtator', 2], ['pubmind', 1]].forEach(function (p) {
        var st = dbState(r, p[0]);
        if (st) checked = true;
        if (st === 'in') score += p[1];
      });
      return checked ? score : -1;
    },
    generated: function (r) { return r.generated || ''; }
  };

  function sorted(rows) {
    var key = SORTERS[state.sort] || SORTERS.pos;
    return rows.slice().sort(function (a, b) {
      var ka = key(a), kb = key(b);
      if (ka < kb) return -1 * state.dir;
      if (ka > kb) return 1 * state.dir;
      return a.pos - b.pos;
    });
  }

  function tissueBar(r) {
    var total = TISSUES.reduce(function (s, k) { return s + (r.tissue[k] || 0); }, 0);
    if (!total) return '<span class="bar-legend">no articles</span>';
    var segs = TISSUES.map(function (k) {
      var n = r.tissue[k] || 0;
      return n ? '<span class="' + k + '" style="width:' + (n / total * 100).toFixed(2) + '%" title="' +
        TISSUE_LABELS[k] + ': ' + n + '"></span>' : '';
    }).join('');
    var legend = TISSUES.filter(function (k) { return r.tissue[k]; })
      .map(function (k) { return TISSUE_SHORT[k] + ' ' + r.tissue[k]; }).join(' · ');
    return '<div class="bar">' + segs + '</div><div class="bar-legend">' + esc(legend) + '</div>';
  }

  function dbCell(r) {
    var dbChecked = r.db && r.db.checked;
    if (!dbChecked && !r.pubmind) return '<span class="mention-zero">not checked</span>';
    var pill = function (source, label, list, extra, tip) {
      if (!list.length) {
        return '<span class="pill db-miss" title="No record for this variant">' + label + '</span>';
      }
      return '<span class="pill db-' + source + '" title="' + esc(tip || list.join(' · ')) + '">' +
        label + (extra || '') + '</span>';
    };
    var hgmd = dbChecked ? r.db.hgmd || [] : [];
    var pubtator = dbChecked ? r.db.pubtator || [] : [];
    var pm = r.pubmind ? r.pubmind.diseases || [] : [];
    var pmHits = r.pubmind ? r.pubmind.hits || [] : [];
    var pmTip = pmHits.map(function (h) {
      return h.label + (h.level === 'protein' ? ' (protein-level match)' : '') +
        (h.sig ? ' · ' + h.sig : '') + ' · ' + h.papers + ' paper' + (h.papers === 1 ? '' : 's') +
        ' · ' + h.pvid;
    }).concat(pm.length ? ['Diseases: ' + pm.join(' · ')] : []).join('\n');

    var html = '';
    if (dbChecked) {
      html += pill('hgmd', 'HGMD', hgmd) +
        pill('pubtator', 'PubTator', pubtator, r.db.n_pubtator > 1 ? ' ' + r.db.n_pubtator : '');
    }
    if (r.pubmind) {
      html += pmHits.length
        ? '<span class="pill db-pubmind" title="' + esc(pmTip) + '">PubMind' +
          (pmHits.length > 1 ? ' ' + pmHits.length : '') + '</span>'
        : pill('pubmind', 'PubMind', []);
    }
    // 一行只放一個疾病（HGMD 優先，因為是人工審編的紀錄），完整清單放 tooltip
    var primary = hgmd[0] || pubtator[0] || pm[0] || '';
    var full = [];
    if (hgmd.length) full.push('HGMD: ' + hgmd.join(' · '));
    if (pubtator.length) full.push('PubTator: ' + pubtator.join(' · '));
    if (pm.length) full.push('PubMind: ' + pm.join(' · '));
    return '<div class="db-cell">' + html +
      (primary ? '<span class="db-disease" title="' + esc(full.join('\n')) + '">' + esc(primary) + '</span>' : '') +
      '</div>';
  }

  function row(r) {
    var evo2 = r.evo2 || {};
    var pred = evo2.prediction
      ? '<span class="pill ' + evo2.prediction + '">' + evo2.prediction.toUpperCase() + '</span>'
      : '<span class="pill na">N/A</span>';
    var delta = evo2.delta != null ? evo2.delta.toFixed(4) : '—';
    var inf = '<span class="mention-zero">—</span>';
    if (r.inference) {
      var conf = (r.inference.confidence || '').toLowerCase();
      inf = '<div class="inference-cell"><span class="dis">' + esc(r.inference.disease || 'Uncertain') + '</span>' +
        '<span class="pill conf-' + (conf || 'low') + '">' + esc(r.inference.confidence || 'Low') + '</span></div>';
    }
    var mention = r.n_mention
      ? '<span class="mention-strong">' + r.n_mention + '</span>'
      : '<span class="mention-zero">0</span>';
    return '<tr>' +
      '<td class="variant"><a href="reports/' + esc(r.file) + '">' + esc(r.variant_id) + '</a>' +
      '<span class="sub">chr' + esc(r.chrom) + ':' + r.pos.toLocaleString() + ' ' + esc(r.ref) + '&rarr;' + esc(r.alt) + '</span></td>' +
      '<td class="hide-sm"><span class="pill model">' + esc(r.model) + '</span>' +
        (r.local ? ' <span class="pill local" title="Analyzed on your DGX Spark">Spark</span>' : '') + '</td>' +
      '<td>' + pred + '</td>' +
      '<td class="num">' + delta + '</td>' +
      '<td class="num">' + r.n_articles + '</td>' +
      '<td class="num">' + mention + '</td>' +
      '<td>' + tissueBar(r) + '</td>' +
      '<td>' + dbCell(r) + '</td>' +
      '<td class="hide-sm">' + inf + '</td>' +
      '<td class="diseases hide-sm" title="' + esc((r.diseases || []).join('\n')) + '">' +
      esc((r.diseases || []).slice(0, 3).join(' · ')) + '</td>' +
      '</tr>';
  }

  function renderStats(rows) {
    var path = 0, benign = 0, articles = 0, mention = 0, withInf = 0, withMention = 0;
    var inHgmd = 0, inPubtator = 0, inPubmind = 0;
    rows.forEach(function (r) {
      if (r.evo2 && r.evo2.prediction === 'pathogenic') path++;
      if (r.evo2 && r.evo2.prediction === 'benign') benign++;
      articles += r.n_articles;
      mention += r.n_mention;
      if (r.inference) withInf++;
      if (r.n_mention) withMention++;
      if (dbState(r, 'hgmd') === 'in') inHgmd++;
      if (dbState(r, 'pubtator') === 'in') inPubtator++;
      if (dbState(r, 'pubmind') === 'in') inPubmind++;
    });
    var cards = [
      { k: 'Reports', v: rows.length, sub: rows.length === DATA.length ? 'all variants' : 'of ' + DATA.length + ' total' },
      { k: 'Evo2 pathogenic', v: path, sub: benign + ' benign' },
      { k: 'Articles', v: articles.toLocaleString(), sub: rows.length ? (articles / rows.length).toFixed(1) + ' per variant' : '—' },
      { k: 'Variant-specific', v: mention.toLocaleString(), sub: withMention + ' variants with ≥1' },
      { k: 'In HGMD', v: inHgmd, sub: inPubtator + ' in PubTator · ' + inPubmind + ' in PubMind' },
      { k: 'Mechanistic inference', v: withInf, sub: 'reports without literature' }
    ];
    els.stats.innerHTML = cards.map(function (c) {
      return '<div class="stat"><div class="k">' + c.k + '</div><div class="v">' + c.v +
        '</div><div class="sub">' + c.sub + '</div></div>';
    }).join('');
  }

  function syncHead() {
    els.head.forEach(function (th) {
      var on = th.dataset.sort === state.sort;
      th.classList.toggle('sorted', on);
      var arrow = th.querySelector('.arrow');
      if (arrow) arrow.textContent = on ? (state.dir === 1 ? '▲' : '▼') : '▲';
    });
  }

  function writeHash() {
    var p = new URLSearchParams();
    if (state.q) p.set('q', state.q);
    if (state.pred) p.set('pred', state.pred);
    if (state.tissue) p.set('tissue', state.tissue);
    if (state.band) p.set('band', state.band);
    if (state.hgmd) p.set('hgmd', state.hgmd);
    if (state.pubtator) p.set('pubtator', state.pubtator);
    if (state.pubmind) p.set('pubmind', state.pubmind);
    if (state.op1 !== 'and') p.set('op1', state.op1);
    if (state.op2 !== 'and') p.set('op2', state.op2);
    if (state.group !== 'left') p.set('group', state.group);
    if (state.mention) p.set('mention', '1');
    if (state.inference) p.set('inference', '1');
    if (state.sort !== 'pos' || state.dir !== 1) p.set('sort', state.sort + ':' + (state.dir === 1 ? 'asc' : 'desc'));
    var s = p.toString();
    history.replaceState(null, '', s ? '#' + s : location.pathname);
  }

  function readHash() {
    var p = new URLSearchParams(location.hash.replace(/^#/, ''));
    state.q = p.get('q') || '';
    state.pred = p.get('pred') || '';
    state.tissue = p.get('tissue') || '';
    state.band = p.get('band') || '';
    state.hgmd = p.get('hgmd') || '';
    state.pubtator = p.get('pubtator') || '';
    state.pubmind = p.get('pubmind') || '';
    // 舊版連結只有單一 dbmode，兩個符號都套用它
    var legacy = p.get('dbmode') === 'or' ? 'or' : 'and';
    state.op1 = p.has('op1') ? (p.get('op1') === 'or' ? 'or' : 'and') : legacy;
    state.op2 = p.has('op2') ? (p.get('op2') === 'or' ? 'or' : 'and') : legacy;
    state.group = p.get('group') === 'right' ? 'right' : 'left';
    state.mention = p.get('mention') === '1';
    state.inference = p.get('inference') === '1';
    var sort = p.get('sort');
    if (sort) {
      var bits = sort.split(':');
      if (SORTERS[bits[0]]) { state.sort = bits[0]; state.dir = bits[1] === 'desc' ? -1 : 1; }
    }
    els.q.value = state.q;
    els.pred.value = state.pred;
    els.tissue.value = state.tissue;
    els.band.value = state.band;
    els.hgmd.value = state.hgmd;
    els.pubtator.value = state.pubtator;
    els.pubmind.value = state.pubmind;
    els.mention.checked = state.mention;
    els.inference.checked = state.inference;
  }

  // 符號只有在右邊的資料庫有設定、且左邊還有其他已設定的資料庫時才會生效
  function syncDbMode() {
    var active = activeDbs();
    els.dbops.forEach(function (btn) {
      var key = btn.getAttribute('data-op');
      var right = key === 'op1' ? 'pubtator' : 'pubmind';
      var leftActive = active.indexOf(right) > 0;
      btn.textContent = state[key] === 'or' ? '∪' : '∩';
      btn.disabled = !leftActive;
      btn.title = !leftActive
        ? 'Set ' + DB_LABELS[right] + ' and a database to its left to use this operator'
        : (state[key] === 'or'
          ? 'Union (∪): click to switch to intersection'
          : 'Intersection (∩): click to switch to union');
    });
    els.dbexpr.textContent = active.length > 1 ? 'Filter: ' + dbExpression(active) : '';
    // 括號位置只在三個資料庫都設定時才有差別
    els.dbgroup.disabled = active.length < 3;
    els.dbgroup.textContent = state.group === 'right' ? 'A ∘ (B ∘ C)' : '(A ∘ B) ∘ C';
    els.dbgroup.title = active.length < 3
      ? 'Set all three database filters to choose where the parentheses go'
      : 'Click to move the parentheses to the ' + (state.group === 'right' ? 'left' : 'right') + ' pair';
  }

  function render() {
    syncDbMode();
    var rows = sorted(filtered());
    renderStats(rows);
    els.count.textContent = rows.length === DATA.length
      ? 'Showing all ' + DATA.length + ' reports'
      : 'Showing ' + rows.length + ' of ' + DATA.length + ' reports';
    els.body.innerHTML = rows.length
      ? rows.map(row).join('')
      : '<tr><td colspan="10" class="empty">No report matches these filters.</td></tr>';
    syncHead();
    writeHash();
  }

  var timer;
  els.q.addEventListener('input', function () {
    clearTimeout(timer);
    timer = setTimeout(function () { state.q = els.q.value.trim(); render(); }, 120);
  });
  ['pred', 'tissue', 'band', 'hgmd', 'pubtator', 'pubmind'].forEach(function (k) {
    els[k].addEventListener('change', function () { state[k] = els[k].value; render(); });
  });
  ['mention', 'inference'].forEach(function (k) {
    els[k].addEventListener('change', function () { state[k] = els[k].checked; render(); });
  });
  els.dbgroup.addEventListener('click', function () {
    state.group = state.group === 'right' ? 'left' : 'right';
    render();
  });
  els.dbops.forEach(function (btn) {
    btn.addEventListener('click', function () {
      var key = btn.getAttribute('data-op');
      state[key] = state[key] === 'or' ? 'and' : 'or';
      render();
    });
  });
  els.reset.addEventListener('click', function () {
    state = freshState();
    els.q.value = ''; els.pred.value = ''; els.tissue.value = '';
    els.band.value = '';
    els.hgmd.value = ''; els.pubtator.value = ''; els.pubmind.value = '';
    els.mention.checked = false; els.inference.checked = false;
    render();
  });
  els.head.forEach(function (th) {
    th.addEventListener('click', function () {
      var key = th.dataset.sort;
      if (state.sort === key) state.dir = -state.dir;
      else { state.sort = key; state.dir = key === 'delta' || key === 'articles' || key === 'mention' ? -1 : 1; }
      render();
    });
  });

  // 給 assets/launcher.js 用：合併新報告後重畫表格
  window.TTN = {
    addReports: function (records) { mergeReports(records); render(); },
    find: function (id) {
      for (var i = 0; i < DATA.length; i++) if (DATA[i].variant_id === id) return DATA[i];
      return null;
    }
  };

  readHash();
  render();
  var built = document.getElementById('built');
  if (built && META.built) built.textContent = META.built;
})();
"""

LAUNCHER_JS = r"""(function () {
  'use strict';
  // 只有透過 TTN-Site 啟動程式開啟（http://127.0.0.1）時才有後端；直接開檔案時顯示提示
  var panel = document.getElementById('analyze');
  var offline = document.getElementById('an-offline');
  if (!panel) return;
  if (location.protocol === 'file:') { offline.hidden = false; return; }

  var VARIANT_RE = /(?:chr)?2-\d+-[ACGT]-[ACGT]/gi;
  var $ = function (id) { return document.getElementById(id); };
  var FIELDS = { host: 'an-host', user: 'an-user', port: 'an-port', key_path: 'an-key', repo: 'an-repo', python: 'an-python' };
  var known = {};      // job_id -> 是否已合併（用來偵測「剛完成」）
  var firstLoad = true;

  function api(method, path, body) {
    return fetch(path, {
      method: method,
      headers: { 'X-TTN-Site': '1', 'Content-Type': 'application/json' },
      body: body ? JSON.stringify(body) : undefined
    }).then(function (r) { return r.json(); });
  }

  function esc(s) {
    return String(s == null ? '' : s).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }

  function msg(id, text, cls) {
    var el = $(id);
    el.textContent = text || '';
    el.className = 'an-msg' + (cls ? ' ' + cls : '');
  }

  // ── 連線設定 ──
  function readForm() {
    var c = {};
    Object.keys(FIELDS).forEach(function (k) { c[k] = $(FIELDS[k]).value.trim(); });
    c.port = parseInt(c.port, 10) || 22;
    return c;
  }

  function handleConnResult(res) {
    $('an-trust').hidden = res.need !== 'trust';
    $('an-pass').hidden = res.need !== 'passphrase';
    if (res.need === 'setup') $('an-setup').hidden = false;
    if (res.ok) $('an-setup').hidden = true;
    if (res.need === 'trust') $('an-fp').textContent = res.fingerprint;
    if (res.ok) {
      msg('an-test-msg', 'Connected. The Spark is ready.', 'ok');
      msg('an-conn', '● Spark connected', 'ok');
      $('an-settings').open = false;
      return true;
    }
    var text = res.error || 'Connection failed';
    if (res.check && res.check.issues && res.check.issues.length) text = res.check.issues.join('; ');
    msg('an-test-msg', text, 'err');
    msg('an-conn', '● Spark not connected', 'err');
    $('an-settings').open = true;
    return false;
  }

  $('an-save').addEventListener('click', function () {
    msg('an-test-msg', 'Connecting…');
    api('POST', '/api/config', readForm())
      .then(function () { return api('POST', '/api/test'); })
      .then(handleConnResult)
      .catch(function (e) { msg('an-test-msg', String(e), 'err'); });
  });
  $('an-trust-btn').addEventListener('click', function () {
    msg('an-test-msg', 'Connecting…');
    api('POST', '/api/trust', { fingerprint: $('an-fp').textContent }).then(handleConnResult);
  });
  $('an-setup-show').addEventListener('click', function () {
    $('an-setup').hidden = !$('an-setup').hidden;
  });
  $('an-setup-btn').addEventListener('click', function () {
    var p = $('an-password');
    msg('an-test-msg', 'Setting up the SSH key…');
    api('POST', '/api/config', readForm())
      .then(function () { return api('POST', '/api/setup-key', { password: p.value }); })
      .then(function (res) {
        p.value = '';
        if (handleConnResult(res)) {
          msg('an-test-msg', 'SSH key ' + (res.key_created ? 'created and ' : '') +
            'installed on the Spark. Connected.', 'ok');
        }
      })
      .catch(function (e) { msg('an-test-msg', String(e), 'err'); });
  });
  $('an-pass-btn').addEventListener('click', function () {
    msg('an-test-msg', 'Connecting…');
    var p = $('an-passphrase');
    api('POST', '/api/passphrase', { passphrase: p.value }).then(function (res) { p.value = ''; handleConnResult(res); });
  });

  // ── 輸入變異 ──
  function parseVariants(text) {
    var seen = {}, out = [];
    (text.match(VARIANT_RE) || []).forEach(function (v) {
      v = v.replace(/^chr/i, '').toUpperCase();
      if (!seen[v]) { seen[v] = true; out.push(v); }
    });
    return out;
  }

  // 已在 site 裡的變異預設只給連結；勾選「Re-run」時也一併送出重新分析
  function preview() {
    var vs = parseVariants($('an-variants').value);
    var rerun = $('an-rerun').checked;
    var have = [], todo = [], fresh = 0;
    vs.forEach(function (v) {
      var r = window.TTN && window.TTN.find(v);
      if (r) have.push('<a href="reports/' + esc(r.file) + '">' + esc(v) + '</a>');
      if (!r) fresh++;
      if (!r || rerun) todo.push(v);
    });
    var parts = [];
    if (fresh) parts.push(fresh + ' new variant(s) will be analyzed');
    if (have.length) {
      parts.push((rerun ? have.length + ' already in this site will be re-run: '
                        : 'already in this site (tick “Re-run” to analyze again): ') + have.join(', '));
    }
    $('an-preview').innerHTML = parts.join(' · ');
    $('an-submit').disabled = !todo.length;
    return todo;
  }
  $('an-variants').addEventListener('input', preview);
  $('an-rerun').addEventListener('change', preview);

  $('an-file').addEventListener('change', function (e) {
    var file = e.target.files[0];
    if (!file) return;
    var reader = new FileReader();
    reader.onload = function () {
      var found = parseVariants(String(reader.result));
      var box = $('an-variants');
      box.value = (box.value.trim() ? box.value.trim() + '\n' : '') + found.join('\n');
      preview();
      msg('an-submit-msg', found.length + ' variant(s) read from ' + file.name, found.length ? 'ok' : 'err');
    };
    reader.readAsText(file);
    e.target.value = '';
  });

  $('an-submit').addEventListener('click', function () {
    var todo = preview();
    if (!todo.length) return;
    $('an-submit').disabled = true;
    msg('an-submit-msg', 'Sending to the Spark…');
    var rerun = $('an-rerun').checked;
    api('POST', '/api/submit', { variants: todo, rerun: rerun }).then(function (res) {
      $('an-submit').disabled = false;
      // 成功時被 Spark 擋下的變異（例如 ref 不符）放在 job.rejected
      var rejected = ((res.ok && res.job ? res.job.rejected : res.rejected) || [])
        .map(function (r) { return r.variant + ' (' + r.reason + ')'; });
      if (res.ok) {
        msg('an-submit-msg', 'Queued ' + res.job.variants.length + ' variant(s)' + (rerun ? ' (re-run)' : '') + '.' +
          (rejected.length ? ' Skipped: ' + rejected.join('; ') : ''), rejected.length ? 'err' : 'ok');
        $('an-variants').value = '';
        preview();
        if (window.Notification && Notification.permission === 'default') Notification.requestPermission();
        refreshJobs();
      } else {
        handleConnResult(res);
        msg('an-submit-msg', (res.error || 'Submission failed') + (rejected.length ? ': ' + rejected.join('; ') : ''), 'err');
      }
    }).catch(function (e) { $('an-submit').disabled = false; msg('an-submit-msg', String(e), 'err'); });
  });

  // ── 任務清單 ──
  function stateText(j) {
    if (j.state === 'done' && !j.merged) return j.merge_error ? 'downloading failed' : 'downloading…';
    if (j.state === 'running') return 'running (' + (j.completed || []).length + '/' + j.variants.length + ')';
    if (j.state === 'waiting') return 'waiting for memory';
    return j.state;
  }

  function resultCell(j) {
    if (j.state === 'failed') return '<span class="an-msg err">' + esc(j.error || 'failed') + '</span>';
    if (j.state === 'waiting') return '<span class="an-msg">' + esc(j.note || 'Waiting for the Spark to free memory') + '</span>';
    if (j.merge_error) return '<span class="an-msg err">' + esc(j.merge_error) + '</span>';
    if (!j.merged) return '—';
    var links = (j.reports || []).map(function (f) {
      return '<a href="reports/' + esc(f) + '">' + esc(f.split('__')[0]) + '</a>';
    });
    var failed = (j.failed_variants || []).length
      ? ' <span class="an-msg err">no report: ' + esc(j.failed_variants.join(', ')) + '</span>' : '';
    return links.join(', ') + failed;
  }

  function refreshJobs() {
    return api('GET', '/api/jobs').then(function (res) {
      var jobs = (res.jobs || []).slice().reverse();
      $('an-jobs').hidden = !jobs.length;
      $('an-jobs').querySelector('tbody').innerHTML = jobs.map(function (j) {
        return '<tr><td>' + esc((j.created || '').replace('T', ' ')) + '<br><span class="hint">' + esc(j.host) + '</span></td>' +
          '<td class="vlist">' + j.variants.map(esc).join('<br>') +
            (j.rerun ? '<br><span class="pill local">re-run</span>' : '') + '</td>' +
          '<td><span class="state ' + esc(j.state) + '">' + esc(stateText(j)) + '</span></td>' +
          '<td>' + resultCell(j) + '</td></tr>';
      }).join('');
      if (res.poll_error) msg('an-conn', '● ' + res.poll_error, 'err');

      var newlyMerged = jobs.filter(function (j) { return j.merged && known[j.job_id] === false; });
      jobs.forEach(function (j) { known[j.job_id] = !!j.merged; });
      if (newlyMerged.length || firstLoad) reloadReports(newlyMerged);
      firstLoad = false;
    });
  }

  function reloadReports(finished) {
    api('GET', '/api/local-reports').then(function (res) {
      if (res.ok && window.TTN) window.TTN.addReports(res.reports || []);
      finished.forEach(function (j) {
        var text = (j.reports || []).length + ' variant(s) added to the site';
        if (window.Notification && Notification.permission === 'granted') {
          new Notification('TTN analysis finished', { body: text });
        }
      });
      preview();
    });
  }

  // ── 啟動 ──
  api('GET', '/api/ping').then(function (res) {
    if (!res.ok) return;
    panel.hidden = false;
    return api('GET', '/api/config').then(function (res) {
      var c = res.config || {};
      Object.keys(FIELDS).forEach(function (k) { if (c[k] != null) $(FIELDS[k]).value = c[k]; });
      if (!c.host || !c.user) {
        $('an-settings').open = true;
        msg('an-conn', '● set up the Spark connection first', 'err');
      } else {
        api('POST', '/api/test').then(handleConnResult);
      }
      preview();
      refreshJobs();
      setInterval(refreshJobs, 5000);
    });
  }).catch(function () { offline.hidden = false; });
})();
"""

INDEX_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>__TITLE__</title>
<link rel="stylesheet" href="assets/site.css">
</head>
<body>
<header class="site">
    <div class="inner">
        <h1>__TITLE__</h1>
        <p>__SUBTITLE__</p>
        __NAV_EXTRA__
    </div>
</header>
<main>
    <section class="stats" id="stats"></section>

    <section class="analyze" id="analyze" hidden>
        <h2>Analyze variants on your DGX Spark <span class="an-msg" id="an-conn"></span></h2>
        <p class="hint">Variants are sent to your DGX Spark over SSH, run through the full pipeline
            (Evo2, literature, database lookup, inference), and added to this site when they finish.
            Keep the TTN-Site window open until then.</p>
        <details id="an-settings">
            <summary>DGX Spark connection</summary>
            <div class="an-grid">
                <label>Spark host / IP<input type="text" id="an-host" placeholder="192.168.1.20" autocomplete="off"></label>
                <label>User name on the Spark<input type="text" id="an-user" autocomplete="off"></label>
                <label>SSH port<input type="number" id="an-port" value="22" min="1" max="65535"></label>
                <label>Private key (blank = ~/.ssh/id_ed25519, id_ecdsa, id_rsa)<input type="text" id="an-key" autocomplete="off"></label>
                <label>ttn-agent folder on the Spark<input type="text" id="an-repo" autocomplete="off"></label>
                <label>Python on the Spark (gm4 environment)<input type="text" id="an-python" autocomplete="off"></label>
            </div>
            <div class="an-actions">
                <button type="button" id="an-save">Save &amp; test connection</button>
                <button type="button" class="secondary" id="an-setup-show">Set up SSH key…</button>
                <span class="an-msg" id="an-test-msg"></span>
            </div>
            <div class="an-box" id="an-setup" hidden>
                <b>One-time SSH key setup.</b> Enter the password you use to log in to the Spark.
                TTN-Site creates an SSH key on this computer if there is none yet
                (<code>%USERPROFILE%&#92;.ssh&#92;id_ed25519</code>) and installs it on the Spark, so you
                will not need the password again. The password is used once and never saved.
                <div class="an-actions">
                    <input type="password" id="an-password" autocomplete="off" placeholder="Spark password">
                    <button type="button" id="an-setup-btn">Set up key</button>
                </div>
            </div>
            <div class="an-box" id="an-trust" hidden>
                First connection to this Spark. Check that its SSH key fingerprint is
                <code id="an-fp"></code>
                (on the Spark: <code>ssh-keygen -lf /etc/ssh/ssh_host_ed25519_key.pub</code>), then
                <button type="button" id="an-trust-btn">Trust this Spark</button>
            </div>
            <div class="an-box" id="an-pass" hidden>
                Your SSH key has a passphrase. It is kept in memory only until TTN-Site is closed.
                <div class="an-actions">
                    <input type="password" id="an-passphrase" autocomplete="off">
                    <button type="button" id="an-pass-btn">Unlock key</button>
                </div>
            </div>
        </details>
        <textarea id="an-variants" placeholder="One variant per line, e.g.&#10;2-178527098-C-T&#10;2-178612477-T-A"></textarea>
        <div class="an-preview" id="an-preview"></div>
        <div class="an-actions">
            <label class="file">Load CSV / text file…<input type="file" id="an-file" accept=".csv,.txt,.tsv"></label>
            <label class="an-check" title="Variants that already have a report are normally just linked. Tick this to run the full pipeline on them again; the new report replaces the old one in this site.">
                <input type="checkbox" id="an-rerun"> Re-run variants already in this site</label>
            <button type="button" id="an-submit">Run on Spark</button>
            <span class="an-msg" id="an-submit-msg"></span>
        </div>
        <table class="an-jobs" id="an-jobs" hidden>
            <thead><tr><th>Submitted</th><th>Variants</th><th>Status</th><th>Result</th></tr></thead>
            <tbody></tbody>
        </table>
    </section>
    <p class="an-offline" id="an-offline" hidden>To analyze variants that are not in this site,
        open the site with <b>TTN-Site.exe</b> (in this folder) instead of opening index.html directly.</p>


    <div class="toolbar">
        <div class="toolbar-row">
            <input type="search" id="q" placeholder="Search variant, position, disease, band, ClinVar ID…" autocomplete="off">
            <select id="f-pred">
                <option value="">Evo2: any</option>
                <option value="pathogenic">Pathogenic</option>
                <option value="benign">Benign</option>
            </select>
            <select id="f-tissue">
                <option value="">Tissue: any</option>
                <option value="cardiac">Cardiac</option>
                <option value="skeletal">Skeletal</option>
                <option value="both">Both</option>
                <option value="not_specified">Not specified</option>
            </select>
            <select id="f-band" title="Structural band of titin that the variant falls in (Protein Domain Localization)">
                <option value="">Band: any</option>
                <option value="Z-disk">Z-disk</option>
                <option value="I-band">I-band</option>
                <option value="A-band">A-band</option>
                <option value="M-band">M-band</option>
            </select>
            <label class="toggle"><input type="checkbox" id="f-mention"> Mentions variant</label>
            <label class="toggle"><input type="checkbox" id="f-inference"> Has inference</label>
            <button class="reset" id="reset">Reset</button>
        </div>
        <div class="toolbar-row db-row">
            <span class="db-row-label">Databases</span>
            <select id="f-hgmd" title="Whether this exact variant is catalogued in HGMD">
                <option value="">HGMD: any</option>
                <option value="in">In HGMD</option>
                <option value="out">Not in HGMD</option>
            </select>
            <button type="button" class="db-op" data-op="op1">∩</button>
            <select id="f-pubtator" title="Whether this exact variant is mentioned in the PubTator-mined corpus">
                <option value="">PubTator: any</option>
                <option value="in">In PubTator</option>
                <option value="out">Not in PubTator</option>
            </select>
            <button type="button" class="db-op" data-op="op2">∩</button>
            <select id="f-pubmind" title="Whether this variant is recorded in PubMind (DNA-level, or protein-level for records given only as an amino-acid change)">
                <option value="">PubMind: any</option>
                <option value="in">In PubMind</option>
                <option value="out">Not in PubMind</option>
            </select>
            <button type="button" class="db-group" id="db-group">(A ∘ B) ∘ C</button>
            <span class="db-expr" id="db-expr"></span>
        </div>
    </div>

    <p class="result-count" id="count"></p>

    <div class="table-wrap">
        <table>
            <thead>
                <tr>
                    <th class="sortable" data-sort="pos">Variant<span class="arrow"></span></th>
                    <th class="sortable hide-sm" data-sort="model">Model<span class="arrow"></span></th>
                    <th class="sortable" data-sort="pred">Evo2<span class="arrow"></span></th>
                    <th class="sortable num" data-sort="delta" title="Evo2 delta score; sorted by absolute value">&Delta; score<span class="arrow"></span></th>
                    <th class="sortable num" data-sort="articles">Articles<span class="arrow"></span></th>
                    <th class="sortable num" data-sort="mention">Mention<span class="arrow"></span></th>
                    <th class="sortable" data-sort="tissue">Tissue split<span class="arrow"></span></th>
                    <th class="sortable" data-sort="db" title="Curated database record for this exact variant">Databases<span class="arrow"></span></th>
                    <th class="hide-sm">Inference</th>
                    <th class="hide-sm">Top diseases</th>
                </tr>
            </thead>
            <tbody id="rows"></tbody>
        </table>
    </div>

    <div class="legend">
        <span><i style="background:#e74c3c"></i>Cardiac</span>
        <span><i style="background:#f39c12"></i>Skeletal</span>
        <span><i style="background:#9b59b6"></i>Both</span>
        <span><i style="background:#b9c2cf"></i>Not specified</span>
    </div>
</main>
<footer class="site">
    <p>Generated by TTN Variant AI Agent · site built <span id="built"></span> · __N_REPORTS__ reports</p>
</footer>
<script id="report-data" type="application/json">__DATA__</script>
<script>
window.__TTN_REPORTS__ = JSON.parse(document.getElementById('report-data').textContent);
window.__TTN_META__ = __META__;
</script>
<script src="data/local-reports.js"></script>
<script src="assets/site.js"></script>
<script src="assets/launcher.js"></script>
</body>
</html>
"""


# ── build ────────────────────────────────────────────────────────────────────


def has_inference(path: Path) -> bool:
    """報告裡是否有 Track B 推理區塊。"""
    try:
        return 'class="section inference-section"' in path.read_text(
            encoding="utf-8", errors="replace"
        )
    except OSError:
        return False


# 同一顆變異可能有多份報告（不同模型 tag）。預設只收一份，依此優先序挑：
# NVFP4 量化的 Gemma 4（正式部署用的）→ 未量化的 Gemma 4 → no-llm（Track B 推理）。
# 清單以外的 tag 排在中間，no-llm 一律最後。
MODEL_PRIORITY = ["Gemma-4-31B-IT-NVFP4", "gemma-4-31B-it"]


def model_rank(tag: str) -> int:
    if tag in MODEL_PRIORITY:
        return MODEL_PRIORITY.index(tag)
    return len(MODEL_PRIORITY) + (1 if tag == "no-llm" else 0)


def dedupe_by_variant(files: List[Path]) -> List[Path]:
    """一顆變異只留一份報告，依 MODEL_PRIORITY 決定留哪份。"""
    best: Dict[str, Path] = {}
    for f in files:
        if "_" not in f.stem:
            continue
        variant, tag = f.stem.rsplit("_", 1)
        cur = best.get(variant)
        if cur is None or model_rank(tag) < model_rank(cur.stem.rsplit("_", 1)[1]):
            best[variant] = f
    return sorted(best.values())


def select_reports(outputs: Path, include: str) -> List[Path]:
    files = sorted(outputs.glob("*.html"))
    picked = []
    for f in files:
        if "smoketest" in f.stem:
            continue
        is_no_llm = f.stem.endswith("_no-llm")
        # llm 模式的用意是「只發佈有模型分析內容的報告」。無文獻的變異檔名雖然是
        # no-llm，只要補上了 Track B 推理就同樣有內容，不該被排除。
        if include == "llm" and is_no_llm and not has_inference(f):
            continue
        if include == "no-llm" and not is_no_llm:
            continue
        picked.append(f)
    return picked


INDEX_FIELDS = (
    "variant_id", "model", "chrom", "pos", "ref", "alt", "generated",
    "clinvar_id", "evo2", "n_articles", "n_mention", "tissue",
    "mention_tissue", "dominant_tissue", "diseases", "inference", "db", "pubmind", "band",
    "file",
)


def index_record(meta: Dict) -> Dict:
    return {k: meta[k] for k in INDEX_FIELDS}


def build_patch(report_paths: List[Path], patch_dir: Path,
                pubmind_path: Optional[Path] = None) -> List[Dict]:
    """把幾份報告轉成能直接合併進「已發佈 site」的檔案，給 spark_job.py 用。

    產出的目錄結構與 site/ 相同（reports/、assets/、assets/img/），外加
    records.json（索引列）。使用者端的啟動程式只要把檔案複製過去、把索引列併進
    data/local-reports.js，不需要 Python，也不需要重建整個 site。
    前後頁連結只指回索引，因為新報告在使用者那份 site 裡的鄰居無從得知。
    """
    reports_dir = patch_dir / "reports"
    assets_dir = patch_dir / "assets"
    (assets_dir / "img").mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    pubmind = PubMindLookup(pubmind_path) if pubmind_path else None

    records: List[Dict] = []
    css_names: Dict[str, str] = {}
    img_names: Dict[str, str] = {}
    for path in report_paths:
        meta = parse_report(path)
        if meta is None:
            continue
        meta["pubmind"] = pubmind_for(pubmind, meta["pos"], meta["ref"], meta["alt"])
        meta["file"] = report_filename(meta)
        source = path.read_text(encoding="utf-8", errors="replace")
        out = rewrite_report(source, assets_dir, css_names, img_names, make_nav(meta, None, None))
        (reports_dir / meta["file"]).write_text(out, encoding="utf-8")
        records.append(index_record(meta))
    (patch_dir / "records.json").write_text(
        json.dumps(records, ensure_ascii=False, separators=(",", ":")), encoding="utf-8"
    )
    return records


def build(outputs: Path, site: Path, include: str, inline_images: bool,
          dedupe: bool = True, pubmind_path: Optional[Path] = None) -> Dict:
    from datetime import datetime

    reports_dir = site / "reports"
    assets_dir = site / "assets"
    img_dir = assets_dir / "img"
    data_dir = site / "data"
    for d in (reports_dir, assets_dir, img_dir, data_dir):
        d.mkdir(parents=True, exist_ok=True)

    files = select_reports(outputs, include)
    if dedupe:
        before = len(files)
        files = dedupe_by_variant(files)
        logger.info("deduplicated %d reports to %d variants", before, len(files))
    if not files:
        raise SystemExit(f"No reports matched --include {include} in {outputs}")
    logger.info("found %d report(s) in %s", len(files), outputs)

    pubmind = PubMindLookup(pubmind_path) if pubmind_path else None

    parsed: List[Dict] = []
    sources: Dict[str, str] = {}
    for path in files:
        meta = parse_report(path)
        if meta is None:
            continue
        meta["pubmind"] = pubmind_for(pubmind, meta["pos"], meta["ref"], meta["alt"])
        meta["file"] = report_filename(meta)
        parsed.append(meta)
        sources[meta["file"]] = path.read_text(encoding="utf-8", errors="replace")

    parsed.sort(key=lambda m: (m["pos"], m["variant_id"], m["model"]))

    css_names: Dict[str, str] = {}
    img_names: Dict[str, str] = {}
    (assets_dir / "report-nav.css").write_text(REPORT_NAV_CSS, encoding="utf-8")

    written_bytes = 0
    for i, meta in enumerate(parsed):
        nav = make_nav(meta, parsed[i - 1] if i else None, parsed[i + 1] if i + 1 < len(parsed) else None)
        source = sources[meta["file"]]
        if inline_images:
            out = source.replace('<div class="container">', nav + '<div class="container">', 1)
            out = out.replace(
                "<style>",
                '<link rel="stylesheet" href="../assets/report-nav.css">\n    <style>',
                1,
            )
        else:
            out = rewrite_report(source, assets_dir, css_names, img_names, nav)
        dest = reports_dir / meta["file"]
        dest.write_text(out, encoding="utf-8")
        written_bytes += len(out.encode("utf-8"))
        if (i + 1) % 100 == 0:
            logger.info("  rewrote %d/%d reports", i + 1, len(parsed))

    (assets_dir / "site.css").write_text(SITE_CSS, encoding="utf-8")
    (assets_dir / "site.js").write_text(SITE_JS, encoding="utf-8")
    (assets_dir / "launcher.js").write_text(LAUNCHER_JS, encoding="utf-8")
    # 使用者端啟動程式（launcher/ 用 Go 編譯）；有編譯好的就一起放進 site
    for exe in sorted((Path(__file__).resolve().parent / "launcher" / "dist").glob("TTN-Site*.exe")):
        shutil.copy2(exe, site / exe.name)

    index_records = [index_record(meta) for meta in parsed]
    payload = json.dumps(index_records, ensure_ascii=False, separators=(",", ":"))
    (data_dir / "reports.json").write_text(payload, encoding="utf-8")

    built = datetime.now().strftime("%Y-%m-%d %H:%M")
    site_meta = {"built": built, "n_reports": len(parsed), "include": include}
    # 驗證摘要頁是另外產生的；重建索引時若不主動補連結就會把入口弄丟。
    nav_extra = (
        '<p class="site-nav"><a href="validation.html">'
        "Batch inference validation &rarr;</a></p>"
        if (site / "validation.html").exists()
        else ""
    )
    index_html = (
        INDEX_HTML.replace("__TITLE__", "TTN Variant AI Reports")
        .replace("__NAV_EXTRA__", nav_extra)
        .replace(
            "__SUBTITLE__",
            f"{len(parsed)} variant reports · Evo2 pathogenicity, literature evidence "
            f"and mechanistic inference · built {built}",
        )
        .replace("__N_REPORTS__", str(len(parsed)))
        .replace("__DATA__", payload.replace("</", "<\\/"))
        .replace("__META__", json.dumps(site_meta, ensure_ascii=False))
    )
    (site / "index.html").write_text(index_html, encoding="utf-8")

    img_bytes = sum(p.stat().st_size for p in img_dir.glob("*.png"))
    return {
        "n_reports": len(parsed),
        "n_images": len(img_names),
        "html_mb": written_bytes / 1024 / 1024,
        "img_mb": img_bytes / 1024 / 1024,
        "source_mb": sum(m["source_bytes"] for m in parsed) / 1024 / 1024,
        "index_kb": (site / "index.html").stat().st_size / 1024,
        "site": site,
    }


def main(argv: Optional[List[str]] = None) -> int:
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", type=Path, default=root / "outputs", help="directory holding the generated reports")
    ap.add_argument("--site", type=Path, default=root / "site", help="output directory for the static site")
    ap.add_argument(
        "--include",
        choices=("llm", "all", "no-llm"),
        default="llm",
        help="which reports to publish (default: llm — reports with LLM-derived content, "
             "i.e. everything except *_no-llm.html that also lacks a Track B inference section)",
    )
    ap.add_argument(
        "--inline-images",
        action="store_true",
        help="keep the base64 schematics inline instead of extracting them to assets/img",
    )
    ap.add_argument(
        "--keep-duplicates",
        action="store_true",
        help=("保留同一顆變異的多份報告。預設只留一份，優先序 "
              "Gemma-4-31B-IT-NVFP4 > gemma-4-31B-it > no-llm"),
    )
    ap.add_argument(
        "--pubmind",
        type=Path,
        default=root / "data" / "pubmind_ttn.json",
        help="PubMind TTN export to match report variants against",
    )
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(message)s",
    )

    stats = build(args.outputs, args.site, args.include, args.inline_images,
                  dedupe=not args.keep_duplicates, pubmind_path=args.pubmind)
    logger.info(
        "\nsite: %s\n  %d reports (%.1f MB html, was %.1f MB inline)\n"
        "  %d schematics (%.1f MB, lazy-loaded)\n  index.html %.0f KB\n"
        "\nopen %s directly, or serve it:\n  python -m http.server -d %s 8000",
        stats["site"],
        stats["n_reports"],
        stats["html_mb"],
        stats["source_mb"],
        stats["n_images"],
        stats["img_mb"],
        stats["index_kb"],
        stats["site"] / "index.html",
        stats["site"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
