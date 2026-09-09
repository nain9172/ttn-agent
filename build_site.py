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
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger("build_site")

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
RE_DB_COUNT = re.compile(r"^(\d+) report")
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
    display: flex; flex-wrap: wrap; gap: 10px; align-items: center;
}
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
"""

SITE_JS = r"""(function () {
  'use strict';

  var DATA = window.__TTN_REPORTS__ || [];
  var META = window.__TTN_META__ || {};
  var TISSUES = ['cardiac', 'skeletal', 'both', 'not_specified'];
  var TISSUE_LABELS = { cardiac: 'Cardiac', skeletal: 'Skeletal', both: 'Both', not_specified: 'Not specified' };
  var TISSUE_SHORT = { cardiac: 'Card', skeletal: 'Skel', both: 'Both', not_specified: 'NS' };

  DATA.forEach(function (r) {
    var bits = [r.variant_id, r.model, r.ref + '>' + r.alt, String(r.pos), r.clinvar_id];
    bits = bits.concat(r.diseases || []);
    if (r.inference) bits.push(r.inference.disease, r.inference.tissue, r.inference.band, r.inference.differential);
    if (r.evo2 && r.evo2.prediction) bits.push(r.evo2.prediction);
    if (r.db) bits = bits.concat(r.db.hgmd || [], r.db.pubtator || []);
    bits.push(TISSUE_LABELS[r.dominant_tissue]);
    r._q = bits.filter(Boolean).join(' ').toLowerCase();
  });

  // 資料庫命中狀態：'in' 有紀錄、'out' 查過但沒有、'' 這份報告沒跑過查表
  function dbState(r, source) {
    if (!r.db || !r.db.checked) return '';
    return (r.db[source] || []).length ? 'in' : 'out';
  }

  var state = {
    q: '', model: '', pred: '', tissue: '', mention: false, inference: false,
    hgmd: '', pubtator: '', sort: 'pos', dir: 1
  };

  var els = {
    q: document.getElementById('q'),
    model: document.getElementById('f-model'),
    pred: document.getElementById('f-pred'),
    tissue: document.getElementById('f-tissue'),
    hgmd: document.getElementById('f-hgmd'),
    pubtator: document.getElementById('f-pubtator'),
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

  function fillModels() {
    var models = {};
    DATA.forEach(function (r) { models[r.model] = (models[r.model] || 0) + 1; });
    Object.keys(models).sort().forEach(function (m) {
      var opt = document.createElement('option');
      opt.value = m;
      opt.textContent = m + ' (' + models[m] + ')';
      els.model.appendChild(opt);
    });
  }

  function filtered() {
    var terms = state.q.toLowerCase().split(/\s+/).filter(Boolean);
    return DATA.filter(function (r) {
      if (state.model && r.model !== state.model) return false;
      if (state.pred && (!r.evo2 || r.evo2.prediction !== state.pred)) return false;
      if (state.tissue && r.dominant_tissue !== state.tissue) return false;
      if (state.hgmd && dbState(r, 'hgmd') !== state.hgmd) return false;
      if (state.pubtator && dbState(r, 'pubtator') !== state.pubtator) return false;
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
    // 兩邊都有 > 只有 HGMD > 只有 PubTator > 都沒有 > 沒查過
    db: function (r) {
      if (!r.db || !r.db.checked) return -1;
      return (r.db.hgmd || []).length * 2 + ((r.db.pubtator || []).length ? 1 : 0);
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
    if (!r.db || !r.db.checked) return '<span class="mention-zero">not checked</span>';
    var pill = function (source, label) {
      var list = r.db[source] || [];
      if (!list.length) {
        return '<span class="pill db-miss" title="No record for this variant">' + label + '</span>';
      }
      var extra = source === 'pubtator' && r.db.n_pubtator > 1 ? ' ' + r.db.n_pubtator : '';
      return '<span class="pill db-' + source + '" title="' + esc(list.join(' · ')) + '">' +
        label + extra + '</span>';
    };
    // 一行只放一個疾病（HGMD 優先，因為是人工審編的紀錄），完整清單放 tooltip
    var primary = (r.db.hgmd || [])[0] || (r.db.pubtator || [])[0] || '';
    var full = [];
    if ((r.db.hgmd || []).length) full.push('HGMD: ' + r.db.hgmd.join(' · '));
    if ((r.db.pubtator || []).length) full.push('PubTator: ' + r.db.pubtator.join(' · '));
    return '<div class="db-cell">' + pill('hgmd', 'HGMD') + pill('pubtator', 'PubTator') +
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
      '<td class="hide-sm"><span class="pill model">' + esc(r.model) + '</span></td>' +
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
    var inHgmd = 0, inPubtator = 0;
    rows.forEach(function (r) {
      if (r.evo2 && r.evo2.prediction === 'pathogenic') path++;
      if (r.evo2 && r.evo2.prediction === 'benign') benign++;
      articles += r.n_articles;
      mention += r.n_mention;
      if (r.inference) withInf++;
      if (r.n_mention) withMention++;
      if (dbState(r, 'hgmd') === 'in') inHgmd++;
      if (dbState(r, 'pubtator') === 'in') inPubtator++;
    });
    var cards = [
      { k: 'Reports', v: rows.length, sub: rows.length === DATA.length ? 'all variants' : 'of ' + DATA.length + ' total' },
      { k: 'Evo2 pathogenic', v: path, sub: benign + ' benign' },
      { k: 'Articles', v: articles.toLocaleString(), sub: rows.length ? (articles / rows.length).toFixed(1) + ' per variant' : '—' },
      { k: 'Variant-specific', v: mention.toLocaleString(), sub: withMention + ' variants with ≥1' },
      { k: 'In HGMD', v: inHgmd, sub: inPubtator + ' in PubTator' },
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
    if (state.model) p.set('model', state.model);
    if (state.pred) p.set('pred', state.pred);
    if (state.tissue) p.set('tissue', state.tissue);
    if (state.hgmd) p.set('hgmd', state.hgmd);
    if (state.pubtator) p.set('pubtator', state.pubtator);
    if (state.mention) p.set('mention', '1');
    if (state.inference) p.set('inference', '1');
    if (state.sort !== 'pos' || state.dir !== 1) p.set('sort', state.sort + ':' + (state.dir === 1 ? 'asc' : 'desc'));
    var s = p.toString();
    history.replaceState(null, '', s ? '#' + s : location.pathname);
  }

  function readHash() {
    var p = new URLSearchParams(location.hash.replace(/^#/, ''));
    state.q = p.get('q') || '';
    state.model = p.get('model') || '';
    state.pred = p.get('pred') || '';
    state.tissue = p.get('tissue') || '';
    state.hgmd = p.get('hgmd') || '';
    state.pubtator = p.get('pubtator') || '';
    state.mention = p.get('mention') === '1';
    state.inference = p.get('inference') === '1';
    var sort = p.get('sort');
    if (sort) {
      var bits = sort.split(':');
      if (SORTERS[bits[0]]) { state.sort = bits[0]; state.dir = bits[1] === 'desc' ? -1 : 1; }
    }
    els.q.value = state.q;
    els.model.value = state.model;
    els.pred.value = state.pred;
    els.tissue.value = state.tissue;
    els.hgmd.value = state.hgmd;
    els.pubtator.value = state.pubtator;
    els.mention.checked = state.mention;
    els.inference.checked = state.inference;
  }

  function render() {
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
  ['model', 'pred', 'tissue', 'hgmd', 'pubtator'].forEach(function (k) {
    els[k].addEventListener('change', function () { state[k] = els[k].value; render(); });
  });
  ['mention', 'inference'].forEach(function (k) {
    els[k].addEventListener('change', function () { state[k] = els[k].checked; render(); });
  });
  els.reset.addEventListener('click', function () {
    state = {
      q: '', model: '', pred: '', tissue: '', mention: false, inference: false,
      hgmd: '', pubtator: '', sort: 'pos', dir: 1
    };
    els.q.value = ''; els.model.value = ''; els.pred.value = ''; els.tissue.value = '';
    els.hgmd.value = ''; els.pubtator.value = '';
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

  fillModels();
  readHash();
  render();
  var built = document.getElementById('built');
  if (built && META.built) built.textContent = META.built;
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

    <div class="toolbar">
        <input type="search" id="q" placeholder="Search variant, position, disease, band, ClinVar ID…" autocomplete="off">
        <select id="f-model"><option value="">All models</option></select>
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
        <select id="f-hgmd" title="Whether this exact variant is catalogued in HGMD">
            <option value="">HGMD: any</option>
            <option value="in">In HGMD</option>
            <option value="out">Not in HGMD</option>
        </select>
        <select id="f-pubtator" title="Whether this exact variant is mentioned in the PubTator-mined corpus">
            <option value="">PubTator: any</option>
            <option value="in">In PubTator</option>
            <option value="out">Not in PubTator</option>
        </select>
        <label class="toggle"><input type="checkbox" id="f-mention"> Mentions variant</label>
        <label class="toggle"><input type="checkbox" id="f-inference"> Has inference</label>
        <button class="reset" id="reset">Reset</button>
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
<script src="assets/site.js"></script>
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


def build(outputs: Path, site: Path, include: str, inline_images: bool) -> Dict:
    from datetime import datetime

    reports_dir = site / "reports"
    assets_dir = site / "assets"
    img_dir = assets_dir / "img"
    data_dir = site / "data"
    for d in (reports_dir, assets_dir, img_dir, data_dir):
        d.mkdir(parents=True, exist_ok=True)

    files = select_reports(outputs, include)
    if not files:
        raise SystemExit(f"No reports matched --include {include} in {outputs}")
    logger.info("found %d report(s) in %s", len(files), outputs)

    parsed: List[Dict] = []
    sources: Dict[str, str] = {}
    for path in files:
        meta = parse_report(path)
        if meta is None:
            continue
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

    index_records = [
        {
            k: meta[k]
            for k in (
                "variant_id", "model", "chrom", "pos", "ref", "alt", "generated",
                "clinvar_id", "evo2", "n_articles", "n_mention", "tissue",
                "mention_tissue", "dominant_tissue", "diseases", "inference", "db", "file",
            )
        }
        for meta in parsed
    ]
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
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(message)s",
    )

    stats = build(args.outputs, args.site, args.include, args.inline_images)
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
