"""
HTML Report Generator Module
Compiles all results into a comprehensive HTML report (Updated for Tissue Support)
"""

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from datetime import datetime
import base64
import html

from config import (
    REPORT_TITLE,
    REPORT_HEADER_COLOR,
    REPORT_ACCENT_COLOR,
    ENABLE_TISSUE_STATISTICS
)

logger = logging.getLogger(__name__)

# Tissue 對應顏色（用於 article border 與 statistics block 顏色塊）
TISSUE_COLORS = {
    'Cardiac':       '#e74c3c',  # 紅
    'Skeletal':      '#f39c12',  # 橘
    'Both':          '#9b59b6',  # 紫
    'Not specified': '#95a5a6',  # 灰
}

# 用來推斷文獻主要疾病的 disease patterns (與 local_clinical_extractor 同步)
# 格式: (regex, canonical_disease_name, tissue)
DISEASE_PATTERNS_FOR_INFERENCE: List[Tuple[str, str, str]] = [
    (r'\bDCM\b|dilated cardiomyopathy',                  'Dilated Cardiomyopathy (DCM)',           'Cardiac'),
    (r'\bHCM\b|hypertrophic cardiomyopathy',             'Hypertrophic Cardiomyopathy (HCM)',      'Cardiac'),
    (r'\bARVC\b|\bACM\b|arrhythmogenic.*cardiomyopathy', 'Arrhythmogenic Cardiomyopathy (ARVC)',   'Cardiac'),
    (r'\bLVNC\b|left ventricular non.?compaction',       'Left Ventricular Noncompaction (LVNC)',  'Cardiac'),
    (r'\bRCM\b|restrictive cardiomyopathy',              'Restrictive Cardiomyopathy (RCM)',       'Cardiac'),
    (r'\bheart failure\b',                               'Heart Failure',                          'Cardiac'),
    (r'\batrial fibrillation\b',                         'Atrial Fibrillation',                    'Cardiac'),
    (r'\bcardiomyopathy\b',                              'Cardiomyopathy (unspecified)',           'Cardiac'),
    (r'\bTMD\b|tibial muscular dystrophy',               'Tibial Muscular Dystrophy (TMD)',        'Skeletal'),
    (r'\bLGMD\b|limb.girdle muscular dystrophy',         'Limb-Girdle Muscular Dystrophy (LGMD)',  'Skeletal'),
    (r'\bHMERF\b|hereditary myopathy with early respiratory failure',
                                                          'HMERF',                                 'Skeletal'),
    (r'\bcentronuclear myopathy\b|\bCNM\b',               'Centronuclear Myopathy (CNM)',          'Skeletal'),
    (r'\bcongenital myopathy\b',                          'Congenital Myopathy',                    'Skeletal'),
    (r'\bmuscular dystrophy\b',                           'Muscular Dystrophy',                     'Skeletal'),
    (r'\bmyopathy\b',                                     'Myopathy',                                'Skeletal'),
    (r'\bcardioskeletal myopathy\b',                      'Cardioskeletal Myopathy',                 'Both'),
    (r'\btitinopath',                                     'Titinopathy',                             'Both'),
]

# Tissue 排序順序（Literature Review 中文獻顯示順序）
TISSUE_SORT_ORDER = {
    'Cardiac': 0,
    'Skeletal': 1,
    'Both': 2,
    'Not specified': 3,
}


class HTMLReportGenerator:
    """Generate comprehensive HTML reports"""
    
    def __init__(self):
        pass
    
    def generate_report(
        self,
        variant_info: Dict,
        evo2_result: Optional[Dict],
        pubmed_results: List[Dict],
        image_path: Path,
        output_path: Path = None,
        transcript_intervals_path: Optional[Path] = None,
        clinvar_info: Optional[Dict] = None,
        aggregated_stats: Optional[Dict] = None, # Kept for older main.py compatibility but ignored below
        mechanistic_result: Optional[Dict] = None,  # Track B：無文獻時的類比推理結果
        db_annotations: Optional[Any] = None,  # HGMD / PubTator 查表命中（DatabaseHits 或 dict）
    ) -> Path:
        logger.info("Generating HTML report...")
        
        image_base64 = self._encode_image(image_path)
        
        html_content = self._generate_html(
            variant_info,
            evo2_result,
            pubmed_results,
            image_base64,
            clinvar_info,
            mechanistic_result,
            db_annotations,
        )
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        logger.info(f"HTML report saved to: {output_path}")
        return output_path
    
    def _encode_image(self, image_path: Path) -> str:
        try:
            with open(image_path, 'rb') as f:
                image_data = f.read()
            return base64.b64encode(image_data).decode('utf-8')
        except Exception as e:
            logger.error(f"Error encoding image: {e}")
            return ""
    
    def _generate_html(
        self,
        variant_info: Dict,
        evo2_result: Optional[Dict],
        pubmed_results: List[Dict],
        image_base64: str,
        clinvar_info: Optional[Dict] = None,
        mechanistic_result: Optional[Dict] = None,
        db_annotations: Optional[Any] = None,
    ) -> str:
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{REPORT_TITLE} - {variant_info['variant_id']}</title>
    <style>
        {self._get_css()}
    </style>
</head>
<body>
    <div class="container">
        {self._generate_header(variant_info, clinvar_info)}
        {self._generate_summary_section(variant_info, evo2_result)}
        {self._generate_database_section(db_annotations)}
        {self._generate_evo2_section(evo2_result)}
        {self._generate_image_section(image_base64, variant_info)}
        {self._generate_inference_section(mechanistic_result)}
        {self._generate_pubmed_section(pubmed_results, clinvar_info)}
        {self._generate_footer()}
    </div>
</body>
</html>"""
    
    def _get_css(self) -> str:
        return f"""
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0, 0, 0, 0.2);
            overflow: hidden;
        }}
        header {{ background: {REPORT_HEADER_COLOR}; color: white; padding: 30px; text-align: center; }}
        header h1 {{ font-size: 2.5em; margin-bottom: 10px; }}
        header .subtitle {{ font-size: 1.2em; opacity: 0.9; }}
        header .clinvar-link {{ color: #7ecfff; text-decoration: underline dotted; }}
        header .clinvar-link:hover {{ color: #ffffff; text-decoration: underline; }}
        .section {{ padding: 30px; border-bottom: 1px solid #eee; }}
        .section h2 {{
            color: {REPORT_HEADER_COLOR};
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid {REPORT_ACCENT_COLOR};
            font-size: 1.8em;
        }}
        .info-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin: 20px 0; }}
        .info-card {{ background: #f8f9fa; padding: 20px; border-radius: 8px; border-left: 4px solid {REPORT_ACCENT_COLOR}; }}
        .info-card .label {{ font-weight: bold; color: #666; font-size: 0.9em; }}
        .info-card .value {{ font-size: 1.2em; color: #333; }}
        .prediction-box {{ background: #f8f9fa; padding: 25px; border-radius: 10px; margin: 20px 0; border: 2px solid #ddd; }}
        .prediction-pathogenic {{ border-color: #e74c3c; background: #fee; }}
        .prediction-benign {{ border-color: #27ae60; background: #efe; }}
        .prediction-result {{ font-size: 2em; font-weight: bold; text-align: center; margin-bottom: 15px; }}
        .prediction-pathogenic .prediction-result {{ color: #e74c3c; }}
        .prediction-benign .prediction-result {{ color: #27ae60; }}
        .score-display {{ display: flex; justify-content: space-around; margin-top: 20px; }}
        .score-value {{ font-size: 1.5em; font-weight: bold; color: {REPORT_ACCENT_COLOR}; }}
        .image-container {{ text-align: center; margin: 20px 0; }}
        .image-container img {{ max-width: 100%; height: auto; border-radius: 8px; box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1); }}
        .pubmed-article {{
            background: #f8f9fa; padding: 20px; margin-bottom: 20px; border-radius: 8px;
            border-left: 6px solid {REPORT_ACCENT_COLOR};
        }}
        .mention-badge {{
            display: inline-block; padding: 2px 8px; border-radius: 4px;
            font-size: 0.8em; font-weight: bold; margin-left: 8px;
            vertical-align: middle;
        }}
        .mention-badge.mention-true {{ background: #27ae60; color: white; }}
        .mention-badge.mention-false {{ background: #95a5a6; color: white; }}
        .source-badge {{
            display: inline-block; padding: 2px 6px; border-radius: 4px;
            font-size: 0.75em; font-weight: bold; margin-left: 6px;
            vertical-align: middle; background: #34495e; color: white;
        }}
        .article-title {{ font-size: 1.2em; font-weight: bold; color: #2c3e50; margin-bottom: 10px; }}
        .article-title a {{ color: {REPORT_ACCENT_COLOR}; text-decoration: none; }}
        .article-meta {{ color: #666; font-size: 0.9em; margin-bottom: 10px; }}
        .article-info {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; margin-top: 15px; }}
        .article-info-item {{ background: white; padding: 10px; border-radius: 5px; }}
        .article-info-label {{ font-weight: bold; color: #666; font-size: 0.85em; }}
        .alert {{ padding: 15px; margin: 20px 0; border-radius: 5px; background: #fff3cd; border: 1px solid #ffc107; color: #856404; }}
        .alert-info {{ background: #d1ecf1; border-color: #bee5eb; color: #0c5460; }}
        footer {{ background: #f8f9fa; padding: 20px; text-align: center; color: #666; font-size: 0.9em; }}
        /* Track B: mechanistic inference（琥珀色調，明顯區隔於 evidence 軌）*/
        .inference-section {{ background: #fffdf7; border-left: 8px solid #f0ad4e; }}
        .inference-section h2 {{ color: #8a5a00; border-bottom-color: #f0ad4e; }}
        .inference-badge {{
            display: inline-block; padding: 3px 10px; border-radius: 4px;
            font-size: 0.5em; font-weight: bold; vertical-align: middle;
            background: #f0ad4e; color: #4a3200; letter-spacing: 0.5px; margin-left: 10px;
        }}
        .inference-note {{
            background: #fcf3e3; border: 1px dashed #e0a44a; color: #7a5a1a;
            padding: 14px 16px; border-radius: 6px; margin-bottom: 18px; font-size: 0.95em;
        }}
        .inference-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 12px 0; }}
        .inference-card {{ background: #fff; padding: 16px; border-radius: 8px; border-left: 4px solid #f0ad4e; }}
        .inference-card .label {{ font-weight: bold; color: #a07a30; font-size: 0.85em; }}
        .inference-card .value {{ font-size: 1.15em; color: #333; }}
        .inference-rationale {{ background: #fff; padding: 16px; border-radius: 8px; margin-top: 16px; border: 1px solid #f0e0c0; }}
        .inference-rationale .label {{ font-weight: bold; color: #a07a30; font-size: 0.85em; margin-bottom: 6px; }}
        .inference-rationale p {{ color: #444; font-style: italic; }}
        .analog-table {{ width: 100%; border-collapse: collapse; margin-top: 8px; font-size: 0.9em; }}
        .analog-table th {{ background: #f7e9cf; color: #7a5a1a; text-align: left; padding: 8px 10px; }}
        .analog-table td {{ padding: 7px 10px; border-bottom: 1px solid #f0e6d2; }}
        .analog-table tr:nth-child(even) {{ background: #fffbf2; }}
        .analog-table a {{ color: {REPORT_ACCENT_COLOR}; text-decoration: none; }}
        {self._get_database_css()}
        """

    def _get_database_css(self) -> str:
        """Curated Database Annotations 區塊的樣式。

        獨立成一個方法，讓 backfill_database_annotations.py 能把同一份樣式
        注入既有報告，不必重跑整個 pipeline。
        """
        return f"""
        /* Curated database lookup（綠/藍調：查表得到的硬紀錄，區隔於琥珀色的推論）*/
        .dbann-section {{ background: #f7fbf8; border-left: 8px solid #2e7d32; }}
        .dbann-section h2 {{ color: #1b5e20; border-bottom-color: #2e7d32; }}
        .dbann-badge {{
            display: inline-block; padding: 3px 10px; border-radius: 4px;
            font-size: 0.5em; font-weight: bold; vertical-align: middle;
            background: #2e7d32; color: #fff; letter-spacing: 0.5px; margin-left: 10px;
        }}
        .dbann-note {{
            background: #e8f5e9; border: 1px solid #a5d6a7; color: #1b5e20;
            padding: 14px 16px; border-radius: 6px; margin-bottom: 18px; font-size: 0.95em;
        }}
        /* 兩邊都沒命中時整段轉為中性灰，避免綠色被誤讀成「有紀錄」 */
        .dbann-section-empty {{ background: #fafafa; border-left-color: #bdbdbd; }}
        .dbann-section-empty h2 {{ color: #555; border-bottom-color: #bdbdbd; }}
        .dbann-section-empty .dbann-badge {{ background: #9e9e9e; }}
        .dbann-section-empty .dbann-note {{ background: #f5f5f5; border-color: #e0e0e0; color: #555; }}
        .dbann-headline {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; margin-bottom: 8px; }}
        .dbann-db {{ background: #fff; padding: 18px; border-radius: 8px; border-left: 6px solid #2e7d32; }}
        .dbann-hgmd {{ border-left-color: #2e7d32; }}
        .dbann-pubtator {{ border-left-color: #1565c0; }}
        .dbann-empty {{ border-left-color: #bdbdbd; background: #fafafa; }}
        .dbann-db-head {{ font-weight: bold; color: #555; font-size: 0.9em; margin-bottom: 8px; }}
        .dbann-disease {{ font-size: 1.35em; font-weight: bold; color: #222; line-height: 1.35; }}
        .dbann-disease.dbann-none {{ color: #9e9e9e; font-weight: normal; }}
        .dbann-db-meta {{ color: #777; font-size: 0.85em; margin-top: 8px; }}
        .dbann-count {{ color: #777; font-size: 0.7em; font-weight: normal; }}
        .dbann-hit {{
            display: inline-block; padding: 2px 7px; border-radius: 3px; margin-left: 6px;
            font-size: 0.85em; background: #2e7d32; color: #fff; letter-spacing: 0.3px;
        }}
        .dbann-pubtator .dbann-hit {{ background: #1565c0; }}
        .dbann-miss {{
            display: inline-block; padding: 2px 7px; border-radius: 3px; margin-left: 6px;
            font-size: 0.85em; background: #bdbdbd; color: #fff; letter-spacing: 0.3px;
        }}
        .dbann-subhead {{ color: #1b5e20; margin: 22px 0 8px; }}
        .dbann-table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; background: #fff; }}
        .dbann-table th {{ background: #e8f5e9; color: #1b5e20; text-align: left; padding: 8px 10px; }}
        .dbann-table td {{ padding: 7px 10px; border-bottom: 1px solid #e6efe7; vertical-align: top; }}
        .dbann-table a {{ color: {REPORT_ACCENT_COLOR}; text-decoration: none; }}
        .dbann-table tr.dbann-evidence td {{
            background: #fbfdfb; color: #555; font-size: 0.92em; font-style: italic;
            border-bottom: 2px solid #e6efe7; padding-top: 0;
        }}
        .dbann-evidence-label {{ font-style: normal; font-weight: bold; color: #1b5e20; }}
        .dbann-details {{ margin-top: 20px; }}
        .dbann-details summary {{ cursor: pointer; color: #33691e; font-size: 0.9em; padding: 8px 0; }}
        """

    def _generate_header(self, variant_info: Dict, clinvar_info: Optional[Dict] = None) -> str:
        variant_id = variant_info['variant_id']
        clinvar_url = clinvar_info.get('clinvar_url') if clinvar_info else None

        if clinvar_url:
            variant_display = (
                f'Variant: <a href="{html.escape(clinvar_url)}" '
                f'target="_blank" rel="noopener" class="clinvar-link">'
                f'{html.escape(variant_id)}</a>'
            )
        else:
            variant_display = f"Variant: {html.escape(variant_id)}"

        return (
            f'<header>'
            f'<h1>{REPORT_TITLE}</h1>'
            f'<div class="subtitle">'
            f'{variant_display}<br>'
            f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
            f'</div>'
            f'</header>'
        )

    def _generate_summary_section(self, variant_info: Dict, evo2_result: Optional[Dict]) -> str:
        pred = evo2_result['prediction'].upper() if evo2_result and evo2_result.get('success') else "N/A"
        return f"""<section class="section"><h2>Variant Summary</h2><div class="info-grid">
            <div class="info-card"><div class="label">Chromosome</div><div class="value">{variant_info['chrom']}</div></div>
            <div class="info-card"><div class="label">Position</div><div class="value">{variant_info['pos']:,}</div></div>
            <div class="info-card"><div class="label">Ref/Alt</div><div class="value">{variant_info['ref']}/{variant_info['alt']}</div></div>
            <div class="info-card"><div class="label">Prediction</div><div class="value">{pred}</div></div>
        </div></section>"""

    # ── Curated database annotations (HGMD / PubTator) ────────────────────

    @staticmethod
    def _pubmed_link(pmid: Any, label: Optional[str] = None) -> str:
        pmid_str = str(pmid or "").strip()
        if not pmid_str or pmid_str.lower() in ("none", "nan", "null", ""):
            return html.escape(str(label or "—"))
        text = html.escape(str(label or f"PMID {pmid_str}"))
        return (
            f'<a href="https://pubmed.ncbi.nlm.nih.gov/{html.escape(pmid_str)}/" '
            f'target="_blank" rel="noopener">{text}</a>'
        )

    @staticmethod
    def _db_cell(value: Any) -> str:
        text = str(value or "").strip()
        if not text or text.lower() in ("none", "nan", "null"):
            return "—"
        return html.escape(text)

    @staticmethod
    def _hgmd_hgvs_cell(record: Dict) -> str:
        """HGMD 的 hgvs 欄位是不帶前綴的 cDNA 寫法（如 82381delG），補上 c. 並註明轉錄本。"""
        raw = str(record.get("hgvs_c") or "").strip()
        if not raw or raw.lower() in ("none", "nan", "null"):
            return "—"
        label = raw if raw[:2].lower() in ("c.", "p.", "g.", "m.", "n.") else f"c.{raw}"
        refseq = str(record.get("refseq") or "").strip()
        if refseq and refseq.lower() not in ("none", "nan", "null"):
            return f'<span title="{html.escape(refseq)}">{html.escape(label)}</span>'
        return html.escape(label)

    @staticmethod
    def _clinvar_sig_text(value: Any) -> str:
        """ClinVar 的分類是底線串接（Conflicting_classifications_of_pathogenicity），改成可讀文字。"""
        text = str(value or "").strip()
        if not text or text.lower() in ("none", "nan", "null"):
            return "—"
        return html.escape(text.replace("_", " ").replace("/", " / "))

    def _generate_database_section(self, db_annotations: Optional[Any]) -> str:
        """把 HGMD / PubTator 對這顆變異登記的疾病直接列出來。

        這一軌是「查表」而非推論：命中就照實呈現資料庫原文，
        所以視覺上用綠/藍（hard record），與 Track B 的琥珀色推論明確區隔。
        """
        if db_annotations is None:
            return ""

        data = db_annotations.to_dict() if hasattr(db_annotations, "to_dict") else dict(db_annotations)

        hgmd = data.get("hgmd_exact") or []
        pubtator = data.get("pubtator_exact") or []
        hgmd_pos = data.get("hgmd_same_position") or []
        pubtator_pos = data.get("pubtator_same_position") or []

        # ── HGMD 摘要卡
        if hgmd:
            diseases = data.get("hgmd_diseases") or []
            disease_text = " · ".join(html.escape(str(d)) for d, _ in diseases) or "Not specified"
            tags = sorted({str(r.get("tag")) for r in hgmd if r.get("tag")})
            mutypes = sorted({str(r.get("mutype")) for r in hgmd if r.get("mutype")})
            inherit = sorted({str(r.get("inheritance")) for r in hgmd if r.get("inheritance")})
            meta_bits = [b for b in (
                "/".join(tags), "/".join(mutypes), "/".join(inherit)
            ) if b]
            hgmd_card = f"""
                <div class="dbann-db dbann-hgmd">
                    <div class="dbann-db-head">HGMD <span class="dbann-hit">EXACT VARIANT MATCH</span></div>
                    <div class="dbann-disease">{disease_text}</div>
                    <div class="dbann-db-meta">{html.escape(' · '.join(meta_bits)) if meta_bits else '&nbsp;'}</div>
                </div>"""
        else:
            note = (
                f"No entry for this variant ({len(hgmd_pos)} other HGMD variant(s) at the same position)"
                if hgmd_pos else "This variant is not catalogued in HGMD"
            )
            hgmd_card = f"""
                <div class="dbann-db dbann-empty">
                    <div class="dbann-db-head">HGMD <span class="dbann-miss">NO MATCH</span></div>
                    <div class="dbann-disease dbann-none">—</div>
                    <div class="dbann-db-meta">{html.escape(note)}</div>
                </div>"""

        # ── PubTator 摘要卡
        if pubtator:
            diseases = data.get("pubtator_diseases") or []
            disease_text = " · ".join(
                f"{html.escape(str(d))}" + (f" <span class=\"dbann-count\">×{n}</span>" if n > 1 else "")
                for d, n in diseases
            ) or "Not specified"
            pmids = sorted({str(r.get("pmid")) for r in pubtator if r.get("pmid")})
            paths = sorted({str(r.get("pathogenicity")) for r in pubtator if r.get("pathogenicity")})
            meta = f"{len(pubtator)} report(s) across {len(pmids)} publication(s)"
            if paths:
                meta += f" · {'/'.join(paths)}"
            pubtator_card = f"""
                <div class="dbann-db dbann-pubtator">
                    <div class="dbann-db-head">PubTator (literature-mined) <span class="dbann-hit">EXACT VARIANT MATCH</span></div>
                    <div class="dbann-disease">{disease_text}</div>
                    <div class="dbann-db-meta">{html.escape(meta)}</div>
                </div>"""
        else:
            note = (
                f"No mention of this variant ({len(pubtator_pos)} other variant(s) reported at the same position)"
                if pubtator_pos else "This variant is not mentioned in the PubTator-mined TTN corpus"
            )
            pubtator_card = f"""
                <div class="dbann-db dbann-empty">
                    <div class="dbann-db-head">PubTator (literature-mined) <span class="dbann-miss">NO MATCH</span></div>
                    <div class="dbann-disease dbann-none">—</div>
                    <div class="dbann-db-meta">{html.escape(note)}</div>
                </div>"""

        # ── HGMD 明細表
        hgmd_table = ""
        if hgmd:
            rows = ""
            for r in hgmd:
                clinvar_cell = self._clinvar_sig_text(r.get("clinvar_clnsig"))
                cid = str(r.get("clinvar_id") or "").split(".")[0]
                if cid and cid.lower() not in ("none", "nan", "null", ""):
                    clinvar_cell = (
                        f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{html.escape(cid)}/" '
                        f'target="_blank" rel="noopener">{clinvar_cell}</a>'
                    )
                ref_label = None
                if r.get("author"):
                    ref_label = f"{r.get('author')} {r.get('year') or ''}".strip()
                rows += f"""
                <tr>
                    <td><b>{self._db_cell(r.get('disease'))}</b></td>
                    <td>{self._hgmd_hgvs_cell(r)}</td>
                    <td>{self._db_cell(r.get('amino') or r.get('descr'))}</td>
                    <td>{self._db_cell(r.get('mutype'))}</td>
                    <td>{self._db_cell(r.get('tag'))}</td>
                    <td>{self._db_cell(r.get('inheritance'))}</td>
                    <td>{clinvar_cell}</td>
                    <td>{self._pubmed_link(r.get('pmid'), ref_label)}</td>
                </tr>"""
            hgmd_table = f"""
            <h4 class="dbann-subhead">HGMD entries for this variant</h4>
            <table class="dbann-table">
                <thead><tr>
                    <th>Disease (as recorded)</th><th>HGVS (c.)</th><th>Protein / change</th>
                    <th>Type</th><th>Class</th><th>Inheritance</th><th>ClinVar</th><th>Reference</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>"""

        # ── PubTator 明細表（新的排前面）
        pubtator_table = ""
        if pubtator:
            def _year_key(rec):
                try:
                    return -int(str(rec.get("year") or 0))
                except ValueError:
                    return 0

            rows = ""
            for r in sorted(pubtator, key=_year_key):
                evidence = str(r.get("evidence_sentence") or "").strip()
                evidence_row = ""
                if evidence and evidence.lower() not in ("none", "nan", "not specified"):
                    if len(evidence) > 400:
                        evidence = evidence[:400] + "…"
                    evidence_row = (
                        f'<tr class="dbann-evidence"><td colspan="7">'
                        f'<span class="dbann-evidence-label">Evidence:</span> '
                        f'&ldquo;{html.escape(evidence)}&rdquo;</td></tr>'
                    )
                rows += f"""
                <tr>
                    <td><b>{self._db_cell(r.get('disease'))}</b></td>
                    <td>{self._db_cell(r.get('hgvs_c') or r.get('variant_text'))}</td>
                    <td>{self._db_cell(r.get('pathogenicity'))}</td>
                    <td>{self._db_cell(r.get('tissue'))}</td>
                    <td>{self._db_cell(r.get('age_onset'))}</td>
                    <td>{self._db_cell(r.get('patient_count'))}</td>
                    <td>{self._pubmed_link(r.get('pmid'), f"{r.get('journal') or 'PubMed'} {r.get('year') or ''}".strip())}</td>
                </tr>{evidence_row}"""
            pubtator_table = f"""
            <h4 class="dbann-subhead">PubTator-mined reports for this variant</h4>
            <table class="dbann-table">
                <thead><tr>
                    <th>Disease (as reported)</th><th>Variant as written</th><th>Pathogenicity</th>
                    <th>Tissue</th><th>Onset</th><th>Patients</th><th>Reference</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>"""

        # ── 同位置但不同 allele（折疊，避免誤讀為同一顆變異）
        same_pos_block = ""
        if hgmd_pos or pubtator_pos:
            rows = ""
            for r in hgmd_pos:
                rows += f"""
                <tr>
                    <td>HGMD</td>
                    <td>{self._hgmd_hgvs_cell(r)}</td>
                    <td>{self._db_cell(r.get('mutype'))}</td>
                    <td><b>{self._db_cell(r.get('disease'))}</b></td>
                    <td>{self._pubmed_link(r.get('pmid'))}</td>
                </tr>"""
            for r in pubtator_pos:
                rows += f"""
                <tr>
                    <td>PubTator</td>
                    <td>{self._db_cell(r.get('hgvs_c') or r.get('variant_text'))}</td>
                    <td>{self._db_cell(r.get('pathogenicity'))}</td>
                    <td><b>{self._db_cell(r.get('disease'))}</b></td>
                    <td>{self._pubmed_link(r.get('pmid'))}</td>
                </tr>"""
            same_pos_block = f"""
            <details class="dbann-details">
                <summary>Other variants recorded at this same genomic position ({len(hgmd_pos) + len(pubtator_pos)}) — different allele, shown for context only</summary>
                <table class="dbann-table">
                    <thead><tr><th>Source</th><th>Variant</th><th>Type / call</th><th>Disease</th><th>Reference</th></tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </details>"""

        has_hit = bool(hgmd or pubtator)
        if has_hit:
            note = (
                "The diseases below are <b>recorded verbatim</b> in the source databases for this exact "
                "variant (same chromosome, position, reference and alternate allele). This is a database "
                "record, not a prediction."
            )
        else:
            note = (
                "This variant was <b>not found</b> in either curated source, so no database-recorded disease "
                "can be reported for it. Any disease shown elsewhere in this report comes from literature "
                "screening or inference."
            )

        section_class = "dbann-section" if has_hit else "dbann-section dbann-section-empty"
        return f"""
        <section class="section {section_class}">
            <h2>Curated Database Annotations <span class="dbann-badge">HGMD · PUBTATOR LOOKUP</span></h2>
            <div class="dbann-note">{note}</div>
            <div class="dbann-headline">
                {hgmd_card}
                {pubtator_card}
            </div>
            {pubtator_table}
            {hgmd_table}
            {same_pos_block}
        </section>"""

    def _generate_evo2_section(self, evo2_result: Optional[Dict]) -> str:
        if not evo2_result or not evo2_result.get('success'): return ""
        cls = 'prediction-pathogenic' if evo2_result['prediction'] == 'pathogenic' else 'prediction-benign'
        return f"""<section class="section"><h2>Evo2 Prediction</h2>
            <div class="prediction-box {cls}"><div class="prediction-result">{evo2_result['prediction'].upper()}</div>
            <div class="score-display">
                <div><div class="score-value">{evo2_result['delta_score']:.4f}</div><div>Delta Score</div></div>
                <div><div class="score-value">{evo2_result['ref_score']:.4f}</div><div>Ref Score</div></div>
                <div><div class="score-value">{evo2_result['var_score']:.4f}</div><div>Var Score</div></div>
            </div></div></section>"""

    def _generate_image_section(self, image_base64: str, variant_info: Dict) -> str:
        if not image_base64: return ""
        return f"""<section class="section"><h2>Protein Domain Localization</h2>
            <div class="image-container"><img src="data:image/png;base64,{image_base64}" alt="Schematic"></div></section>"""

    def _generate_inference_section(self, mechanistic_result: Optional[Dict]) -> str:
        """Track B：無文獻時的類比推理區塊。

        用明顯不同的視覺（琥珀色 + INFERENCE badge）與 evidence 軌區隔，
        清楚標示「推論，非文獻證據」。只有 Track B 觸發時才顯示。
        """
        if not mechanistic_result:
            return ""

        feats = mechanistic_result.get("features", {}) or {}
        inf = mechanistic_result.get("inference", {}) or {}

        disease = html.escape(str(inf.get("inferred_disease", "Uncertain")))
        tissue = html.escape(str(inf.get("inferred_tissue", "Uncertain")))
        confidence = html.escape(str(inf.get("confidence", "Low")))
        rationale = html.escape(str(inf.get("rationale", "")))
        source = html.escape(str(inf.get("source", "")))

        exact = feats.get("hgmd_exact_match") or {}
        exact_banner = ""
        if exact:
            ex_label = html.escape(str(exact.get("prot") or exact.get("dna") or ""))
            ex_phen = html.escape(str(exact.get("phen") or ""))
            ex_tissue = html.escape(str(exact.get("phenotype") or ""))
            exact_banner = (
                f'<div class="inference-note" style="background:#e6f4ea;border-color:#34a853;color:#1e6b34;">'
                f'<b>Direct HGMD match:</b> this exact variant ({ex_label}) is catalogued in HGMD as '
                f'pathogenic for <b>{ex_phen}</b> ({ex_tissue}). The prediction below is anchored to this '
                f'HGMD annotation (it reached this step only because ClinVar had no cited article).</div>'
            )

        band = html.escape(str(feats.get("band", "N/A")))
        vclass = html.escape(str(feats.get("vclass", "N/A")))
        vgroup = html.escape(str(feats.get("vgroup", "N/A")))
        isoforms = html.escape(", ".join(feats.get("isoforms", []) or []) or "N/A")
        iso_hint = feats.get("isoform_tissue_hint", {}) or {}
        iso_hint_text = html.escape(str(iso_hint.get("text", "")))

        # differential diagnoses
        differential = inf.get("differential_diagnoses") or []
        differential_html = ""
        if differential:
            items = " · ".join(html.escape(str(d)) for d in differential if str(d).strip())
            if items:
                differential_html = (
                    f'<div class="inference-card"><div class="label">Differential</div>'
                    f'<div class="value">{items}</div></div>'
                )

        evo2 = feats.get("evo2", {}) or {}
        evo2_disp = "N/A"
        if evo2:
            ds = evo2.get("delta_score")
            ds_str = f"{ds:.4f}" if isinstance(ds, (int, float)) else str(ds)
            level = evo2.get("confidence_level")
            level_str = f", {html.escape(str(level))} conf." if level else ""
            evo2_disp = f"{html.escape(str(evo2.get('prediction', 'N/A')).upper())} (Δ={html.escape(ds_str)}{level_str})"

        # 類比變異表
        analog_rows = ""
        for a in feats.get("analogs", []) or []:
            cid = a.get("clinvar_id")
            variant_label = html.escape(str(a.get("prot") or a.get("dna") or f"{a.get('ref')}>{a.get('alt')}"))
            if cid and str(cid).lower() not in ("none", "nan"):
                cid_clean = str(cid).split(".")[0]
                variant_cell = (
                    f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{html.escape(cid_clean)}/" '
                    f'target="_blank" rel="noopener">{variant_label}</a>'
                )
            else:
                variant_cell = variant_label
            band_only = ' <span style="color:#b26a00;font-size:0.8em;">(band-only)</span>' if a.get("band_only_match") else ""
            analog_rows += f"""
                <tr>
                    <td>{variant_cell}</td>
                    <td>{html.escape(str(a.get('vclass', '')))}</td>
                    <td>{html.escape(str(a.get('phen', '')))}{band_only}</td>
                    <td>{html.escape(str(a.get('phenotype', '')))}</td>
                    <td>{a.get('distance', '')} bp</td>
                </tr>"""

        analog_table = ""
        if analog_rows:
            analog_table = f"""
            <h4 style="color:#8a5a00; margin:18px 0 8px;">Analogous known pathogenic variants (HGMD · same band + same type · nearest first)</h4>
            <table class="analog-table">
                <thead><tr>
                    <th>Variant</th><th>Type</th><th>Disease</th><th>Tissue</th><th>Distance</th>
                </tr></thead>
                <tbody>{analog_rows}</tbody>
            </table>"""

        # band 統計先驗
        prior = feats.get("band_prior", {}) or {}
        counts = prior.get("counts", {}) or {}
        prior_line = ""
        if prior.get("total"):
            prior_line = (
                f'<p style="color:#8a6d3b; margin-top:12px; font-size:0.92em;">'
                f'Band-level prior for <b>{band}</b>: '
                f'Cardiac {counts.get("cardiac", 0)}, Skeletal {counts.get("skeletal", 0)}, '
                f'Both {counts.get("both", 0)} '
                f'(n={prior.get("total", 0)} known pathogenic variants).</p>'
            )

        return f"""
        <section class="section inference-section">
            <h2>Mechanistic Inference <span class="inference-badge">INFERENCE · NO DIRECT LITERATURE</span></h2>
            <div class="inference-note">
                No published article was found that directly reports this variant.
                The following is an <b>analogy-based hypothesis</b> derived from Evo2, the variant's structural
                band and type, and known pathogenic TTN variants in the same region — <b>not evidence-based</b>.
            </div>
            {exact_banner}
            <div class="inference-grid">
                <div class="inference-card" style="border-left-width:6px;"><div class="label">Predicted Disease (inference)</div><div class="value" style="font-weight:bold;">{disease}</div></div>
                <div class="inference-card"><div class="label">Predicted Tissue</div><div class="value">{tissue}</div></div>
                <div class="inference-card"><div class="label">Confidence</div><div class="value">{confidence}</div></div>
                {differential_html}
            </div>
            <div class="inference-grid">
                <div class="inference-card"><div class="label">Structural Band</div><div class="value">{band}</div></div>
                <div class="inference-card"><div class="label">Variant Type</div><div class="value">{vclass} ({vgroup})</div></div>
                <div class="inference-card"><div class="label">Evo2</div><div class="value">{evo2_disp}</div></div>
                <div class="inference-card"><div class="label">Affected Isoforms</div><div class="value">{isoforms}</div></div>
            </div>
            {f'<div class="inference-note" style="margin:0 0 12px;"><b>Isoform → tissue:</b> {iso_hint_text}</div>' if iso_hint_text else ''}
            <div class="inference-rationale">
                <div class="label">Rationale</div>
                <p>{rationale}</p>
            </div>
            {analog_table}
            {prior_line}
            <p style="color:#999; font-size:0.8em; margin-top:12px;">Inference source: {source}</p>
        </section>"""

    # Removed _generate_clinical_stats_section completely.

    # ── Helpers ────────────────────────────────────────────────────────────

    @staticmethod
    def _normalize_tissue(tissue_raw: str) -> str:
        """把任意 tissue 字串標準化為 Cardiac / Skeletal / Both / Not specified."""
        if not tissue_raw:
            return 'Not specified'
        t = str(tissue_raw).strip().lower()
        if 'both' in t:
            return 'Both'
        if 'cardiac' in t:
            return 'Cardiac'
        if 'skeletal' in t:
            return 'Skeletal'
        return 'Not specified'

    @staticmethod
    def _article_mentions_variant(article: Dict) -> bool:
        """判斷文獻是否直接提到輸入變異
        
        判斷依據（任一條成立即視為「直接提到」）:
        - clinical_info.evidence_sentence 不是空、不是 "Not specified"
        - clinical_info.disease 不是 "Not specified"（且不是 JSON parsing failed）
        - patient_count > 0
        """
        info = article.get('clinical_info', {}) or {}
        
        evidence = (info.get('evidence_sentence') or '').strip()
        if evidence and evidence.lower() not in ('not specified', 'n/a', 'none', ''):
            return True
        
        disease = (info.get('disease') or '').strip()
        if disease and disease.lower() not in ('not specified', 'n/a', 'none', '', 'json parsing failed'):
            return True
        
        try:
            if int(info.get('patient_count') or 0) > 0:
                return True
        except (ValueError, TypeError):
            pass
        
        return False

    @staticmethod
    def _infer_main_disease_from_article(article: Dict) -> Tuple[str, str]:
        """從文獻 title + abstract 推斷主要疾病 / tissue (regex match)
        
        Returns: (disease_name, tissue) — 若無匹配回傳 ('Not specified', 'Not specified')
        """
        title = article.get('title') or ''
        abstract = article.get('abstract') or ''
        text = f"{title} {abstract}"
        if not text.strip():
            return ('Not specified', 'Not specified')
        
        for pattern, disease, tissue in DISEASE_PATTERNS_FOR_INFERENCE:
            if re.search(pattern, text, re.IGNORECASE):
                return (disease, tissue)
        
        return ('Not specified', 'Not specified')

    # ── Sections ──────────────────────────────────────────────────────────

    def _generate_pubmed_section(self, pubmed_results: List[Dict], clinvar_info: Optional[Dict] = None) -> str:
        if not pubmed_results:
            return """<section class="section"><h2>Literature Review</h2><div class="alert">No PubMed articles found.</div></section>"""
        
        # Step 1: 對每篇 article 計算 mention_input_variant、最終疾病/tissue
        # 若有提到變異 → 使用 LLM 提取結果
        # 若沒提到變異 → 從 title+abstract 用 regex 推斷主要疾病
        enriched_articles: List[Dict] = []
        for article in pubmed_results:
            info = article.get('clinical_info', {}) or {}
            mention = self._article_mentions_variant(article)
            
            if mention:
                final_disease = info.get('disease') or 'Not specified'
                final_tissue = self._normalize_tissue(info.get('tissue_affected'))
                disease_source = 'extracted_from_variant'
            else:
                # 沒直接提到變異 → 推斷文獻主要疾病
                inferred_disease, inferred_tissue = self._infer_main_disease_from_article(article)
                final_disease = inferred_disease
                final_tissue = inferred_tissue
                disease_source = 'inferred_from_article'
            
            enriched_articles.append({
                'article': article,
                'mention': mention,
                'final_disease': final_disease,
                'final_tissue': final_tissue,
                'disease_source': disease_source,
            })
        
        # Step 2: 按 tissue 排序 (Cardiac → Skeletal → Both → Not specified)
        # 排序的次要 key 用 mention（提到的優先），再來用 year（新的優先）
        def _sort_key(item):
            tissue_order = TISSUE_SORT_ORDER.get(item['final_tissue'], 99)
            mention_order = 0 if item['mention'] else 1  # 提到的優先
            year_str = str(item['article'].get('year') or '0')
            try:
                year_neg = -int(year_str)
            except ValueError:
                year_neg = 0
            return (tissue_order, mention_order, year_neg)
        
        enriched_articles.sort(key=_sort_key)
        
        html_out = f'<section class="section"><h2>Literature Review ({len(pubmed_results)} articles)</h2>'
        
        # Step 3: Tissue Involvement Statistics（雙排）
        if ENABLE_TISSUE_STATISTICS:
            tissue_stats = self._calculate_tissue_statistics(enriched_articles)
            html_out += self._generate_tissue_statistics_block(tissue_stats)
        
        # Step 4: 渲染每篇文獻
        for idx, item in enumerate(enriched_articles, 1):
            article = item['article']
            info = article.get('clinical_info', {}) or {}
            
            final_disease = item['final_disease']
            final_tissue = item['final_tissue']
            mention = item['mention']
            disease_source = item['disease_source']
            
            age = info.get('age_onset', 'Not specified')
            inheritance = info.get('inheritance', 'Not specified')
            reasoning = info.get('reasoning', '')
            evidence = info.get('evidence_sentence', '')

            # AI Analysis Block
            ai_block = ""
            if reasoning or evidence:
                ai_block = f"""
                <div style="margin-top: 15px; padding: 15px; background: white; border-radius: 5px; border: 1px solid #eee;">
                    <div style="font-weight: bold; color: {REPORT_ACCENT_COLOR}; margin-bottom: 5px;">AI Analysis</div>
                    <div style="font-style: italic; color: #555; margin-bottom: 8px;">"{html.escape(reasoning)}"</div>
                    {f'<div style="background: #eef; padding: 5px 10px; border-radius: 4px; font-family: monospace; font-size: 0.9em; color: #333;">Evidence: "{html.escape(evidence)}"</div>' if evidence else ''}
                </div>
                """

            # 來源 badge（保留 LitVar 標籤，但用獨立 badge）
            source = article.get('source') or 'PubMed'
            source_badge = f'<span class="source-badge">{html.escape(str(source))}</span>'
            
            # Mention input variant badge
            mention_class = 'mention-true' if mention else 'mention-false'
            mention_label = 'Mention input variant: True' if mention else 'Mention input variant: False'
            mention_badge = f'<span class="mention-badge {mention_class}">{mention_label}</span>'
            
            # 左側 border 顏色 → 按 tissue
            border_color = TISSUE_COLORS.get(final_tissue, TISSUE_COLORS['Not specified'])
            
            # 如果是推斷的疾病，顯示 (inferred) 註記
            disease_display = html.escape(final_disease)
            if disease_source == 'inferred_from_article' and final_disease != 'Not specified':
                disease_display = f'{html.escape(final_disease)} <span style="color:#888; font-size:0.85em;">(inferred from article)</span>'
            
            html_out += f"""
            <div class="pubmed-article" style="border-left-color: {border_color};">
                <div class="article-title">{idx}. <a href="{article['pubmed_link']}" target="_blank">{html.escape(article['title'])}</a>{source_badge}{mention_badge}</div>
                <div class="article-meta">{html.escape(str(article['journal']))} ({article['year']}) | PMID: {article['pmid']}</div>
                
                <div class="article-info">
                    <div class="article-info-item"><div class="article-info-label">Disease</div><div>{disease_display}</div></div>
                    <div class="article-info-item"><div class="article-info-label">Tissue Affected</div><div>{html.escape(final_tissue)}</div></div>
                    <div class="article-info-item"><div class="article-info-label">Age</div><div>{html.escape(str(age))}</div></div>
                    <div class="article-info-item"><div class="article-info-label">Inheritance</div><div>{html.escape(str(inheritance))}</div></div>
                </div>
                {ai_block}
            </div>
            """
        
        html_out += "</section>"
        return html_out

    def _calculate_tissue_statistics(self, enriched_articles: List[Dict]) -> Dict:
        """計算雙排組織影響統計
        
        Args:
            enriched_articles: List of dicts with keys 'mention', 'final_tissue', 'article', ...
        
        Returns:
            {
              'mentioned':     {'cardiac': N, 'skeletal': N, 'both': N, 'not_specified': N, 'total': N},
              'not_mentioned': {'cardiac': N, 'skeletal': N, 'both': N, 'not_specified': N, 'total': N},
              'total_articles': N
            }
        """
        def _empty():
            return {'cardiac': 0, 'skeletal': 0, 'both': 0, 'not_specified': 0, 'total': 0}
        
        stats = {
            'mentioned': _empty(),
            'not_mentioned': _empty(),
            'total_articles': len(enriched_articles),
        }
        
        tissue_key_map = {
            'Cardiac': 'cardiac',
            'Skeletal': 'skeletal',
            'Both': 'both',
            'Not specified': 'not_specified',
        }
        
        for item in enriched_articles:
            bucket = 'mentioned' if item['mention'] else 'not_mentioned'
            tissue_key = tissue_key_map.get(item['final_tissue'], 'not_specified')
            stats[bucket][tissue_key] += 1
            stats[bucket]['total'] += 1
        
        return stats
    
    def _generate_tissue_statistics_block(self, stats: Dict) -> str:
        """生成雙排 Tissue Involvement Statistics 區塊。
        
        - 上排：文獻直接提到輸入變異（disease 取自 LLM extraction）
        - 下排：文獻沒提到輸入變異（disease 取自文獻主要討論的疾病推斷）
        """
        if stats['total_articles'] == 0:
            return """<div class="alert alert-info">No articles available for tissue statistics.</div>"""
        
        def _row(title: str, subtitle: str, bucket: Dict) -> str:
            total = bucket['total']
            if total == 0:
                return f"""
                <div style="margin-bottom: 12px;">
                    <h4 style="color: {REPORT_HEADER_COLOR}; margin-bottom: 8px;">{title}</h4>
                    <p style="color: #999; font-style: italic;">{subtitle} — No articles in this category.</p>
                </div>
                """
            
            def _pct(n):
                return (n / total * 100) if total > 0 else 0
            
            return f"""
            <div style="margin-bottom: 16px;">
                <h4 style="color: {REPORT_HEADER_COLOR}; margin-bottom: 6px;">{title} <span style="font-weight: normal; color: #666;">— {total} articles</span></h4>
                <p style="color: #777; font-size: 0.9em; margin-bottom: 10px;">{subtitle}</p>
                <div class="info-grid">
                    <div class="info-card" style="border-left-color: {TISSUE_COLORS['Cardiac']};">
                        <div class="label">Cardiac</div>
                        <div class="value">{bucket['cardiac']} ({_pct(bucket['cardiac']):.1f}%)</div>
                    </div>
                    <div class="info-card" style="border-left-color: {TISSUE_COLORS['Skeletal']};">
                        <div class="label">Skeletal</div>
                        <div class="value">{bucket['skeletal']} ({_pct(bucket['skeletal']):.1f}%)</div>
                    </div>
                    <div class="info-card" style="border-left-color: {TISSUE_COLORS['Both']};">
                        <div class="label">Both</div>
                        <div class="value">{bucket['both']} ({_pct(bucket['both']):.1f}%)</div>
                    </div>
                    <div class="info-card" style="border-left-color: {TISSUE_COLORS['Not specified']};">
                        <div class="label">Not Specified</div>
                        <div class="value">{bucket['not_specified']} ({_pct(bucket['not_specified']):.1f}%)</div>
                    </div>
                </div>
            </div>
            """
        
        mentioned_row = _row(
            'Articles that directly mention the input variant',
            'Disease / tissue extracted from the article for the target variant.',
            stats['mentioned'],
        )
        not_mentioned_row = _row(
            'Articles that do NOT directly mention the input variant',
            'Disease / tissue inferred from the main topic discussed in the article.',
            stats['not_mentioned'],
        )
        
        return f"""
        <div style="margin: 20px 0; padding: 20px; background: linear-gradient(135deg, #667eea15 0%, #764ba215 100%); border-radius: 10px; border: 2px solid {REPORT_ACCENT_COLOR};">
            <h3 style="color: {REPORT_HEADER_COLOR}; margin-bottom: 15px; font-size: 1.3em;">Tissue Involvement Statistics</h3>
            <p style="color: #666; margin-bottom: 15px;">Total: {stats['total_articles']} articles
            ({stats['mentioned']['total']} mention the variant, {stats['not_mentioned']['total']} do not)</p>
            {mentioned_row}
            {not_mentioned_row}
        </div>
        """

    def _generate_footer(self) -> str:
        return f"<footer><p>Generated by TTN Variant AI Agent</p></footer>"

    def _get_javascript(self) -> str:
        return ""