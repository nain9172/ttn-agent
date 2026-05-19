"""
HTML Report Generator Module
Compiles all results into a comprehensive HTML report (Updated for Tissue Support)
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
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
        aggregated_stats: Optional[Dict] = None # Kept for older main.py compatibility but ignored below
    ) -> Path:
        logger.info("Generating HTML report...")
        
        image_base64 = self._encode_image(image_path)
        
        html_content = self._generate_html(
            variant_info,
            evo2_result,
            pubmed_results,
            image_base64,
            clinvar_info,
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
        # Removed aggregated_stats: Optional[Dict] = None,
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
        {self._generate_evo2_section(evo2_result)}
        {self._generate_image_section(image_base64, variant_info)}
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