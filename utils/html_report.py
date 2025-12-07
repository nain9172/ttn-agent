"""
HTML Report Generator Module
Compiles all results into a comprehensive HTML report (Updated for Tissue Support)
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import base64
import html

from config import (
    REPORT_TITLE,
    REPORT_HEADER_COLOR,
    REPORT_ACCENT_COLOR
)

logger = logging.getLogger(__name__)


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
            clinvar_info
            # aggregated_stats is no longer passed to _generate_html
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
        {self._generate_header(variant_info)}
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
            border-left: 4px solid {REPORT_ACCENT_COLOR};
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

    def _generate_header(self, variant_info: Dict) -> str:
        return f"""<header><h1>{REPORT_TITLE}</h1><div class="subtitle">Variant: {variant_info['variant_id']}<br>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div></header>"""

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

    def _generate_pubmed_section(self, pubmed_results: List[Dict], clinvar_info: Optional[Dict] = None) -> str:
        if not pubmed_results:
            return """<section class="section"><h2>Literature Review</h2><div class="alert">No PubMed articles found.</div></section>"""
        
        html_out = f'<section class="section"><h2>Literature Review ({len(pubmed_results)} articles)</h2>'
        
        for idx, article in enumerate(pubmed_results, 1):
            info = article.get('clinical_info', {})
            
            # Extract fields
            disease = info.get('disease', 'Not specified')
            tissue = info.get('tissue_affected', 'Not specified')
            age = info.get('age_onset', 'Not specified')
            inheritance = info.get('inheritance', 'Not specified')
            reasoning = info.get('reasoning', '')
            evidence = info.get('evidence_sentence', '')

            # AI Analysis Block - kept because the user only asked to remove the *Synthesis* part
            ai_block = ""
            if reasoning or evidence:
                ai_block = f"""
                <div style="margin-top: 15px; padding: 15px; background: white; border-radius: 5px; border: 1px solid #eee;">
                    <div style="font-weight: bold; color: {REPORT_ACCENT_COLOR}; margin-bottom: 5px;">🤖 AI Analysis</div>
                    <div style="font-style: italic; color: #555; margin-bottom: 8px;">"{reasoning}"</div>
                    {f'<div style="background: #eef; padding: 5px 10px; border-radius: 4px; font-family: monospace; font-size: 0.9em; color: #333;">Evidence: "{evidence}"</div>' if evidence else ''}
                </div>
                """

            html_out += f"""
            <div class="pubmed-article">
                <div class="article-title">{idx}. <a href="{article['pubmed_link']}" target="_blank">{article['title']}</a></div>
                <div class="article-meta">{article['journal']} ({article['year']}) | PMID: {article['pmid']}</div>
                
                <div class="article-info">
                    <div class="article-info-item"><div class="article-info-label">Disease</div><div>{disease}</div></div>
                    <div class="article-info-item"><div class="article-info-label">Tissue Affected</div><div>{tissue}</div></div>
                    <div class="article-info-item"><div class="article-info-label">Age Onset</div><div>{age}</div></div>
                    <div class="article-info-item"><div class="article-info-label">Inheritance</div><div>{inheritance}</div></div>
                </div>
                {ai_block}
            </div>
            """
        
        html_out += "</section>"
        return html_out

    def _generate_footer(self) -> str:
        return f"<footer><p>Generated by TTN Variant AI Agent</p></footer>"

    def _get_javascript(self) -> str:
        return ""