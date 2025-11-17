"""
HTML Report Generator Module
Compiles all results into a comprehensive HTML report
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import base64

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
        aggregated_stats: Optional[Dict] = None
    ) -> Path:
        """
        Generate complete HTML report
        
        Args:
            variant_info: Variant information
            evo2_result: Evo2 prediction results
            pubmed_results: PubMed search results
            image_path: Path to generated protein schematic (now includes transcript intervals)
            output_path: Output HTML file path
            transcript_intervals_path: Deprecated, kept for backward compatibility
            clinvar_info: ClinVar information including conditions and variant impact
            aggregated_stats: Aggregated clinical statistics from LLM extraction
        
        Returns:
            Path to generated HTML report
        """
        logger.info("Generating HTML report...")
        
        # Convert image to base64 for embedding
        image_base64 = self._encode_image(image_path)
        
        # Generate HTML content
        html = self._generate_html(
            variant_info,
            evo2_result,
            pubmed_results,
            image_base64,
            clinvar_info,
            aggregated_stats
        )
        
        # Write to file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        
        logger.info(f"HTML report saved to: {output_path}")
        return output_path
    
    def _encode_image(self, image_path: Path) -> str:
        """Encode image to base64 for embedding in HTML"""
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
        aggregated_stats: Optional[Dict] = None
    ) -> str:
        """Generate complete HTML content"""
        
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
        {self._generate_clinical_stats_section(aggregated_stats) if aggregated_stats else ""}
        {self._generate_pubmed_section(pubmed_results, clinvar_info)}
        {self._generate_footer()}
    </div>
    
    <script>
        {self._get_javascript()}
    </script>
</body>
</html>"""
    
    def _get_css(self) -> str:
        """Generate CSS styles"""
        return f"""
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
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
        
        header {{
            background: {REPORT_HEADER_COLOR};
            color: white;
            padding: 30px;
            text-align: center;
        }}
        
        header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        header .subtitle {{
            font-size: 1.2em;
            opacity: 0.9;
        }}
        
        .section {{
            padding: 30px;
            border-bottom: 1px solid #eee;
        }}
        
        .section:last-child {{
            border-bottom: none;
        }}
        
        .section h2 {{
            color: {REPORT_HEADER_COLOR};
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid {REPORT_ACCENT_COLOR};
            font-size: 1.8em;
        }}
        
        .section h3 {{
            color: {REPORT_ACCENT_COLOR};
            margin-top: 20px;
            margin-bottom: 15px;
            font-size: 1.3em;
        }}
        
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .info-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid {REPORT_ACCENT_COLOR};
        }}
        
        .info-card .label {{
            font-weight: bold;
            color: #666;
            font-size: 0.9em;
            margin-bottom: 5px;
        }}
        
        .info-card .value {{
            font-size: 1.2em;
            color: #333;
        }}
        
        .prediction-box {{
            background: #f8f9fa;
            padding: 25px;
            border-radius: 10px;
            margin: 20px 0;
            border: 2px solid #ddd;
        }}
        
        .prediction-pathogenic {{
            border-color: #e74c3c;
            background: #fee;
        }}
        
        .prediction-benign {{
            border-color: #27ae60;
            background: #efe;
        }}
        
        .prediction-result {{
            font-size: 2em;
            font-weight: bold;
            text-align: center;
            margin-bottom: 15px;
        }}
        
        .prediction-pathogenic .prediction-result {{
            color: #e74c3c;
        }}
        
        .prediction-benign .prediction-result {{
            color: #27ae60;
        }}
        
        .score-display {{
            display: flex;
            justify-content: space-around;
            margin-top: 20px;
        }}
        
        .score-item {{
            text-align: center;
        }}
        
        .score-value {{
            font-size: 1.5em;
            font-weight: bold;
            color: {REPORT_ACCENT_COLOR};
        }}
        
        .score-label {{
            font-size: 0.9em;
            color: #666;
            margin-top: 5px;
        }}
        
        .image-container {{
            text-align: center;
            margin: 20px 0;
        }}
        
        .image-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }}
        
        .pubmed-article {{
            background: #f8f9fa;
            padding: 20px;
            margin-bottom: 20px;
            border-radius: 8px;
            border-left: 4px solid {REPORT_ACCENT_COLOR};
            transition: transform 0.2s, box-shadow 0.2s;
        }}
        
        .pubmed-article:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
        }}
        
        .article-title {{
            font-size: 1.2em;
            font-weight: bold;
            color: #2c3e50;
            margin-bottom: 10px;
        }}
        
        .article-title a {{
            color: {REPORT_ACCENT_COLOR};
            text-decoration: none;
        }}
        
        .article-title a:hover {{
            text-decoration: underline;
        }}
        
        .article-meta {{
            color: #666;
            font-size: 0.9em;
            margin-bottom: 10px;
        }}
        
        .article-info {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 10px;
            margin-top: 15px;
        }}
        
        .article-info-item {{
            background: white;
            padding: 10px;
            border-radius: 5px;
        }}
        
        .article-info-label {{
            font-weight: bold;
            color: #666;
            font-size: 0.85em;
        }}
        
        .article-info-value {{
            color: #333;
            margin-top: 3px;
        }}
        
        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: bold;
            margin-right: 5px;
        }}
        
        .badge-cardiac {{
            background: #e74c3c;
            color: white;
        }}
        
        .badge-skeletal {{
            background: #3498db;
            color: white;
        }}
        
        .badge-both {{
            background: #9b59b6;
            color: white;
        }}
        
        .badge-pathogenic {{
            background: #c0392b;
            color: white;
        }}
        
        .badge-likely-pathogenic {{
            background: #e74c3c;
            color: white;
        }}
        
        .badge-uncertain {{
            background: #95a5a6;
            color: white;
        }}
        
        .badge-likely-benign {{
            background: #27ae60;
            color: white;
        }}
        
        .badge-benign {{
            background: #1e8449;
            color: white;
        }}
        
        .badge-default {{
            background: #7f8c8d;
            color: white;
        }}
        
        footer {{
            background: #f8f9fa;
            padding: 20px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
        }}
        
        .alert {{
            padding: 15px;
            margin: 20px 0;
            border-radius: 5px;
            background: #fff3cd;
            border: 1px solid #ffc107;
            color: #856404;
        }}
        
        .alert-info {{
            background: #d1ecf1;
            border-color: #bee5eb;
            color: #0c5460;
        }}
        
        @media print {{
            body {{
                background: white;
                padding: 0;
            }}
            
            .container {{
                box-shadow: none;
            }}
        }}
        """
    
    def _generate_header(self, variant_info: Dict) -> str:
        """Generate header section"""
        return f"""
        <header>
            <h1>{REPORT_TITLE}</h1>
            <div class="subtitle">
                Variant: {variant_info['variant_id']}<br>
                Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
            </div>
        </header>
        """
    
    def _generate_summary_section(
        self,
        variant_info: Dict,
        evo2_result: Optional[Dict]
    ) -> str:
        """Generate summary section"""
        prediction_text = "N/A"
        if evo2_result and evo2_result.get('success'):
            prediction_text = evo2_result['prediction'].upper()
        
        return f"""
        <section class="section">
            <h2>Variant Summary</h2>
            <div class="info-grid">
                <div class="info-card">
                    <div class="label">Chromosome</div>
                    <div class="value">{variant_info['chrom']}</div>
                </div>
                <div class="info-card">
                    <div class="label">Position (GRCh38)</div>
                    <div class="value">{variant_info['pos']:,}</div>
                </div>
                <div class="info-card">
                    <div class="label">Reference Allele</div>
                    <div class="value">{variant_info['ref']}</div>
                </div>
                <div class="info-card">
                    <div class="label">Alternate Allele</div>
                    <div class="value">{variant_info['alt']}</div>
                </div>
                <div class="info-card">
                    <div class="label">Gene</div>
                    <div class="value">TTN (Titin)</div>
                </div>
                <div class="info-card">
                    <div class="label">Prediction</div>
                    <div class="value">{prediction_text}</div>
                </div>
            </div>
        </section>
        """
    
    def _generate_evo2_section(self, evo2_result: Optional[Dict]) -> str:
        """Generate Evo2 prediction section"""
        if not evo2_result:
            return """
            <section class="section">
                <h2>Evo2 Pathogenicity Prediction</h2>
                <div class="alert">Evo2 prediction was skipped.</div>
            </section>
            """
        
        if not evo2_result.get('success'):
            error_msg = evo2_result.get('error', 'Unknown error')
            return f"""
            <section class="section">
                <h2>Evo2 Pathogenicity Prediction</h2>
                <div class="alert">
                    Error in Evo2 prediction: {error_msg}
                </div>
            </section>
            """
        
        prediction = evo2_result['prediction']
        delta_score = evo2_result['delta_score']
        box_class = 'prediction-pathogenic' if prediction == 'pathogenic' else 'prediction-benign'
        
        return f"""
        <section class="section">
            <h2>Evo2 Pathogenicity Prediction</h2>
            
            <div class="prediction-box {box_class}">
                <div class="prediction-result">
                    {prediction.upper()}
                </div>
                
                <div class="score-display">
                    <div class="score-item">
                        <div class="score-value">{delta_score:.6f}</div>
                        <div class="score-label">Delta Score</div>
                    </div>
                    <div class="score-item">
                        <div class="score-value">{evo2_result['ref_score']:.6f}</div>
                        <div class="score-label">Reference Score</div>
                    </div>
                    <div class="score-item">
                        <div class="score-value">{evo2_result['var_score']:.6f}</div>
                        <div class="score-label">Variant Score</div>
                    </div>
                </div>
            </div>
            
            <div class="alert alert-info">
                <strong>Interpretation:</strong> 
                {'Delta score < 0 indicates the variant is likely pathogenic (deleterious).' 
                 if prediction == 'pathogenic' 
                 else 'Delta score ≥ 0 indicates the variant is likely benign (tolerated).'}
                <br><br>
                <strong>Confidence:</strong> {abs(delta_score):.6f} 
                (larger absolute value = higher confidence)
            </div>
        </section>
        """
    
    def _generate_image_section(self, image_base64: str, variant_info: Dict) -> str:
        """Generate protein schematic section"""
        if not image_base64:
            return """
            <section class="section">
                <h2>Protein Domain Localization</h2>
                <div class="alert">Image generation failed.</div>
            </section>
            """
        
        return f"""
        <section class="section">
            <h2>Protein Domain Localization & Transcript Intervals</h2>
            <p>The figure below shows the TTN protein domain structure (Z-disk, I-band, A-band, M-band) 
               and the genomic intervals for different transcript isoforms, with variant {variant_info['variant_id']} 
               position marked.</p>
            
            <div class="image-container">
                <img src="data:image/png;base64,{image_base64}" 
                     alt="TTN Protein Schematic">
            </div>
            
            <div class="alert alert-info">
                <strong>Note:</strong> The top panel shows the protein domain structure. 
                Each subsequent row shows the genomic intervals (colored blocks) for a specific transcript isoform.
                The red arrow indicates the variant position across all views. All visualizations are aligned 
                by genomic coordinates on chromosome 2 (negative strand).
            </div>
        </section>
        """
    
    def _generate_transcript_intervals_section(self, transcript_intervals_base64: str, variant_info: Dict) -> str:
        """Generate transcript intervals section"""
        if not transcript_intervals_base64:
            return ""
        
        return f"""
        <section class="section">
            <h2>Transcript Intervals (Genomic Coordinates)</h2>
            <p>The figure below shows the genomic intervals for different TTN transcript isoforms, 
               with the variant {variant_info['variant_id']} position marked.</p>
            
            <div class="image-container">
                <img src="data:image/png;base64,{transcript_intervals_base64}" 
                     alt="TTN Transcript Intervals">
            </div>
            
            <div class="alert alert-info">
                <strong>Note:</strong> Each colored block represents a genomic interval (exon or exonic region) 
                for the transcript. The intervals are shown on the genomic coordinate scale (chr2), 
                with the TTN gene located on the negative strand.
            </div>
        </section>
        """
    
    def _generate_pubmed_section(self, pubmed_results: List[Dict], clinvar_info: Optional[Dict] = None) -> str:
        """Generate PubMed literature section"""
        if not pubmed_results:
            return """
            <section class="section">
                <h2>Literature Review (PubMed)</h2>
                <div class="alert alert-info">
                    <strong>無 PubMed 文獻</strong>
                    <p>本報告僅包含 ClinVar 提供的精確文獻。該變異在 ClinVar 中未關聯任何 PubMed 文章。</p>
                    <p style="margin-top: 10px;">這可能表示：</p>
                    <ul style="margin-left: 20px;">
                        <li>該變異為新發現的變異</li>
                        <li>該變異尚未被詳細研究</li>
                        <li>相關文獻未被 ClinVar 收錄</li>
                    </ul>
                </div>
            </section>
            """
        
        # Add summary statistics for ages if available
        age_summary = ""
        if clinvar_info:
            age_summary = self._generate_age_summary(pubmed_results)
        
        articles_html = ""
        for idx, article in enumerate(pubmed_results, 1):
            phenotype_badge = self._get_phenotype_badge(article.get('phenotype', 'Not specified'))
            
            # Get clinical info if available
            clinical_info = article.get('clinical_info', {})
            
            # Format inheritance with confidence
            inheritance_display = article.get('inheritance', 'Not specified')
            if clinical_info.get('inheritance'):
                inheritance_data = clinical_info['inheritance']
                pattern = inheritance_data.get('pattern', 'not_specified')
                confidence = inheritance_data.get('confidence', '')
                if pattern != 'not_specified':
                    pattern_names = {
                        'autosomal_dominant': 'Autosomal Dominant (顯性)',
                        'autosomal_recessive': 'Autosomal Recessive (隱性)',
                        'x_linked': 'X-linked (性聯)',
                        'de_novo': 'De novo (新發)',
                    }
                    inheritance_display = pattern_names.get(pattern, pattern.replace('_', ' ').title())
                    if confidence:
                        inheritance_display += f' <span style="font-size: 0.85em; opacity: 0.8;">[{confidence} confidence]</span>'
            
            # Format age onset and distribution
            age_onset_display = article.get('age_onset', 'Not specified')
            if clinical_info.get('age_distribution', {}).get('age_onset'):
                age_onset = clinical_info['age_distribution']['age_onset']
                if age_onset != 'not_specified':
                    age_names = {
                        'congenital': 'Congenital (先天)',
                        'infantile': 'Infantile (嬰兒)',
                        'childhood': 'Childhood (兒童)',
                        'adolescent': 'Adolescent (青少年)',
                        'adult': 'Adult (成人)',
                        'late-onset': 'Late-onset (晚發)',
                    }
                    age_onset_display = age_names.get(age_onset, age_onset.title())
                    # Add mean age if available
                    mean_age = clinical_info['age_distribution'].get('mean_age')
                    median_age = clinical_info['age_distribution'].get('median_age')
                    age_ranges = clinical_info['age_distribution'].get('age_ranges', [])
                    
                    age_details = []
                    if mean_age:
                        age_details.append(f'平均 {mean_age} 歲')
                    if median_age:
                        age_details.append(f'中位數 {median_age} 歲')
                    if age_ranges:
                        age_details.append(f'範圍 {", ".join(age_ranges)}')
                    
                    if age_details:
                        age_onset_display += f' <span style="font-size: 0.85em; opacity: 0.8;">({", ".join(age_details)})</span>'
            
            # Format affected tissue
            affected_tissue_display = 'Not specified'
            if clinical_info.get('affected_tissue'):
                tissue = clinical_info['affected_tissue']
                cardiac = tissue.get('cardiac', False)
                skeletal = tissue.get('skeletal_muscle', False)
                
                tissue_parts = []
                if cardiac:
                    phenotypes = tissue.get('cardiac_phenotype', [])
                    if phenotypes:
                        tissue_parts.append(f'Cardiac (心肌): {", ".join(phenotypes)}')
                    else:
                        tissue_parts.append('Cardiac (心肌)')
                
                if skeletal:
                    phenotypes = tissue.get('skeletal_phenotype', [])
                    if phenotypes:
                        tissue_parts.append(f'Skeletal muscle (骨骼肌): {", ".join(phenotypes)}')
                    else:
                        tissue_parts.append('Skeletal muscle (骨骼肌)')
                
                if tissue_parts:
                    affected_tissue_display = '<br>'.join(tissue_parts)
                elif cardiac or skeletal:
                    affected_tissue_display = 'Specified but details unclear'
            
            # Additional clinical info row
            extra_info_html = ""
            sample_size = clinical_info.get('sample_size', {})
            n_patients = sample_size.get('patients')
            n_families = sample_size.get('families')
            severity = clinical_info.get('severity')
            
            if n_patients or n_families or (severity and severity != 'not_specified'):
                extra_items = []
                if n_patients:
                    extra_items.append(f'<div class="article-info-item"><div class="article-info-label">Sample Size:</div><div class="article-info-value">{n_patients} patients' + (f', {n_families} families' if n_families else '') + '</div></div>')
                if severity and severity != 'not_specified':
                    severity_badge = f'<span class="badge badge-{severity}">{severity.title()}</span>'
                    extra_items.append(f'<div class="article-info-item"><div class="article-info-label">Severity:</div><div class="article-info-value">{severity_badge}</div></div>')
                
                if extra_items:
                    extra_info_html = '<div class="article-info" style="margin-top: 10px;">' + ''.join(extra_items) + '</div>'
            
            # Key findings if available
            key_findings_html = ""
            key_findings = clinical_info.get('key_findings')
            if key_findings:
                key_findings_html = f'''
                <div style="margin-top: 15px; padding: 12px; background: #f0f7ff; border-left: 4px solid #667eea; border-radius: 4px;">
                    <strong style="color: #667eea;">Key Findings:</strong>
                    <div style="margin-top: 5px; color: #555;">{key_findings}</div>
                </div>
                '''
            
            articles_html += f"""
            <div class="pubmed-article">
                <div class="article-title">
                    {idx}. <a href="{article['pubmed_link']}" target="_blank">{article['title']}</a>
                </div>
                <div class="article-meta">
                    {article['authors']} | {article['journal']} ({article['year']}) | 
                    PMID: <a href="{article['pubmed_link']}" target="_blank">{article['pmid']}</a>
                </div>
                
                <div class="article-info">
                    <div class="article-info-item">
                        <div class="article-info-label">Age Distribution (年齡分布):</div>
                        <div class="article-info-value">{age_onset_display}</div>
                    </div>
                    <div class="article-info-item">
                        <div class="article-info-label">Affected Tissue (影響組織):</div>
                        <div class="article-info-value">{affected_tissue_display}</div>
                    </div>
                    <div class="article-info-item">
                        <div class="article-info-label">Inheritance (遺傳模式):</div>
                        <div class="article-info-value">{inheritance_display}</div>
                    </div>
                </div>
                
                {extra_info_html}
                {key_findings_html}
            </div>
            """
        
        return f"""
        <section class="section">
            <h2>Literature Review (PubMed)</h2>
            <p>Found {len(pubmed_results)} relevant articles in PubMed:</p>
            {age_summary}
            {articles_html}
        </section>
        """
    
    def _generate_age_summary(self, pubmed_results: List[Dict]) -> str:
        """Generate age distribution summary"""
        age_categories = {}
        for article in pubmed_results:
            age_onset = article.get('age_onset', 'Not specified')
            if age_onset != 'Not specified':
                age_categories[age_onset] = age_categories.get(age_onset, 0) + 1
        
        if not age_categories:
            return ""
        
        categories_html = ""
        for age, count in sorted(age_categories.items(), key=lambda x: x[1], reverse=True):
            categories_html += f"""
                <div style="padding: 8px; background: #f8f9fa; border-radius: 4px; margin: 5px 0;">
                    <strong>{age}:</strong> {count} article(s)
                </div>
            """
        
        return f"""
        <div class="alert alert-info" style="margin: 15px 0;">
            <strong>年齡分布統計 (Age Distribution Summary):</strong>
            <div style="margin-top: 10px;">
                {categories_html}
            </div>
        </div>
        """
    
    def _generate_clinvar_section(self, clinvar_info: Dict, pubmed_results: List[Dict] = None) -> str:
        """Generate ClinVar information section"""
        if not clinvar_info:
            return ""
        
        # Format conditions
        conditions_html = ""
        if clinvar_info.get('conditions'):
            for condition in clinvar_info['conditions']:
                conditions_html += f'<li>{condition}</li>'
        else:
            conditions_html = '<li>Not specified</li>'
        
        # Format variant impact
        impact_badge = self._get_impact_badge(clinvar_info.get('variant_impact', 'Not specified'))
        
        # Format clinical significance
        significance = clinvar_info.get('clinical_significance', 'Not specified')
        significance_class = self._get_significance_class(significance)
        
        # PubMed IDs from ClinVar - 使用實際獲取的文章 PMID（與 Literature Review 一致）
        pmid_info = ""
        if pubmed_results:
            # 從 pubmed_results 中提取 PMID（這些是實際獲取到的文章）
            actual_pmids = [article['pmid'] for article in pubmed_results]
            pmid_count = len(actual_pmids)
            pmids_formatted = ', '.join([f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">{pmid}</a>' 
                                        for pmid in actual_pmids[:10]])
            if pmid_count > 10:
                pmids_formatted += f' ... and {pmid_count - 10} more'
            pmid_info = f"""
            <div class="info-row">
                <div class="info-label">PubMed IDs (已獲取的文章):</div>
                <div class="info-value">{pmid_count} articles ({pmids_formatted})</div>
            </div>
            """
        
        return f"""
        <section class="section">
            <h2>ClinVar Information</h2>
            
            <div class="info-grid" style="margin-top: 20px;">
                <div class="info-row">
                    <div class="info-label">臨床意義 (Clinical Significance):</div>
                    <div class="info-value">
                        <span class="badge badge-{significance_class}">{significance}</span>
                    </div>
                </div>
                
                <div class="info-row">
                    <div class="info-label">審查狀態 (Review Status):</div>
                    <div class="info-value">{clinvar_info.get('review_status', 'Not specified')}</div>
                </div>
                
                <div class="info-row">
                    <div class="info-label">變異影響 (Variant Impact):</div>
                    <div class="info-value">{impact_badge}</div>
                </div>
                
                {pmid_info}
            </div>
            
            <div style="margin-top: 20px;">
                <h3 style="font-size: 1.1em; margin-bottom: 10px;">相關疾病 (Associated Conditions):</h3>
                <ul style="margin-left: 20px; line-height: 1.8;">
                    {conditions_html}
                </ul>
            </div>
        </section>
        """
    
    def _get_impact_badge(self, impact: str) -> str:
        """Get HTML badge for variant impact"""
        impact_lower = impact.lower()
        if 'cardiac' in impact_lower and 'skeletal' in impact_lower:
            return '<span class="badge badge-both">Cardiac & Skeletal</span>'
        elif 'cardiac' in impact_lower:
            return '<span class="badge badge-cardiac">Cardiac</span>'
        elif 'skeletal' in impact_lower:
            return '<span class="badge badge-skeletal">Skeletal Muscle</span>'
        else:
            return f'<span class="badge">{impact}</span>'
    
    def _get_significance_class(self, significance: str) -> str:
        """Get CSS class for clinical significance"""
        sig_lower = significance.lower()
        if 'pathogenic' in sig_lower and 'likely' not in sig_lower:
            return 'pathogenic'
        elif 'likely pathogenic' in sig_lower:
            return 'likely-pathogenic'
        elif 'benign' in sig_lower and 'likely' not in sig_lower:
            return 'benign'
        elif 'likely benign' in sig_lower:
            return 'likely-benign'
        elif 'uncertain' in sig_lower:
            return 'uncertain'
        else:
            return 'default'
    
    def _get_phenotype_badge(self, phenotype: str) -> str:
        """Get HTML badge for phenotype"""
        phenotype_lower = phenotype.lower()
        if 'cardiac' in phenotype_lower or 'heart' in phenotype_lower:
            return '<span class="badge badge-cardiac">Cardiac</span>'
        elif 'skeletal' in phenotype_lower:
            return '<span class="badge badge-skeletal">Skeletal Muscle</span>'
        elif 'both' in phenotype_lower:
            return '<span class="badge badge-both">Cardiac & Skeletal</span>'
        else:
            return phenotype
    
    def _generate_clinical_stats_section(self, aggregated_stats: Dict) -> str:
        """Generate clinical statistics summary section from LLM extraction"""
        if not aggregated_stats:
            return ""
        
        # Format inheritance patterns
        inheritance_html = ""
        inheritance_dist = aggregated_stats.get('inheritance_distribution', {})
        if inheritance_dist:
            inheritance_items = []
            pattern_names = {
                'autosomal_dominant': 'Autosomal Dominant (顯性遺傳)',
                'autosomal_recessive': 'Autosomal Recessive (隱性遺傳)',
                'x_linked': 'X-linked (性聯遺傳)',
                'de_novo': 'De novo (新發突變)',
                'compound_heterozygous': 'Compound Heterozygous (複合雜合)',
            }
            for pattern, count in sorted(inheritance_dist.items(), key=lambda x: x[1], reverse=True):
                if count > 0:
                    display_name = pattern_names.get(pattern, pattern.replace('_', ' ').title())
                    inheritance_items.append(f"<li>{display_name}: {count} 篇文章</li>")
            if inheritance_items:
                inheritance_html = f"<ul style='margin: 10px 0; padding-left: 25px;'>{''.join(inheritance_items)}</ul>"
        
        # Format age onset distribution
        age_onset_html = ""
        age_onset_dist = aggregated_stats.get('age_onset_distribution', {})
        if age_onset_dist:
            age_items = []
            age_names = {
                'congenital': 'Congenital (先天)',
                'infantile': 'Infantile (嬰兒期)',
                'childhood': 'Childhood (兒童期)',
                'adolescent': 'Adolescent (青少年期)',
                'adult': 'Adult (成人)',
                'late-onset': 'Late-onset (晚發)',
            }
            for age, count in sorted(age_onset_dist.items(), key=lambda x: x[1], reverse=True):
                if count > 0:
                    display_name = age_names.get(age, age.replace('_', ' ').title())
                    age_items.append(f"<li>{display_name}: {count} 篇文章</li>")
            if age_items:
                age_onset_html = f"<ul style='margin: 10px 0; padding-left: 25px;'>{''.join(age_items)}</ul>"
        
        # Format tissue affected
        tissue_html = ""
        tissue_affected = aggregated_stats.get('tissue_affected', {})
        if tissue_affected:
            cardiac_count = tissue_affected.get('cardiac', 0)
            skeletal_count = tissue_affected.get('skeletal_muscle', 0)
            both_count = tissue_affected.get('both', 0)
            tissue_html = f"""
                <ul style='margin: 10px 0; padding-left: 25px;'>
                    <li>心臟 (Cardiac): {cardiac_count} 篇文章</li>
                    <li>骨骼肌 (Skeletal Muscle): {skeletal_count} 篇文章</li>
                    <li>兩者皆有 (Both): {both_count} 篇文章</li>
                </ul>
            """
        
        # Format phenotypes
        cardiac_phenotypes_html = ""
        cardiac_phenotypes = aggregated_stats.get('cardiac_phenotypes', {})
        if cardiac_phenotypes:
            phenotype_items = []
            for phenotype, count in sorted(cardiac_phenotypes.items(), key=lambda x: x[1], reverse=True):
                if count > 0:
                    phenotype_items.append(f"<li>{phenotype}: {count} 篇</li>")
            if phenotype_items:
                cardiac_phenotypes_html = f"<ul style='margin: 10px 0; padding-left: 25px;'>{''.join(phenotype_items)}</ul>"
        
        skeletal_phenotypes_html = ""
        skeletal_phenotypes = aggregated_stats.get('skeletal_phenotypes', {})
        if skeletal_phenotypes:
            phenotype_items = []
            for phenotype, count in sorted(skeletal_phenotypes.items(), key=lambda x: x[1], reverse=True):
                if count > 0:
                    phenotype_items.append(f"<li>{phenotype}: {count} 篇</li>")
            if phenotype_items:
                skeletal_phenotypes_html = f"<ul style='margin: 10px 0; padding-left: 25px;'>{''.join(phenotype_items)}</ul>"
        
        return f"""
        <section class="section" style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 10px; margin: 20px;">
            <h2 style="color: white; border-bottom: 3px solid rgba(255,255,255,0.3);">臨床資訊統計 (Clinical Information Summary)</h2>
            <p style="opacity: 0.9; margin-bottom: 20px;">基於 {aggregated_stats.get('total_articles', 0)} 篇文獻的 AI 分析結果</p>
            
            <div class="info-grid" style="display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 20px; margin-top: 20px;">
                
                <div style="background: rgba(255,255,255,0.15); padding: 20px; border-radius: 10px; backdrop-filter: blur(10px);">
                    <h3 style="color: white; margin-top: 0; font-size: 1.2em;">👥 研究規模</h3>
                    <div style="font-size: 1.1em; line-height: 1.8;">
                        <div><strong>總患者數:</strong> {aggregated_stats.get('total_patients', 0)}</div>
                        <div><strong>總家族數:</strong> {aggregated_stats.get('total_families', 0)}</div>
                    </div>
                </div>
                
                {'<div style="background: rgba(255,255,255,0.15); padding: 20px; border-radius: 10px; backdrop-filter: blur(10px);"><h3 style="color: white; margin-top: 0; font-size: 1.2em;">遺傳模式</h3>' + inheritance_html + '</div>' if inheritance_html else ''}
                
                {'<div style="background: rgba(255,255,255,0.15); padding: 20px; border-radius: 10px; backdrop-filter: blur(10px);"><h3 style="color: white; margin-top: 0; font-size: 1.2em;">發病年齡</h3>' + age_onset_html + '</div>' if age_onset_html else ''}
                
                {'<div style="background: rgba(255,255,255,0.15); padding: 20px; border-radius: 10px; backdrop-filter: blur(10px);"><h3 style="color: white; margin-top: 0; font-size: 1.2em;">受影響組織</h3>' + tissue_html + '</div>' if tissue_html else ''}
                
                {'<div style="background: rgba(255,255,255,0.15); padding: 20px; border-radius: 10px; backdrop-filter: blur(10px);"><h3 style="color: white; margin-top: 0; font-size: 1.2em;">心臟表型</h3>' + cardiac_phenotypes_html + '</div>' if cardiac_phenotypes_html else ''}
                
                {'<div style="background: rgba(255,255,255,0.15); padding: 20px; border-radius: 10px; backdrop-filter: blur(10px);"><h3 style="color: white; margin-top: 0; font-size: 1.2em;">骨骼肌表型</h3>' + skeletal_phenotypes_html + '</div>' if skeletal_phenotypes_html else ''}
                
            </div>
            
            <div style="margin-top: 20px; padding: 15px; background: rgba(255,255,255,0.1); border-radius: 8px; font-size: 0.9em;">
                <strong>說明：</strong> 以上統計資料由 AI 模型自動從文獻中提取，僅供參考。詳細資訊請參閱下方各篇文獻的完整內容。
            </div>
        </section>
        """
    
    def _generate_footer(self) -> str:
        """Generate footer section"""
        return f"""
        <footer>
            <p><strong>Disclaimer:</strong> This report is for research purposes only and should not be used 
            for clinical decision-making without proper validation.</p>
            <p>Generated by TTN Variant AI Agent | Powered by Evo2, PubMed, and Python</p>
            <p>© {datetime.now().year} | Version 1.0</p>
        </footer>
        """
    
    def _get_javascript(self) -> str:
        """Generate JavaScript code"""
        return """
        // Add smooth scrolling
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                document.querySelector(this.getAttribute('href')).scrollIntoView({
                    behavior: 'smooth'
                });
            });
        });
        
        // Print functionality
        function printReport() {
            window.print();
        }
        
        console.log('TTN Variant AI Report loaded successfully');
        """