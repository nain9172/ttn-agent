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
        transcript_intervals_path: Optional[Path] = None
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
            image_base64
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
        image_base64: str
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
        {self._generate_pubmed_section(pubmed_results)}
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
            <h2>📊 Variant Summary</h2>
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
                <h2>🧬 Evo2 Pathogenicity Prediction</h2>
                <div class="alert">Evo2 prediction was skipped.</div>
            </section>
            """
        
        if not evo2_result.get('success'):
            error_msg = evo2_result.get('error', 'Unknown error')
            return f"""
            <section class="section">
                <h2>🧬 Evo2 Pathogenicity Prediction</h2>
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
            <h2>🧬 Evo2 Pathogenicity Prediction</h2>
            
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
                <h2>🔬 Protein Domain Localization</h2>
                <div class="alert">Image generation failed.</div>
            </section>
            """
        
        return f"""
        <section class="section">
            <h2>🔬 Protein Domain Localization & Transcript Intervals</h2>
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
            <h2>🧬 Transcript Intervals (Genomic Coordinates)</h2>
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
    
    def _generate_pubmed_section(self, pubmed_results: List[Dict]) -> str:
        """Generate PubMed literature section"""
        if not pubmed_results:
            return """
            <section class="section">
                <h2>📚 Literature Review (PubMed)</h2>
                <div class="alert">No relevant articles found in PubMed.</div>
            </section>
            """
        
        articles_html = ""
        for idx, article in enumerate(pubmed_results, 1):
            phenotype_badge = self._get_phenotype_badge(article.get('phenotype', 'Not specified'))
            
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
                        <div class="article-info-label">Phenotype:</div>
                        <div class="article-info-value">{phenotype_badge}</div>
                    </div>
                    <div class="article-info-item">
                        <div class="article-info-label">Inheritance:</div>
                        <div class="article-info-value">{article.get('inheritance', 'Not specified')}</div>
                    </div>
                    <div class="article-info-item">
                        <div class="article-info-label">Age of Onset:</div>
                        <div class="article-info-value">{article.get('age_onset', 'Not specified')}</div>
                    </div>
                </div>
                
                <p style="margin-top: 10px; color: #666; font-size: 0.9em;">
                    {article.get('abstract', '')}
                </p>
            </div>
            """
        
        return f"""
        <section class="section">
            <h2>📚 Literature Review (PubMed)</h2>
            <p>Found {len(pubmed_results)} relevant articles in PubMed:</p>
            {articles_html}
        </section>
        """
    
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