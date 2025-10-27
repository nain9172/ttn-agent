#!/usr/bin/env python3
"""
TTN Variant AI Agent - Main Entry Point
Orchestrates the complete pipeline for variant analysis
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime

from config import OUTPUT_DIR, LOG_FILE
from utils.variant_parser import parse_variant
from utils.evo2_predictor import Evo2Predictor
from utils.pubmed_search import PubMedSearcher
from utils.image_generator import ImageGenerator
from utils.html_report import HTMLReportGenerator


def setup_logging():
    """Setup logging configuration"""
    import logging
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(LOG_FILE),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def main():
    """Main pipeline orchestrator"""
    logger = setup_logging()
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='TTN Variant AI Agent - Comprehensive variant analysis'
    )
    parser.add_argument(
        'variant',
        type=str,
        help='Variant in format: chromosome-position-ref-alt (e.g., 2-178612477-T-A)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Output HTML file path (default: outputs/variant_report_TIMESTAMP.html)'
    )
    parser.add_argument(
        '--skip-evo2',
        action='store_true',
        help='Skip Evo2 prediction (faster but less accurate)'
    )
    
    args = parser.parse_args()
    
    logger.info(f"Starting TTN Variant AI Agent for variant: {args.variant}")
    
    try:
        # Step 1: Parse variant
        logger.info("Step 1: Parsing variant...")
        variant_info = parse_variant(args.variant)
        logger.info(f"Parsed variant: {variant_info}")
        
        # Step 2: Evo2 Prediction
        evo2_result = None
        if not args.skip_evo2:
            logger.info("Step 2: Running Evo2 prediction...")
            predictor = Evo2Predictor()
            evo2_result = predictor.predict(variant_info)
            logger.info(f"Evo2 prediction complete: {evo2_result}")
        else:
            logger.info("Step 2: Skipped Evo2 prediction")
        
        # Step 3: PubMed Search
        logger.info("Step 3: Searching PubMed...")
        pubmed_searcher = PubMedSearcher()
        pubmed_results = pubmed_searcher.search(variant_info)
        logger.info(f"Found {len(pubmed_results)} PubMed articles")
        
        # Step 4: Generate Protein Schematic
        logger.info("Step 4: Generating protein schematic...")
        image_gen = ImageGenerator()
        image_path = image_gen.generate_titin_schematic(variant_info)
        logger.info(f"Generated image: {image_path}")
        
        # Step 5: Generate HTML Report
        logger.info("Step 5: Generating HTML report...")
        report_gen = HTMLReportGenerator()
        
        # Determine output path
        if args.output:
            output_path = Path(args.output)
        else:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = OUTPUT_DIR / f"variant_report_{timestamp}.html"
        
        report_path = report_gen.generate_report(
            variant_info=variant_info,
            evo2_result=evo2_result,
            pubmed_results=pubmed_results,
            image_path=image_path,
            output_path=output_path
        )
        
        logger.info(f"✅ Report generated successfully: {report_path}")
        print(f"\n{'='*70}")
        print(f"🎉 SUCCESS! Report generated at:")
        print(f"   {report_path}")
        print(f"{'='*70}\n")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error in pipeline: {str(e)}", exc_info=True)
        print(f"\n Error: {str(e)}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())