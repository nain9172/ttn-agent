#!/usr/bin/env python3
"""
Batch Variant Analysis Script
Process multiple variants from a file
"""

import argparse
import sys
import csv
from pathlib import Path
from datetime import datetime
import logging

from main import setup_logging
from utils.variant_parser import parse_variant
from utils.evo2_predictor import Evo2Predictor
from utils.clinvar_parser import ClinVarParser
from utils.pubmed_search import PubMedSearcher
from utils.image_generator import ImageGenerator
from utils.html_report import HTMLReportGenerator
from config import OUTPUT_DIR


def read_variants_from_file(file_path: Path) -> list:
    """
    Read variants from CSV file
    Expected format: chromosome,position,ref,alt
    Or: variant_id (format: chr-pos-ref-alt)
    """
    variants = []
    
    with open(file_path, 'r') as f:
        # Try to detect format
        first_line = f.readline().strip()
        f.seek(0)
        
        if ',' in first_line:
            # CSV format
            reader = csv.DictReader(f)
            for row in reader:
                if 'variant_id' in row:
                    variants.append(row['variant_id'])
                else:
                    variant = f"{row['chromosome']}-{row['position']}-{row['ref']}-{row['alt']}"
                    variants.append(variant)
        else:
            # Simple list format (one variant per line)
            variants = [line.strip() for line in f if line.strip()]
    
    return variants


def process_batch(
    variants: list,
    output_dir: Path,
    skip_evo2: bool = False
) -> dict:
    """
    Process a batch of variants
    
    Returns:
        Dictionary with summary statistics
    """
    logger = logging.getLogger(__name__)
    
    # Initialize components (reuse for efficiency)
    clinvar_parser = ClinVarParser()
    evo2_predictor = None if skip_evo2 else Evo2Predictor()
    pubmed_searcher = PubMedSearcher()
    image_generator = ImageGenerator()
    report_generator = HTMLReportGenerator()
    
    results = {
        'total': len(variants),
        'successful': 0,
        'failed': 0,
        'pathogenic': 0,
        'benign': 0,
        'errors': []
    }
    
    logger.info(f"Processing {len(variants)} variants...")
    
    for idx, variant_str in enumerate(variants, 1):
        logger.info(f"\n{'='*70}")
        logger.info(f"Processing variant {idx}/{len(variants)}: {variant_str}")
        logger.info(f"{'='*70}")
        
        try:
            # Parse variant
            variant_info = parse_variant(variant_str)
            
            # Parse ClinVar information
            logger.info("Parsing ClinVar information...")
            clinvar_info = clinvar_parser.parse_variant(variant_info)
            if clinvar_info and clinvar_info.get('pmid_list'):
                logger.info(f"Found {len(clinvar_info['pmid_list'])} PubMed IDs from ClinVar")
            else:
                clinvar_info = None
            
            # Evo2 prediction
            evo2_result = None
            if not skip_evo2:
                logger.info("Running Evo2 prediction...")
                evo2_result = evo2_predictor.predict(variant_info)
                
                if evo2_result.get('success'):
                    if evo2_result['prediction'] == 'pathogenic':
                        results['pathogenic'] += 1
                    else:
                        results['benign'] += 1
            
            # PubMed search (只使用 ClinVar PMIDs，不進行廣泛搜索)
            logger.info("Searching PubMed...")
            pubmed_results = []
            
            if clinvar_info and clinvar_info.get('pmid_list'):
                pmid_list = clinvar_info.get('pmid_list')
                logger.info(f"使用 ClinVar 提供的 {len(pmid_list)} 個精確 PubMed IDs")
                pubmed_results = pubmed_searcher.search(variant_info, pmid_list=pmid_list)
                logger.info(f"成功獲取 {len(pubmed_results)} 篇 ClinVar 相關文獻")
            else:
                logger.warning(f"⚠️  變異 {variant_str}: ClinVar 未提供 PubMed IDs，跳過文獻搜索")
            
            # Generate image
            logger.info("Generating protein schematic...")
            image_path = image_generator.generate_titin_schematic(variant_info)
            
            # Generate report
            logger.info("Generating HTML report...")
            report_path = output_dir / f"report_{variant_info['variant_id']}.html"
            report_generator.generate_report(
                variant_info=variant_info,
                evo2_result=evo2_result,
                pubmed_results=pubmed_results,
                image_path=image_path,
                output_path=report_path,
                clinvar_info=clinvar_info
            )
            
            results['successful'] += 1
            logger.info(f"✅ Successfully processed: {variant_str}")
            
        except Exception as e:
            results['failed'] += 1
            error_msg = f"{variant_str}: {str(e)}"
            results['errors'].append(error_msg)
            logger.error(f"Error processing {variant_str}: {e}", exc_info=True)
    
    return results


def generate_summary_report(results: dict, output_path: Path):
    """Generate summary CSV report"""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Metric', 'Value'])
        writer.writerow(['Total Variants', results['total']])
        writer.writerow(['Successful', results['successful']])
        writer.writerow(['Failed', results['failed']])
        writer.writerow(['Pathogenic', results['pathogenic']])
        writer.writerow(['Benign', results['benign']])
        writer.writerow([''])
        writer.writerow(['Errors'])
        for error in results['errors']:
            writer.writerow(['', error])


def main():
    """Main batch processing function"""
    parser = argparse.ArgumentParser(
        description='Batch process multiple TTN variants'
    )
    parser.add_argument(
        'input_file',
        type=str,
        help='Input file with variants (CSV or text list)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory for reports'
    )
    parser.add_argument(
        '--skip-evo2',
        action='store_true',
        help='Skip Evo2 prediction for faster processing'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging()
    
    # Read variants
    input_path = Path(args.input_file)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        return 1
    
    logger.info(f"Reading variants from: {input_path}")
    variants = read_variants_from_file(input_path)
    logger.info(f"Found {len(variants)} variants to process")
    
    if not variants:
        logger.error("No variants found in input file")
        return 1
    
    # Setup output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = OUTPUT_DIR / f"batch_{timestamp}"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Process batch
    logger.info("\nStarting batch processing...")
    results = process_batch(
        variants=variants,
        output_dir=output_dir,
        skip_evo2=args.skip_evo2
    )
    
    # Generate summary
    summary_path = output_dir / "summary.csv"
    generate_summary_report(results, summary_path)
    
    # Print results
    print("\n" + "="*70)
    print("BATCH PROCESSING SUMMARY")
    print("="*70)
    print(f"Total variants:     {results['total']}")
    print(f"Successful:         {results['successful']}")
    print(f"Failed:             {results['failed']}")
    if not args.skip_evo2:
        print(f"Pathogenic:         {results['pathogenic']}")
        print(f"Benign:             {results['benign']}")
    print(f"\nOutput directory:   {output_dir}")
    print(f"Summary report:     {summary_path}")
    print("="*70)
    
    if results['errors']:
        print("\nErrors encountered:")
        for error in results['errors']:
            print(f"  - {error}")
    
    return 0 if results['failed'] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())