#!/usr/bin/env python3
"""
TTN Variant AI Agent - Main Entry Point (Improved)
Fixes for variant-specific extraction
"""

import argparse
import sys
import os
from pathlib import Path
from datetime import datetime

from config import (
    OUTPUT_DIR, 
    LOG_FILE, 
    PROJECT_ROOT,
    ENABLE_LOCAL_CLINICAL_EXTRACTION,
    LOCAL_LLM_BACKEND,
    LOCAL_LLM_MODEL,
    LOCAL_LLM_TENSOR_PARALLEL,
    LOCAL_LLM_MAX_MODEL_LEN,
    LOCAL_LLM_MAX_CONTEXT_LENGTH,
    ENABLE_LITVAR_SEARCH,
    TRITON_PTXAS_PATH
)

# Try to get ENABLE_EASY_PROMPT setting
try:
    from config import ENABLE_EASY_PROMPT
except ImportError:
    ENABLE_EASY_PROMPT = False

# 設定 Triton PTXAS 環境變數（用於 vLLM 編譯）
os.environ['TRITON_PTXAS_PATH'] = TRITON_PTXAS_PATH

# Try to get enhanced PubMed settings
try:
    from config import ENABLE_FULL_TEXT_FETCH, MAX_TEXT_LENGTH
except ImportError:
    ENABLE_FULL_TEXT_FETCH = False
    MAX_TEXT_LENGTH = 8000

# Try to get Docling settings
try:
    from config import ENABLE_DOCLING_PDF
except ImportError:
    ENABLE_DOCLING_PDF = False

from utils.variant_parser import parse_variant
from utils.evo2_predictor import Evo2Predictor
from utils.clinvar_parser import ClinVarParser
from utils.litvar_search import LitVarSearcher
from utils.variant_utils import get_variant_aliases
# Try to import enhanced PubMed searcher, fallback to standard
try:
    from utils.enhanced_pubmed_search import EnhancedPubMedSearcher
    PUBMED_SEARCHER_CLASS = EnhancedPubMedSearcher
    ENHANCED_AVAILABLE = True
except ImportError:
    from utils.pubmed_search import PubMedSearcher
    PUBMED_SEARCHER_CLASS = PubMedSearcher
    ENHANCED_AVAILABLE = False

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
    parser.add_argument(
        '--skip-clinical-extraction',
        action='store_true',
        help='Skip local LLM clinical information extraction'
    )
    parser.add_argument(
        '--llm-backend',
        type=str,
        default=None,
        choices=['ollama', 'transformers', 'vllm'],
        help='Override LLM backend (default: from config.py)'
    )
    parser.add_argument(
        '--llm-model',
        type=str,
        default=None,
        help='Override LLM model (default: from config.py)'
    )
    
    args = parser.parse_args()
    
    logger.info(f"Starting TTN Variant AI Agent for variant: {args.variant}")
    
    try:
        # Step 1: Parse variant
        logger.info("Step 1: Parsing variant...")
        variant_info = parse_variant(args.variant)
        logger.info(f"Parsed variant: {variant_info}")
        
        # Step 2: Parse ClinVar Information
        logger.info("Step 2: Parsing ClinVar information...")
        clinvar_parser = ClinVarParser()
        clinvar_info = clinvar_parser.parse_variant(variant_info)
        if clinvar_info and clinvar_info.get('pmid_list'):
            logger.info(f"Found {len(clinvar_info['pmid_list'])} PubMed IDs from ClinVar")
            logger.info(f"Variant impact: {clinvar_info.get('variant_impact')}")
            logger.info(f"Clinical significance: {clinvar_info.get('clinical_significance')}")
        else:
            logger.warning("No ClinVar information found (this is okay, will use standard PubMed search)")
            clinvar_info = None
        
        # Step 3: Evo2 Prediction
        evo2_result = None
        if not args.skip_evo2:
            logger.info("Step 3: Running Evo2 prediction...")
            predictor = Evo2Predictor()
            evo2_result = predictor.predict(variant_info)
            logger.info(f"Evo2 prediction complete: {evo2_result}")
            
            # 積極釋放 Evo2 模型佔用的 GPU 記憶體
            logger.info("Clearing Evo2 from GPU memory...")
            del predictor
            import torch
            import gc
            
            # 清理 Python 垃圾回收
            gc.collect()
            
            # 清理 PyTorch CUDA 快取
            torch.cuda.empty_cache()
            
            # 重置 CUDA 峰值記憶體統計
            torch.cuda.reset_peak_memory_stats()
            
            # 顯示清理後的記憶體狀態
            for i in range(torch.cuda.device_count()):
                mem_allocated = torch.cuda.memory_allocated(i) / 1024**3
                mem_reserved = torch.cuda.memory_reserved(i) / 1024**3
                logger.info(f"GPU {i} after cleanup: {mem_allocated:.2f}GB allocated, {mem_reserved:.2f}GB reserved")
            
            logger.info("✓ Cleared Evo2 from GPU memory")
            
            # 等待一下確保清理完成
            import time
            time.sleep(2)
        else:
            logger.info("Step 3: Skipped Evo2 prediction")
        
        # Step 4: PubMed Search
        logger.info("Step 4: Searching PubMed...")
        
        # 從 ClinVar 信息生成變異別名（用於表格過濾）
        variant_id = f"{variant_info['chrom']}-{variant_info['pos']}-{variant_info['ref']}-{variant_info['alt']}"
        variant_aliases = get_variant_aliases(variant_id, clinvar_info)
        logger.info(f"從 ClinVar 獲取 {len(variant_aliases)} 個變異別名用於表格過濾")
        
        # Initialize PubMed searcher (enhanced or standard)
        # 如果啟用 EASY_PROMPT，跳過全文下載和 docling 處理
        if ENABLE_EASY_PROMPT:
            logger.info("使用標準 PubMed 搜尋器（EASY_PROMPT 模式：僅摘要，跳過全文下載和 Docling 處理）")
            if ENHANCED_AVAILABLE:
                pubmed_searcher = PUBMED_SEARCHER_CLASS(try_full_text=False, use_docling=False)
            else:
                pubmed_searcher = PUBMED_SEARCHER_CLASS()
        elif ENHANCED_AVAILABLE and ENABLE_FULL_TEXT_FETCH:
            use_docling = ENABLE_DOCLING_PDF
            logger.info("使用增強版 PubMed 搜尋器（支持全文獲取）")
            if use_docling:
                logger.info("  ✓ Docling PDF 處理已啟用（優先提取 Tables/Results/Supplementary）")
            pubmed_searcher = PUBMED_SEARCHER_CLASS(
                try_full_text=True,
                max_text_length=MAX_TEXT_LENGTH,
                use_docling=use_docling
            )
        elif ENHANCED_AVAILABLE:
            logger.info("使用增強版 PubMed 搜尋器（全文獲取已禁用）")
            pubmed_searcher = PUBMED_SEARCHER_CLASS(try_full_text=False, use_docling=False)
        else:
            logger.info("使用標準 PubMed 搜尋器（僅摘要）")
            pubmed_searcher = PUBMED_SEARCHER_CLASS()
        
        pubmed_results = []
        
        # 只使用 ClinVar 提供的 PMIDs，不進行廣泛搜索
        if clinvar_info and clinvar_info.get('pmid_list'):
            pmid_list = clinvar_info.get('pmid_list')
            logger.info(f"使用 ClinVar 提供的 {len(pmid_list)} 個精確 PubMed IDs")
            pubmed_results = pubmed_searcher.search(variant_info, pmid_list=pmid_list, variant_aliases=variant_aliases)
            logger.info(f"成功獲取 {len(pubmed_results)} 篇 ClinVar 相關文獻")
        else:
            logger.warning("ClinVar 未提供 PubMed IDs，報告將不包含文獻（避免不相關文章）")
            logger.info("如需查看文獻，請確認該變異在 ClinVar 資料庫中有記錄")
        
        # Step 4a: LitVar Search (BEFORE clinical extraction so articles are processed together)
        # 如果啟用 EASY_PROMPT，跳過 LitVar 搜尋
        if ENABLE_EASY_PROMPT:
            logger.info("Step 4a: Skipped LitVar search (EASY_PROMPT mode)")
        elif ENABLE_LITVAR_SEARCH:
            logger.info("Step 4a: Searching LitVar using rsID...")
            # 在 EASY_PROMPT 模式下不使用 docling
            use_docling_for_litvar = False if ENABLE_EASY_PROMPT else ENABLE_DOCLING_PDF
            litvar_searcher = LitVarSearcher(use_docling=use_docling_for_litvar)
            
            # LitVar uses rsID to search - this is extracted from clinvar_info
            litvar_results = litvar_searcher.search_multiple_formats(
                variant_info=variant_info,
                clinvar_info=clinvar_info,  # Contains rsid from ClinVar
                variant_aliases=variant_aliases,  # 傳遞別名用於表格過濾
                max_results=20
            )
            
            logger.info(f"Found {len(litvar_results)} articles from LitVar")
            
            # Merge results (remove duplicates by PMID)
            pubmed_pmids = {a['pmid'] for a in pubmed_results}
            new_articles = [a for a in litvar_results if a['pmid'] not in pubmed_pmids]
            
            logger.info(f"Adding {len(new_articles)} new articles from LitVar")
            pubmed_results.extend(new_articles)
            
            logger.info(f"Total articles after merging ClinVar + LitVar: {len(pubmed_results)}")
        else:
            logger.info("Step 4a: LitVar search disabled")
        
        # Step 4.5: Extract clinical information using local LLM
        aggregated_stats = None
        if pubmed_results and ENABLE_LOCAL_CLINICAL_EXTRACTION and not args.skip_clinical_extraction:
            logger.info("=" * 70)
            logger.info("Step 4.5: Extracting clinical information with local LLM...")
            logger.info("=" * 70)
            
            try:
                from utils.local_clinical_extractor import LocalClinicalExtractor
                
                # 使用命令行參數或配置檔案的設定
                backend = args.llm_backend or LOCAL_LLM_BACKEND
                model = args.llm_model or LOCAL_LLM_MODEL
                
                logger.info(f"Using LLM backend: {backend}")
                logger.info(f"Using model: {model}")
                logger.info(f"Tensor parallel size: {LOCAL_LLM_TENSOR_PARALLEL}")
                logger.info(f"Max model length: {LOCAL_LLM_MAX_MODEL_LEN} tokens")
                logger.info(f"Max context length: {LOCAL_LLM_MAX_CONTEXT_LENGTH} chars")
                
                extractor = LocalClinicalExtractor(
                    backend=backend,
                    model_name=model,
                    tensor_parallel_size=LOCAL_LLM_TENSOR_PARALLEL,
                    max_model_len=LOCAL_LLM_MAX_MODEL_LEN,
                    max_context_length=LOCAL_LLM_MAX_CONTEXT_LENGTH
                )
                
                # 批次擷取臨床資訊
                # *** 重要：傳入 variant_info 以確保只提取目標變異的信息 ***
                logger.info(f"Processing {len(pubmed_results)} articles...")
                logger.info(f"Target variant: chr{variant_info['chrom']}:{variant_info['pos']}:{variant_info['ref']}>{variant_info['alt']}")
                pubmed_results = extractor.batch_extract(pubmed_results, variant_info, clinvar_info)
                
                # 統計分析
                aggregated_stats = extractor.aggregate_stats(pubmed_results)
                
                logger.info("Clinical Information Extraction Complete!")
                logger.info(f"  Total articles: {aggregated_stats['total_articles']}")
                logger.info(f"  Successfully extracted: {aggregated_stats['articles_with_extraction']}")
                logger.info(f"  Total patients: {aggregated_stats['total_patients']}")
                logger.info(f"  Total families: {aggregated_stats['total_families']}")
                logger.info(f"  Inheritance patterns: {aggregated_stats['inheritance_distribution']}")
                logger.info(f"  Age onset distribution: {aggregated_stats['age_onset_distribution']}")
                logger.info(f"  Affected tissues: {aggregated_stats['tissue_affected']}")
                logger.info(f"  Cardiac phenotypes: {aggregated_stats['cardiac_phenotypes']}")
                logger.info(f"  Skeletal phenotypes: {aggregated_stats['skeletal_phenotypes']}")
                
            except ImportError as e:
                logger.warning(f"Could not import local clinical extractor: {e}")
                logger.warning("Skipping clinical extraction. Install dependencies: pip install vllm transformers")
            except Exception as e:
                logger.error(f"Error during clinical extraction: {e}", exc_info=True)
                logger.warning("Continuing without clinical extraction...")
        elif not pubmed_results:
            logger.info("Step 4.5: Skipped (no PubMed results)")
        elif not ENABLE_LOCAL_CLINICAL_EXTRACTION:
            logger.info("Step 4.5: Skipped (disabled in config)")
        else:
            logger.info("Step 4.5: Skipped (--skip-clinical-extraction flag)")
        
        # Step 5: Generate Protein Schematic (with integrated transcript intervals)
        logger.info("Step 5: Generating protein schematic...")
        image_gen = ImageGenerator()
        
        # Check for transcript intervals file
        xlsx_file = PROJECT_ROOT / "transcript_interval.xlsx"
        xlsx_path = str(xlsx_file) if xlsx_file.exists() else None
        
        if xlsx_path:
            logger.info(f"Using transcript intervals from: {xlsx_file}")
        else:
            logger.info("Transcript intervals file not found, using default visualization")
        
        image_path = image_gen.generate_titin_schematic(variant_info, xlsx_path)
        logger.info(f"Generated image: {image_path}")
        
        # Step 6: Generate HTML Report
        logger.info("Step 6: Generating HTML report...")
        report_gen = HTMLReportGenerator()
        
        # Determine output path
        if args.output:
            output_path = Path(args.output)
        else:
            # 使用 variant_id 作為檔名（格式：chrom-pos-ref-alt）
            variant_id = f"{variant_info['chrom']}-{variant_info['pos']}-{variant_info['ref']}-{variant_info['alt']}"
            
            # 取得 LLM 模型名稱並簡化（去除路徑前綴）
            if ENABLE_LOCAL_CLINICAL_EXTRACTION and not args.skip_clinical_extraction and pubmed_results:
                model_name = args.llm_model or LOCAL_LLM_MODEL
                model_name_short = model_name.split('/')[-1] if '/' in model_name else model_name
                # 清理模型名稱中的特殊字元，確保檔名合法
                model_name_short = model_name_short.replace(':', '-').replace(' ', '_')
            else:
                model_name_short = "no-llm"
            
            output_path = OUTPUT_DIR / f"{variant_id}_{model_name_short}.html"
        
        # Pass aggregated stats to report generator
        import inspect
        sig = inspect.signature(report_gen.generate_report)
        if 'aggregated_stats' in sig.parameters:
            report_path = report_gen.generate_report(
                variant_info=variant_info,
                evo2_result=evo2_result,
                pubmed_results=pubmed_results,
                image_path=image_path,
                output_path=output_path,
                clinvar_info=clinvar_info,
                aggregated_stats=aggregated_stats
            )
        else:
            # Fallback for older html_report.py version
            report_path = report_gen.generate_report(
                variant_info=variant_info,
                evo2_result=evo2_result,
                pubmed_results=pubmed_results,
                image_path=image_path,
                output_path=output_path,
                clinvar_info=clinvar_info
            )
        
        logger.info(f"Report generated successfully: {report_path}")
        print(f"\n{'='*70}")
        print(f"SUCCESS! Report generated at:")
        print(f"   {report_path}")
        if aggregated_stats:
            print(f"\nClinical Extraction Summary:")
            print(f"   Articles analyzed: {aggregated_stats['total_articles']}")
            print(f"   Successfully extracted: {aggregated_stats['articles_with_extraction']}")
            if aggregated_stats.get('disease_distribution'):
                print(f"   Diseases found: {len(aggregated_stats['disease_distribution'])}")
            if aggregated_stats.get('inheritance_distribution'):
                print(f"   Inheritance patterns: {len(aggregated_stats['inheritance_distribution'])}")
        print(f"{'='*70}\n")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error in pipeline: {str(e)}", exc_info=True)
        print(f"\nError: {str(e)}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())