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
from typing import Optional

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
    TRITON_PTXAS_PATH,
)

# LitVar 結果上限（舊 config 沒有此項時 fallback 到 50）
try:
    from config import LITVAR_MAX_RESULTS
except ImportError:
    LITVAR_MAX_RESULTS = 50

# Try to get ENABLE_EASY_PROMPT setting
try:
    from config import ENABLE_EASY_PROMPT
except ImportError:
    ENABLE_EASY_PROMPT = False

# Track B: mechanistic / analogy-based inference（無文獻時的類比推理）
try:
    from config import (
        ENABLE_MECHANISTIC_INFERENCE,
        INFERENCE_TOP_K_ANALOGS,
        HGMD_KB_PATH,
        TTN_BAND_EDGES,
        TTN_BAND_LABELS,
        TRANSCRIPT_INTERVAL_PATH,
        TTN_TRANSCRIPT_MODEL_PATH,
        HGVS_ALLOW_NETWORK,
        REFERENCE_GENOME_PATH,
        LLM_GPU_MEMORY_UTILIZATION,
        INFERENCE_HGMD_CSV_PATH,
        INFERENCE_CARDIODB_PATH,
        INFERENCE_PREDICT_MODE,
        INFERENCE_M_EXON,
        INFERENCE_M_BUCKET,
        INFERENCE_M_BAND,
    )
except ImportError:
    ENABLE_MECHANISTIC_INFERENCE = False
    INFERENCE_TOP_K_ANALOGS = 8
    HGMD_KB_PATH = None
    INFERENCE_HGMD_CSV_PATH = None
    INFERENCE_CARDIODB_PATH = None
    INFERENCE_PREDICT_MODE = "analogy"
    INFERENCE_M_EXON = 4.0
    INFERENCE_M_BUCKET = 4.0
    INFERENCE_M_BAND = 8.0
    TTN_BAND_EDGES = [178525989, 178535982, 178618584, 178775356, 178807423]
    TTN_BAND_LABELS = ["M-band", "A-band", "I-band", "Z-disk"]
    TRANSCRIPT_INTERVAL_PATH = None
    TTN_TRANSCRIPT_MODEL_PATH = None
    HGVS_ALLOW_NETWORK = False
    REFERENCE_GENOME_PATH = None
    LLM_GPU_MEMORY_UTILIZATION = 0.85

try:
    from config import ENABLE_INHERITANCE_INFERENCE
except ImportError:
    ENABLE_INHERITANCE_INFERENCE = False

# 資料庫查表：HGMD / PubTator / PubMind 對這顆變異登記的疾病（純查表，非推論）
try:
    from config import (
        ENABLE_DATABASE_ANNOTATIONS,
        HGMD_FULL_CSV_PATH,
        PUBTATOR_VARIANTS_CSV_PATH,
        PUBTATOR_COORD_CACHE_PATH,
        PUBMIND_JSON_PATH,
    )
except ImportError:
    ENABLE_DATABASE_ANNOTATIONS = False
    HGMD_FULL_CSV_PATH = None
    PUBTATOR_VARIANTS_CSV_PATH = None
    PUBTATOR_COORD_CACHE_PATH = None
    PUBMIND_JSON_PATH = None

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

# Pre-import vLLM before Evo2/PyTorch to avoid dynamic linker symbol conflicts.
# vLLM's _C.abi3.so must be loaded into the process before torch initialises
# its CUDA runtime; once torch is live the required symbols are shadowed and
# the import fails with "undefined symbol".
_VLLM_PRELOAD_OK = False
if ENABLE_LOCAL_CLINICAL_EXTRACTION and LOCAL_LLM_BACKEND == "vllm":
    try:
        import vllm as _vllm_preload  # noqa: F401
        _VLLM_PRELOAD_OK = True
    except Exception as _e:
        import warnings
        warnings.warn(
            f"vLLM pre-import failed (clinical extraction will be skipped): {_e}"
        )

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


def gpu_headroom() -> tuple:
    """回傳 (可用 GiB, 總量 GiB, 佔用者字串)；查不到時為 (None, None, "")。

    這台機器是 GB10 統一記憶體，nvidia-smi 的 memory.free 會回 [N/A]，
    只有 torch.cuda.mem_get_info 給得出正確數字，所以以它為準；
    佔用者名單再從 nvidia-smi 的 compute-apps 取得，純粹是為了讓錯誤訊息有用。
    """
    try:
        import torch

        if not torch.cuda.is_available():
            return (None, None, "")
        free, total = torch.cuda.mem_get_info(0)
    except Exception:  # noqa: BLE001
        return (None, None, "")

    occupants = ""
    try:
        import subprocess

        out = subprocess.run(
            ["nvidia-smi", "--query-compute-apps=pid,used_memory", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip()
        rows = []
        for line in out.splitlines():
            pid = line.split(",")[0].strip()
            who = subprocess.run(["ps", "-o", "user=,cmd=", "-p", pid],
                                 capture_output=True, text=True, timeout=5).stdout.strip()
            rows.append(f"{line.strip()} ({who[:80] or 'unknown'})")
        occupants = "; ".join(rows)
    except Exception:  # noqa: BLE001
        pass
    return (free / 2**30, total / 2**30, occupants)


# ── 共用 LLM（batch_run_variants.py --in-process）────────────────────────────
# 模型載入一次要 4–5 分鐘；同一個程序內跑多顆變異時，Step 4.5 與 Track B
# 都從這裡拿同一個 extractor，不再逐筆重載。
_SHARED_EXTRACTOR = None
_SHARED_EXTRACTOR_KEY = None

# run_pipeline 的回傳碼：vLLM 引擎已壞，此程序內無法再產出有 LLM 結果的報告
EXIT_LLM_ENGINE_DEAD = 3

_FATAL_LLM_ERROR_NAMES = {"EngineDeadError", "OutOfMemoryError"}


def _is_fatal_llm_error(exc) -> bool:
    """引擎掛掉或 CUDA 出錯時，同一個程序內救不回來（沿 __cause__/__context__ 往上找）。"""
    seen = set()
    while exc is not None and id(exc) not in seen:
        seen.add(id(exc))
        if type(exc).__name__ in _FATAL_LLM_ERROR_NAMES or "CUDA error" in str(exc):
            return True
        exc = exc.__cause__ or exc.__context__
    return False


def _watch_llm_engine(extractor) -> None:
    """包住 llm.generate，記下致命錯誤。

    batch_extract 與 Track B 都會吞掉 generate 的例外並退回規則式結果；
    單筆 subprocess 時只影響一筆，in-process 時卻會讓之後每一筆都默默少了
    LLM 結果，所以要在這裡留下記號，由 run_pipeline 決定中止。
    """
    extractor.engine_dead = False
    llm = getattr(extractor, "llm", None)
    if llm is None:
        return
    original_generate = llm.generate

    def generate(*args, **kwargs):
        try:
            return original_generate(*args, **kwargs)
        except Exception as e:
            if _is_fatal_llm_error(e):
                extractor.engine_dead = True
            raise

    llm.generate = generate


def llm_engine_dead() -> bool:
    return _SHARED_EXTRACTOR is not None and getattr(_SHARED_EXTRACTOR, "engine_dead", False)


def llm_loaded(backend: str, model: str) -> bool:
    return _SHARED_EXTRACTOR is not None and _SHARED_EXTRACTOR_KEY == (backend, model)


def release_gpu_cache(need_gib: float) -> None:
    """釋放檔案快取（scripts/release_page_cache.py），讓 vLLM 的顯存預檢看得到足夠的可用記憶體。"""
    import logging
    import subprocess

    logger = logging.getLogger(__name__)
    env = dict(os.environ, NEED_GIB=f"{need_gib:.1f}")
    try:
        out = subprocess.run([sys.executable, str(Path(__file__).parent / "scripts" / "release_page_cache.py")],
                             env=env, capture_output=True, text=True, timeout=600)
        logger.info((out.stdout or out.stderr).strip() or "page cache released")
    except Exception as e:  # noqa: BLE001
        logger.warning(f"Could not release the page cache: {e}")


def _main_disease(mechanistic_result, pubmed_results) -> Optional[str]:
    """遺傳模式推論用的疾病：Track B 的推論疾病，否則取直接提到變異的文獻中最常擷取到的疾病。"""
    inf = (mechanistic_result or {}).get("inference") or {}
    disease = str(inf.get("inferred_disease") or "").strip()
    if disease and disease.lower() != "uncertain":
        return disease
    counts = {}
    for a in pubmed_results or []:
        d = str(((a.get("clinical_info") or {}).get("disease")) or "").strip()
        if d and d.lower() not in ("not specified", "n/a", "none", "json parsing failed"):
            counts[d] = counts.get(d, 0) + 1
    return max(counts, key=counts.get) if counts else None


_MECHANISTIC_ENGINE = None


def build_mechanistic_engine():
    """Track B 與遺傳模式推論共用的查表引擎（外顯子／band／PSI／HGMD 類比）；同一程序只建一次。"""
    global _MECHANISTIC_ENGINE
    if _MECHANISTIC_ENGINE is None:
        from utils.mechanistic_inference import MechanisticInference

        _MECHANISTIC_ENGINE = MechanisticInference(
            kb_path=HGMD_KB_PATH,
            band_edges=TTN_BAND_EDGES,
            band_labels=TTN_BAND_LABELS,
            top_k=INFERENCE_TOP_K_ANALOGS,
            transcript_interval_path=TRANSCRIPT_INTERVAL_PATH,
            transcript_model_path=TTN_TRANSCRIPT_MODEL_PATH,
            reference_fasta_path=REFERENCE_GENOME_PATH,
            allow_network=HGVS_ALLOW_NETWORK,
            # 新資料來源：cardiodb 的 364 個外顯子（座標／band／實測 PSI）
            # 與 HGMD 完整匯出（2,193 筆）。兩者皆為 None 時退回舊行為。
            cardiodb_path=INFERENCE_CARDIODB_PATH,
            hgmd_csv_path=INFERENCE_HGMD_CSV_PATH,
            predict_mode=INFERENCE_PREDICT_MODE,
            m_exon=INFERENCE_M_EXON,
            m_bucket=INFERENCE_M_BUCKET,
            m_band=INFERENCE_M_BAND,
        )
    return _MECHANISTIC_ENGINE


def get_shared_extractor(backend: str, model: str):
    """同一程序內只載入一次 LLM；之後的呼叫直接回傳同一個 extractor。"""
    global _SHARED_EXTRACTOR, _SHARED_EXTRACTOR_KEY
    if llm_loaded(backend, model):
        if llm_engine_dead():
            raise RuntimeError("shared vLLM engine is dead; restart the process")
        return _SHARED_EXTRACTOR

    from utils.local_clinical_extractor import LocalClinicalExtractor

    extractor = LocalClinicalExtractor(
        backend=backend,
        model_name=model,
        tensor_parallel_size=LOCAL_LLM_TENSOR_PARALLEL,
        max_model_len=LOCAL_LLM_MAX_MODEL_LEN,
        max_context_length=LOCAL_LLM_MAX_CONTEXT_LENGTH,
    )
    _watch_llm_engine(extractor)
    _SHARED_EXTRACTOR, _SHARED_EXTRACTOR_KEY = extractor, (backend, model)
    return extractor


def _any_article_mentions_variant(pubmed_results) -> bool:
    """判斷是否有任何文獻「直接提到」輸入變異（決定是否需要走 Track B 推理）。

    與 html_report._article_mentions_variant 邏輯一致：只要有一篇文章擷取到
    非空的 evidence_sentence、具體 disease、或 patient_count > 0，就視為有文獻證據。
    """
    if not pubmed_results:
        return False
    skip = ("not specified", "n/a", "none", "", "json parsing failed")
    for article in pubmed_results:
        info = (article.get("clinical_info") or {}) if isinstance(article, dict) else {}
        evidence = str(info.get("evidence_sentence") or "").strip().lower()
        if evidence and evidence not in skip:
            return True
        disease = str(info.get("disease") or "").strip().lower()
        if disease and disease not in skip:
            return True
        try:
            if int(info.get("patient_count") or 0) > 0:
                return True
        except (ValueError, TypeError):
            pass
    return False


def _parse_variant_input(raw: str, logger) -> dict:
    """接受 `2-178612477-T-A`，也接受文獻慣用的 HGVS（c./p./g.）。

    chr-pos-ref-alt 用轉錄方向（與 c. 寫法同一股）；基因組正股（VCF）寫法也能解析，會自動轉換。
    由 utils/variant_ids 用參考序列判斷。回傳的 variant_info 一律是轉錄方向
    （variant_id、ref、alt），另附 display_id / vcf_id / hgvs_c_display 給報告顯示。
    HGVS 由 utils/hgvs_resolver 用 TTN 轉錄本外顯子模型離線換算成 GRCh38 座標。
    """
    from utils import variant_ids
    from utils.variant_parser import VariantParseError

    if variant_ids.ID_RE.match(raw or ""):
        try:
            v = variant_ids.normalize(raw)
        except variant_ids.VariantIdError as e:
            raise VariantParseError(str(e)) from e
        if v["orientation"] == "plus":
            logger.warning(f"輸入 {raw} 是基因組正股（VCF）寫法；已轉成 TTN 轉錄方向 {v['internal_id']}")
        info = {
            "chrom": v["chrom"], "pos": v["pos"], "ref": v["tx_ref"], "alt": v["tx_alt"],
            "variant_id": v["internal_id"], "input_orientation": v["orientation"],
        }
        return variant_ids.add_display_fields(info)

    try:
        return variant_ids.add_display_fields(parse_variant(raw))
    except VariantParseError as e:
        from utils.hgvs_resolver import HgvsResolveError, HgvsResolver

        try:
            resolved = HgvsResolver(
                model_path=TTN_TRANSCRIPT_MODEL_PATH,
                fasta_path=REFERENCE_GENOME_PATH,
                allow_network=HGVS_ALLOW_NETWORK,
            ).resolve(raw)
        except HgvsResolveError as he:
            raise VariantParseError(f"{e}；也無法當成 HGVS 解析（{he}）") from he
        info = resolved.to_variant_info()
        if not info.get("ref") or not info.get("alt"):
            raise VariantParseError(
                f"{raw} 只能定位到 codon 起點（{info['variant_id']}），"
                "缺少 ref/alt 無法跑 Evo2；請改用 c.HGVS 或 chr-pos-ref-alt"
            )
        logger.info(
            f"HGVS 解析：{raw} -> {info['variant_id']} "
            f"（{info.get('transcript')} {info.get('hgvs_c') or info.get('hgvs_p')}, "
            f"{info.get('hgvs_vclass')}, confidence={info.get('hgvs_confidence')}）"
        )
        for w in info.get("hgvs_warnings") or []:
            logger.warning(f"HGVS 解析提醒：{w}")
        return variant_ids.add_display_fields(info)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='TTN Variant AI Agent - Comprehensive variant analysis'
    )
    parser.add_argument(
        'variant',
        type=str,
        help=(
            'Variant as chromosome-position-ref-alt (e.g., 2-178612477-T-A) '
            'or HGVS (e.g., c.98671-2A>G, NM_001267550.2:c.29800T>C, p.R31056*)'
        )
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
    parser.add_argument(
        '--require-llm',
        action='store_true',
        help=('Abort with a non-zero exit code if clinical extraction was expected '
              'but the LLM could not run (e.g. the GPU is busy). Without this the run '
              'continues and the report is labelled no-llm.')
    )
    return parser


def main():
    logger = setup_logging()
    args = build_arg_parser().parse_args()
    return run_pipeline(args, logger)


def run_pipeline(args: argparse.Namespace, logger) -> int:
    """Main pipeline orchestrator（batch_run_variants.py --in-process 會直接呼叫）"""
    logger.info(f"Starting TTN Variant AI Agent for variant: {args.variant}")
    
    try:
        # Step 1: Parse variant（chr-pos-ref-alt，或 HGVS c./p./g.）
        logger.info("Step 1: Parsing variant...")
        variant_info = _parse_variant_input(args.variant, logger)
        logger.info(f"Parsed variant: {variant_info}")
        
        # Step 2: Parse ClinVar Information
        logger.info("Step 2: Parsing ClinVar information...")
        clinvar_parser = ClinVarParser()
        clinvar_info = clinvar_parser.parse_variant(variant_info)
        if clinvar_info and clinvar_info.get('pmid_list'):
            logger.info(f"Found {len(clinvar_info['pmid_list'])} PubMed IDs from ClinVar")
            logger.info(f"Variant impact: {clinvar_info.get('variant_impact')}")
            logger.info(f"Clinical significance: {clinvar_info.get('clinical_significance')}")
        elif clinvar_info:
            # ClinVar 有此變異記錄，但沒有引用文章。
            # 保留 clinvar_info（conditions / HGVS / rsid），供 Track B 類比推理與
            # LitVar(rsID) 使用；不再整個丟棄。文獻檢索仍會因無 pmid_list 而略過。
            logger.warning(
                "ClinVar 有記錄但無引用文章 —— 保留 conditions/HGVS/rsID 供類比推理使用"
            )
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
                max_results=LITVAR_MAX_RESULTS
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
        clinical_extractor = None  # 供 Track B（Step 4.6）重用已載入的 LLM
        backend = args.llm_backend or LOCAL_LLM_BACKEND
        model = args.llm_model or LOCAL_LLM_MODEL
        vllm_ok = _VLLM_PRELOAD_OK if backend == "vllm" else True
        # 只有真的跑完 batch_extract + aggregate_stats 才算成功；報告檔名依此命名。
        clinical_extraction_ok = False
        gpu_blocked = False
        if pubmed_results and ENABLE_LOCAL_CLINICAL_EXTRACTION and not args.skip_clinical_extraction:
            logger.info("=" * 70)
            logger.info("Step 4.5: Extracting clinical information with local LLM...")
            logger.info("=" * 70)

            if backend == "vllm" and not vllm_ok:
                logger.error(
                    "vLLM pre-import failed at startup — clinical extraction skipped. "
                    "Check that vLLM is installed and compatible with the current PyTorch build."
                )

            # 顯存不夠時 vLLM 會丟出一長串 CUDA out of memory 堆疊，看不出是誰佔走的。
            # 這裡先擋下來，直接說明還差多少、現在是誰在用。
            # 模型已在本程序載入時，顯存本來就被自己佔著，不能再檢查。
            free_gib, total_gib, occupants = (
                (None, None, "") if llm_loaded(backend, model) else gpu_headroom()
            )
            if free_gib is not None and free_gib < total_gib * LLM_GPU_MEMORY_UTILIZATION:
                # GB10 共用記憶體：Evo2 權重、PDF 等的檔案快取會被算成 GPU 已使用，先清再判斷
                release_gpu_cache(total_gib * LLM_GPU_MEMORY_UTILIZATION)
                free_gib, total_gib, occupants = gpu_headroom()
            if free_gib is not None:
                needed = total_gib * LLM_GPU_MEMORY_UTILIZATION
                logger.info(f"GPU headroom: {free_gib:.1f} / {total_gib:.1f} GiB free "
                            f"(vLLM will request {needed:.1f} GiB)")
                if free_gib < needed:
                    logger.error(
                        f"Not enough free GPU memory to load {LOCAL_LLM_MODEL}: "
                        f"{free_gib:.1f} GiB free but {needed:.1f} GiB required."
                    )
                    if occupants:
                        logger.error(f"GPU currently held by: {occupants}")
                    gpu_blocked = True

            try:
                if gpu_blocked:
                    raise RuntimeError("insufficient free GPU memory (pre-flight check)")

                # 使用命令行參數或配置檔案的設定（backend / model 已在外層定義）
                logger.info(f"Using LLM backend: {backend}")
                logger.info(f"Using model: {model}")
                logger.info(f"Tensor parallel size: {LOCAL_LLM_TENSOR_PARALLEL}")
                logger.info(f"Max model length: {LOCAL_LLM_MAX_MODEL_LEN} tokens")
                logger.info(f"Max context length: {LOCAL_LLM_MAX_CONTEXT_LENGTH} chars")
                
                extractor = get_shared_extractor(backend, model)
                clinical_extractor = extractor  # 保留給 Track B 重用
                
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
                clinical_extraction_ok = True
                
            except ImportError as e:
                logger.warning(f"Could not import local clinical extractor: {e}")
                logger.warning("Skipping clinical extraction. Install dependencies: pip install vllm transformers")
            except Exception as e:
                logger.error(f"Error during clinical extraction: {e}", exc_info=True)
                logger.warning("Continuing without clinical extraction...")

            if not clinical_extraction_ok:
                # 這裡曾經是個沉默的陷阱：擷取失敗後流程照跑，報告卻仍以模型名
                # 命名，看起來像是 LLM 分析過。批次跑時會產出一整批標錯的報告。
                logger.error(
                    "Clinical extraction did not complete — the report will be "
                    "labelled 'no-llm' and will contain no AI Analysis section."
                )
                if args.require_llm:
                    logger.error("--require-llm was set, aborting instead of writing a report.")
                    return 1
        elif not pubmed_results:
            logger.info("Step 4.5: Skipped (no PubMed results)")
        elif not ENABLE_LOCAL_CLINICAL_EXTRACTION:
            logger.info("Step 4.5: Skipped (disabled in config)")
        else:
            logger.info("Step 4.5: Skipped (--skip-clinical-extraction flag)")

        # Step 4.6: Mechanistic / analogy-based inference (Track B)
        # 當沒有任何文獻直接提到輸入變異時，切換到類比推理路徑：
        # Evo2 + 結構 band + 變異類型 + HGMD 同區帶同類型的已知致病變異 → LLM 推論。
        # 輸出明確標示為 inference（非文獻證據），與 evidence 軌在報告上分軌呈現。
        mechanistic_result = None
        literature_mentions_variant = _any_article_mentions_variant(pubmed_results)
        if ENABLE_MECHANISTIC_INFERENCE and not literature_mentions_variant:
            logger.info("=" * 70)
            logger.info("Step 4.6: No literature mentions the variant -> mechanistic inference (Track B)")
            logger.info("=" * 70)
            try:
                mi = build_mechanistic_engine()

                # 若 Step 4.5 沒建過 LLM（例如完全無文獻），這裡按需載入一次
                inference_extractor = clinical_extractor
                if (
                    inference_extractor is None
                    and ENABLE_LOCAL_CLINICAL_EXTRACTION
                    and not args.skip_clinical_extraction
                    and vllm_ok
                ):
                    try:
                        free_gib, total_gib, _ = gpu_headroom()
                        if free_gib is not None and free_gib < total_gib * LLM_GPU_MEMORY_UTILIZATION:
                            release_gpu_cache(total_gib * LLM_GPU_MEMORY_UTILIZATION)
                        logger.info(f"Loading LLM for inference (backend={backend}, model={model})...")
                        inference_extractor = get_shared_extractor(backend, model)
                        clinical_extractor = inference_extractor
                    except Exception as e:
                        logger.warning(f"Could not load LLM for inference, using rule-based fallback: {e}")
                        inference_extractor = None

                mechanistic_result = mi.infer(
                    inference_extractor, variant_info, clinvar_info, evo2_result
                )

                feats = mechanistic_result.get("features", {})
                inf = mechanistic_result.get("inference", {})
                logger.info(
                    f"  Band: {feats.get('band')} | Type: {feats.get('vclass')} "
                    f"({feats.get('vgroup')}) | Analogs: {len(feats.get('analogs', []))}"
                )
                logger.info(
                    f"  Inferred disease: {inf.get('inferred_disease')} "
                    f"| tissue: {inf.get('inferred_tissue')} "
                    f"| confidence: {inf.get('confidence')} "
                    f"| source: {inf.get('source')}"
                )
            except Exception as e:
                logger.error(f"Error during mechanistic inference: {e}", exc_info=True)
                logger.warning("Continuing without mechanistic inference...")
        elif not ENABLE_MECHANISTIC_INFERENCE:
            logger.info("Step 4.6: Skipped (mechanistic inference disabled in config)")
        else:
            logger.info("Step 4.6: Skipped (literature evidence for the variant is available)")

        if llm_engine_dead():
            # 擷取與推論都已吞掉錯誤、退回規則式結果；照常出報告會混進一份
            # 看似有跑 LLM 的結果，所以不寫報告，交給呼叫端重啟程序。
            logger.error("vLLM engine died while processing this variant — not writing a report.")
            return EXIT_LLM_ENGINE_DEAD

        # Step 5: Generate Protein Schematic (with integrated transcript intervals)
        logger.info("Step 5: Generating protein schematic...")
        image_gen = ImageGenerator()
        
        # Check for transcript intervals file
        xlsx_file = TRANSCRIPT_INTERVAL_PATH
        xlsx_path = str(xlsx_file) if xlsx_file.exists() else None
        
        if xlsx_path:
            logger.info(f"Using transcript intervals from: {xlsx_file}")
        else:
            logger.info("Transcript intervals file not found, using default visualization")
        
        image_path = image_gen.generate_titin_schematic(variant_info, xlsx_path)
        logger.info(f"Generated image: {image_path}")
        
        # Step 5.5: Curated database lookup（HGMD / PubTator 有沒有這顆變異、登記什麼疾病）
        db_annotations = None
        if ENABLE_DATABASE_ANNOTATIONS:
            logger.info("Step 5.5: Looking up HGMD / PubTator / PubMind for this variant...")
            try:
                from utils.variant_db_lookup import get_shared_lookup

                db_lookup = get_shared_lookup(
                    hgmd_csv=HGMD_FULL_CSV_PATH,
                    hgmd_kb=HGMD_KB_PATH,
                    pubtator_csv=PUBTATOR_VARIANTS_CSV_PATH,
                    cache_path=PUBTATOR_COORD_CACHE_PATH,
                    pubmind_json=PUBMIND_JSON_PATH,
                )
                db_annotations = db_lookup.lookup(
                    chrom=str(variant_info['chrom']),
                    pos=int(variant_info['pos']),
                    ref=variant_info['ref'],
                    alt=variant_info['alt'],
                )
                logger.info(
                    f"  HGMD: {db_annotations.hgmd_diseases() or 'no match'} "
                    f"| PubTator: {db_annotations.pubtator_diseases() or 'no match'} "
                    f"| PubMind: {db_annotations.pubmind_diseases() or 'no match'}"
                )
            except Exception as e:
                logger.error(f"Error during database lookup: {e}", exc_info=True)
                logger.warning("Continuing without database annotations...")

        # Step 5.6: 遺傳模式推論（AD / AR）。用到 Step 5.5 的資料庫紀錄，所以放在它之後
        inheritance_result = None
        if ENABLE_INHERITANCE_INFERENCE:
            logger.info("Step 5.6: Inferring mode of inheritance (AD / AR)...")
            try:
                from utils.inheritance_inference import InheritanceInference

                mi = build_mechanistic_engine()
                features = (mechanistic_result or {}).get("features") or mi.build_features(
                    variant_info, clinvar_info, evo2_result
                )
                inheritance_result = InheritanceInference(mi, INFERENCE_HGMD_CSV_PATH).infer(
                    clinical_extractor if clinical_extractor is not None and not llm_engine_dead() else None,
                    features,
                    db_annotations=db_annotations,
                    pubmed_results=pubmed_results,
                    predicted_disease=_main_disease(mechanistic_result, pubmed_results),
                )
                logger.info(
                    f"  Inheritance: {inheritance_result['inheritance']} "
                    f"({inheritance_result['confidence']}, {inheritance_result['source']})"
                )
            except Exception as e:
                logger.error(f"Error during inheritance inference: {e}", exc_info=True)
            if llm_engine_dead():
                logger.error("vLLM engine died during inheritance inference — not writing a report.")
                return EXIT_LLM_ENGINE_DEAD

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
            if clinical_extraction_ok:
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
            report_kwargs = dict(
                variant_info=variant_info,
                evo2_result=evo2_result,
                pubmed_results=pubmed_results,
                image_path=image_path,
                output_path=output_path,
                clinvar_info=clinvar_info,
                aggregated_stats=aggregated_stats,
            )
            # 只有在報告產生器支援時才傳入 Track B 結果（向後相容）
            if 'mechanistic_result' in sig.parameters:
                report_kwargs['mechanistic_result'] = mechanistic_result
            if 'db_annotations' in sig.parameters:
                report_kwargs['db_annotations'] = db_annotations
            if 'inheritance_result' in sig.parameters:
                report_kwargs['inheritance_result'] = inheritance_result
            report_path = report_gen.generate_report(**report_kwargs)
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