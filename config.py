"""
Configuration file
"""
import os
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
DATA_DIR = PROJECT_ROOT / "data"

OUTPUT_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)
LOG_FILE = OUTPUT_DIR / "ttn_agent.log"
TRITON_PTXAS_PATH = "/usr/local/cuda/bin/ptxas"
# Evo2
# EVO2_MODEL = "evo2_1b_base"
EVO2_MODEL = "evo2_7b_base"
EVO2_WINDOW_SIZE = 8192
PATHOGENIC_THRESHOLD = 0.0
REFERENCE_GENOME_PATH = DATA_DIR / "sequence.fasta"
TTN_SEQUENCE_START = 178807423
TTN_SEQUENCE_END = 178525989

# PubMed
PUBMED_EMAIL = os.getenv("PUBMED_EMAIL", "your.email@example.com")
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY", "")
PUBMED_MAX_RESULTS = 50
ENABLE_FULL_TEXT_FETCH = True
MAX_TEXT_LENGTH = 80000  # 增加到 80K 字元，充分利用 MedGemma 的長 context 能力（最高 128K tokens）
ENABLE_LITVAR_SEARCH = True

# Docling PDF Processing
ENABLE_DOCLING_PDF = True
DOCLING_MAX_PRIORITY_LENGTH = 120000  # 優先內容最大長度（增加到 120K 以充分利用 MedGemma 的長 context）
DOCLING_ENABLE_OCR = True  # 啟用 OCR 提取圖片中的表格
DOCLING_ENABLE_TABLE_EXTRACTION = True  # 啟用表格提取
DOCLING_OUTPUT_DIR = OUTPUT_DIR / "markdown"  # Markdown 輸出目錄（快取：已轉換的 .md 檔案會保存在此）
PDF_OUTPUT_DIR = OUTPUT_DIR / "pdfs"  # PDF 下載目錄（快取：已下載的 PDF 會保存在此）

# Supplementary Files Processing
DOWNLOAD_SUPPLEMENTARY_FILES = True  # 是否下載並解析 supplementary files（Excel 等）
MAX_SUPPLEMENTARY_FILES = 3  # 最多下載的 supplementary files 數量（避免太慢）
# 允許的 supplementary files 格式（只下載這些格式的檔案）
ALLOWED_SUPPLEMENTARY_FORMATS = ['.doc', '.docx', '.xlsx', '.xls', '.csv', '.pdf']
FORCE_REDOWNLOAD_SUPPLEMENTARY = False  # 是否強制重新下載已存在的 supplementary files

# Cache Settings - 快取設定（節省執行時間）
# 注意：所有下載和轉換功能都已實現快取機制，會自動跳過已處理的文件：
# - PDF 下載：如果 outputs/pdfs/PMCXXXX.pdf 已存在，直接使用
# - Markdown 轉換：如果 outputs/markdown/PMCXXXX.md 已存在，直接讀取
# - Supplementary 下載：如果 outputs/supplementary/pmid_XXXX/ 目錄有文件，直接使用
# 如需強制重新處理，可刪除對應的快取檔案或目錄

# LLM
LOCAL_LLM_BACKEND = "vllm" # or ollama
# LOCAL_LLM_MODEL = "google/medgemma-27b-text-it"
LOCAL_LLM_MODEL = "google/medgemma-27b-it"  # MedGemma-27B instruction-tuned (醫療專用，27B 參數，約需 54GB)
# LOCAL_LLM_MODEL = "google/medgemma-1.5-4b-it"  # 更輕量的選擇（約需 8GB）
# LOCAL_LLM_MODEL = "meta-llama/Llama-3.2-8B-Instruct"  # Llama 3.2 8B（約需 16GB）
# LOCAL_LLM_MODEL = "google/gemma-3-27b-it"
LOCAL_LLM_TENSOR_PARALLEL = 1  # 單 GPU 使用 1，多 GPU 可增加（例如 2 或 4）
LOCAL_LLM_MAX_MODEL_LEN = 131072  # vLLM 最大模型長度（65K tokens，MedGemma 支持最高 128K）
LOCAL_LLM_MAX_CONTEXT_LENGTH = 120000  # 準備給 LLM 的文本最大字元數（增加以利用長 context）
ENABLE_LOCAL_CLINICAL_EXTRACTION = True
ENABLE_EASY_PROMPT = False  # 如果為 True，簡化 LLM prompt，不包含文獻內容，並跳過全文下載、Docling 處理、LitVar 搜尋等步驟（快速測試模式）

# TTN Info
TTN_GENE_INFO = {
    "chromosome": "2", "gene_name": "TTN",
    "start": 178807423, "end": 178525989, "strand": "-"
}

TTN_DOMAINS = {
    "Z-disk": {"start": 1, "end": 28, "color":"#FF6B6B"},
    "I-band": {"start": 29, "end": 252, "color": "#4ECDC4"},
    "A-band": {"start": 253, "end": 358, "color": "#45B7D1"},
    "M-band": {"start": 359, "end": 364, "color": "#96CEB4"}
}

TTN_TRANSCRIPTS = {
    "N2BA": {"id": "ENST00000589042", "length": 34350, "description": "Full-length cardiac"},
    "N2B": {"id": "ENST00000460472", "length": 26926, "description": "Shorter cardiac"}
}

REPORT_TITLE = "TTN Variant AI Report"
REPORT_HEADER_COLOR = "#2C3E50"
REPORT_ACCENT_COLOR = "#3498DB"
ENABLE_TISSUE_STATISTICS = True  # 是否在報告中顯示組織影響統計（cardiac, skeletal, both, not specified）===================

PHENOTYPE_CATEGORIES = {
    "heart": ["cardiomyopathy", "DCM", "HCM", "ARVC", "cardiac"],
    "skeletal_muscle": ["myopathy", "muscular dystrophy", "LGMD"],
    "both": ["heart", "skeletal"]
}