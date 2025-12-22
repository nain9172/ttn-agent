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

# Evo2
EVO2_MODEL = "evo2_1b_base"
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
MAX_TEXT_LENGTH = 15000  # Increased for Llama 3 analysis
ENABLE_LITVAR_SEARCH = True

# Docling PDF Processing
ENABLE_DOCLING_PDF = True
DOCLING_MAX_PRIORITY_LENGTH = 8000  # 優先內容最大長度（節省 GPU 記憶體）
DOCLING_ENABLE_OCR = True  # 啟用 OCR 提取圖片中的表格
DOCLING_ENABLE_TABLE_EXTRACTION = True  # 啟用表格提取
DOCLING_OUTPUT_DIR = OUTPUT_DIR / "markdown"  # Markdown 輸出目錄
PDF_OUTPUT_DIR = OUTPUT_DIR / "pdfs"  # PDF 下載目錄

# LLM
LOCAL_LLM_BACKEND = "vllm" # or ollama
LOCAL_LLM_MODEL = "meta-llama/Llama-3.2-3B-Instruct" # Supports 3B or 8B
LOCAL_LLM_TENSOR_PARALLEL = 2
ENABLE_LOCAL_CLINICAL_EXTRACTION = True

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

PHENOTYPE_CATEGORIES = {
    "heart": ["cardiomyopathy", "DCM", "HCM", "ARVC", "cardiac"],
    "skeletal_muscle": ["myopathy", "muscular dystrophy", "LGMD"],
    "both": ["heart", "skeletal"]
}