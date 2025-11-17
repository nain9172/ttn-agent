"""
Configuration file for TTN Variant AI Agent
Contains constants, API keys, and file paths
"""

from pathlib import Path
import os
CUDA_VISIBLE_DEVICES=1
# Project paths
PROJECT_ROOT = Path(__file__).parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
TEMPLATES_DIR = PROJECT_ROOT / "templates"
DATA_DIR = PROJECT_ROOT / "data"

# Create directories if they don't exist
OUTPUT_DIR.mkdir(exist_ok=True)
TEMPLATES_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)

# Logging
LOG_FILE = OUTPUT_DIR / "ttn_agent.log"

# Evo2 Configuration
EVO2_MODEL = "evo2_1b_base"  # or "evo2_7b_base" for better accuracy
EVO2_WINDOW_SIZE = 8192
PATHOGENIC_THRESHOLD = 0.0  # delta_score < 0 = pathogenic

# Reference genome
# Full genome: DATA_DIR / "GCF_000001405.40_GRCh38.p14_genomic.fna"
# TTN-specific sequence (chr2:178807423-178525989, negative strand)
REFERENCE_GENOME_PATH = DATA_DIR / "sequence.fasta"
TTN_SEQUENCE_START = 178807423  # Genomic coordinate of first base in sequence.fasta
TTN_SEQUENCE_END = 178525989    # Genomic coordinate of last base in sequence.fasta

# PubMed API Configuration
PUBMED_EMAIL = os.getenv("PUBMED_EMAIL", "ryan910702@gmail.com")
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY", '987512860b8bf8e96a09d672b03ff32e2c08')
PUBMED_MAX_RESULTS = 50  # 增加到 50 以涵蓋更多相關文獻
ENABLE_FULL_TEXT_FETCH = True  # 是否嘗試獲取全文
MAX_TEXT_LENGTH = 8000  # LLM 分析的最大文本長度
# Local LLM Configuration
# Backend options: "ollama", "transformers", "vllm"
# 推薦使用 vLLM 以充分利用雙 4090 GPU
LOCAL_LLM_BACKEND = os.getenv("LOCAL_LLM_BACKEND", "vllm")

# Model selection:
# - Llama-3.2-3B: 快速但提取準確度較低（不推薦）
# - Llama-3.1-8B-Instruct: 推薦！較好的指令理解和資訊提取能力
# - Llama-3.1-70B-Instruct: 最佳效果但需要更多記憶體
LOCAL_LLM_MODEL = os.getenv("LOCAL_LLM_MODEL", "meta-llama/Llama-3.2-3B")

# Number of GPUs to use (2 for dual 4090)
LOCAL_LLM_TENSOR_PARALLEL = int(os.getenv("LOCAL_LLM_TENSOR_PARALLEL", "2"))

# Enable/disable local clinical extraction
ENABLE_LOCAL_CLINICAL_EXTRACTION = os.getenv("ENABLE_LOCAL_CLINICAL_EXTRACTION", "true").lower() == "true"

# TTN Gene Information
TTN_GENE_INFO = {
    "chromosome": "2",
    "gene_name": "TTN",
    "gene_id": "7273",
    "start": 178807423,
    "end": 178525989,
    "strand": "-"
}

# TTN Protein Domains (based on exon ranges)
# Position ranges: Z-disk(178807423-178775356), I-band(178775202-178618584), 
#                  A-band(178618491-178535982), M-band(178535849-178525989)
# Exon ranges: Z-disk(1-28), I-band(29-252), A-band(253-358), M-band(359-364)
TTN_DOMAINS = {
    "Z-disk": {"start": 1, "end": 28, "color":"#FF6B6B"},
    "I-band": {"start": 29, "end": 252, "color": "#4ECDC4"},
    "A-band": {"start": 253, "end": 358, "color": "#45B7D1"},
    "M-band": {"start": 359, "end": 364, "color": "#96CEB4"}
}

# TTN Transcripts
TTN_TRANSCRIPTS = {
    "N2BA": {"id": "ENST00000589042", "length": 34350, "description": "Full-length cardiac isoform"},
    "N2B": {"id": "ENST00000460472", "length": 26926, "description": "Shorter cardiac isoform"},
    "Novex-1": {"id": "ENST00000591111", "length": 5604, "description": "Cardiac isoform"},
    "Novex-2": {"id": "ENST00000342992", "length": 5476, "description": "Cardiac isoform"},
    "Novex-3": {"id": "ENST00000360870", "length": 5270, "description": "Skeletal muscle isoform"}
}

# HTML Report Configuration
REPORT_TITLE = "TTN Variant Analysis Report"
REPORT_HEADER_COLOR = "#2C3E50"
REPORT_ACCENT_COLOR = "#3498DB"

# Phenotype categories
PHENOTYPE_CATEGORIES = {
    "heart": ["cardiomyopathy", "DCM", "HCM", "ARVC", "cardiac"],
    "skeletal_muscle": ["myopathy", "muscular dystrophy", "tibial", "LGMD"],
    "both": ["heart", "skeletal"]
}
