"""
Configuration file for TTN Variant AI Agent
Contains constants, API keys, and file paths
"""

from pathlib import Path
import os

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
REFERENCE_GENOME_PATH = DATA_DIR / "GCF_000001405.40_GRCh38.p14_genomic.fna"

# PubMed API Configuration
PUBMED_EMAIL = os.getenv("PUBMED_EMAIL", "ryan910702@gmail.com")
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY", '987512860b8bf8e96a09d672b03ff32e2c08')
PUBMED_MAX_RESULTS = 20

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