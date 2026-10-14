"""
Configuration file
"""
import os
from pathlib import Path

# 自動載入 .env（如果存在）
def _load_dotenv():
    env_path = Path(__file__).parent / ".env"
    if not env_path.exists():
        return
    with open(env_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            # 只在環境變數尚未設定時才套用（讓 shell 的設定優先）
            if key and value and key not in os.environ:
                os.environ[key] = value

_load_dotenv()

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
# 校準後的閾值（從 ClinVar 評估資料集挑出 Youden's J 最佳化的切點）
# delta_score < PATHOGENIC_THRESHOLD → pathogenic
PATHOGENIC_THRESHOLD = -0.000714 # evo2閾值
REFERENCE_GENOME_PATH = DATA_DIR / "sequence.fasta"
TTN_SEQUENCE_START = 178807423
TTN_SEQUENCE_END = 178525989

# PubMed / ClinVar / NCBI API 認證
# email 是 NCBI 必要欄位；API Key 可提高速率限制（3→10 req/s）
# 在 .env 中設定 PUBMED_EMAIL / PUBMED_API_KEY，或透過 shell 環境變數傳入
PUBMED_EMAIL = os.getenv("PUBMED_EMAIL", "")
_raw_api_key = os.getenv("PUBMED_API_KEY", "")
PUBMED_API_KEY = _raw_api_key if _raw_api_key not in ("", "your_api_key_here") else ""
PUBMED_MAX_RESULTS = 50
ENABLE_FULL_TEXT_FETCH = True
MAX_TEXT_LENGTH = 80000  # 增加到 80K 字元，充分利用 MedGemma 的長 context 能力（最高 128K tokens）
ENABLE_LITVAR_SEARCH = True
# LitVar2 一次最多撈幾篇 publication（API 回傳的 PMIDs 上限）
LITVAR_MAX_RESULTS = 50

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
# LOCAL_LLM_MODEL = "unsloth/medgemma-27b-text-it-GGUF"
# LOCAL_LLM_MODEL = "google/medgemma-27b-it"  # MedGemma-27B instruction-tuned (醫療專用，27B 參數，約需 54GB)
LOCAL_LLM_MODEL = "nvidia/Gemma-4-31B-IT-NVFP4"
# LOCAL_LLM_MODEL = "google/gemma-4-31B-it"
# LOCAL_LLM_MODEL = "unsloth/gemma-4-31B-it-GGUF"
# LOCAL_LLM_MODEL = "meta-llama/Llama-3.2-8B-Instruct"  # Llama 3.2 8B（約需 16GB）
# LOCAL_LLM_MODEL = "google/gemma-3-27b-it"
LOCAL_LLM_TENSOR_PARALLEL = 1  # 單 GPU 使用 1，多 GPU 可增加（例如 2 或 4）
# vLLM 會嘗試佔用「總顯存 × 這個比例」，不足就直接 CUDA out of memory。
# main.py 的前置檢查用同一個數字判斷現在載不載得起來，兩邊必須一致。
LLM_GPU_MEMORY_UTILIZATION = 0.85
LOCAL_LLM_MAX_MODEL_LEN = 32768  # vLLM 最大模型長度（單個 prompt 上限，tokens）
# 注意：n=5 self-consistency 時會 sample 5 次，KV cache 是 5 倍。
# DGX Spark 128GB unified memory，MedGemma-27B BF16 weights ≈ 54GB，
# 剩下 ~50GB 給 KV cache + activations，5 個 sample × 8K tokens ≈ 35K KV，剛好。
LOCAL_LLM_MAX_CONTEXT_LENGTH = 30000  # 準備給 LLM 的文本最大字元數（≈ 7.5K tokens）
ENABLE_LOCAL_CLINICAL_EXTRACTION = True
ENABLE_EASY_PROMPT = False  # 如果為 True，簡化 LLM prompt，不包含文獻內容，並跳過全文下載、Docling 處理、LitVar 搜尋等步驟（快速測試模式）

# === Extraction quality knobs ===
# B1: Self-consistency — 對同一個 prompt sample N 次，逐欄位 majority vote
# 註：成本是 1x 推論時間 × N，但會顯著降低 hallucination
LLM_SELF_CONSISTENCY_N = 5
LLM_SAMPLING_TEMPERATURE = 0.2  # self-consistency 需要一些隨機性才有用
# B3: Guided JSON decoding — 透過 vLLM structured_outputs 強制輸出符合 schema
LLM_USE_GUIDED_JSON = True
# B4: Evidence grounding — 驗證 LLM 回傳的 evidence_sentence 是否真的出現在原文
LLM_VERIFY_EVIDENCE = True
LLM_EVIDENCE_NGRAM_OVERLAP_THRESHOLD = 0.5  # 5-gram token overlap 比例下限

# A: Deterministic alias gating —「變異是否出現在文中」由程式做正規化精確比對決定，
# 不交給 LLM 用「看起來像」去猜。任何別名都比對不到 → 直接 Not specified，
# 並跳過 LLM 與 regex，避免抓到表格內位置相近的其他變異（false positive）。
# 正規化會吃掉 docling/OCR 雜訊（c.107867T>C / c.107867T . C / c.107867T | C
# 都收斂成 `c 107867 t c`），但仍要求別名「完整連續」命中，數字相近的不同變異不會誤判。
LLM_REQUIRE_DETERMINISTIC_ALIAS = True
# 別名要被當成「可用識別碼」需具備的特異性：正規化後須含 ≥ 此位數的連續數字，
# 或為 rsID。避免短數字巧合命中。
LLM_ALIAS_MIN_DIGIT_RUN = 5

# === Track B: Mechanistic / analogy-based inference（無文獻時的推論） ===
# 當文獻完全沒有提到輸入變異時，改走「類比推理」路徑：
# Evo2 致病性 + 結構 band + 變異類型（truncating/missense/splice）+ 從 HGMD
# 檢索同 band、同類型的已知致病變異作為類比證據，由 LLM 生成標示為
# inference（非 citation）的推論說明。
ENABLE_MECHANISTIC_INFERENCE = True
INFERENCE_TOP_K_ANALOGS = 8  # 每次檢索的類比變異數量
# 舊的 HGMD 類比資料庫（含 POS / REF / ALT / DNA(c.) / PROT(p.) / PHEN / phenotype）。
# INFERENCE_HGMD_CSV_PATH 設定後就不再作為類比庫，只剩 variant_db_lookup 用它
# 的 acc_num → REF/ALT 對照來做精確比對。
HGMD_KB_PATH = DATA_DIR / "HGMD_TTN_processed.xlsx"

# --- 推論的資料來源（pipeline 用到的資料都放在 data/）-----------------------
INFERENCE_DATA_DIR = DATA_DIR  # 舊名稱，保留給仍引用它的程式
# 類比／計數用的 HGMD 完整匯出（2,193 筆，含 disease / mutype / rankscore）。
# 這份沒有 REF/ALT，所以引擎的 exact 比對會關閉，改由 Step 5.5 的查表負責。
INFERENCE_HGMD_CSV_PATH = INFERENCE_DATA_DIR / "HGMD_latest.csv"
# cardiodb：364 個外顯子的 Hg38 座標、結構 band、心肌／骨骼肌實測 PSI。
# 設定後取代轉錄本模型成為 exon / band / PSI 來源（涵蓋 364 vs 舊表的 323）。
INFERENCE_CARDIODB_PATH = INFERENCE_DATA_DIR / "cardiodb.xlsx"

# 第 3 層的決策方式：
#   "hierarchical" 階層式收縮 全庫 → band → PSI 類別 → exon（預設）
#   "analogy"      距離加權類比投票（原本的做法，保留可切換）
INFERENCE_PREDICT_MODE = "hierarchical"
# 每一層往上一層收縮的強度（m-estimate 的等效樣本數）。值越大代表越不信任
# 該層自己的計數。1,926 顆變異分到 364 個外顯子、平均 5.3 顆，m_exon=4 讓
# 典型外顯子的估計約一半來自自身、一半來自 PSI 類別的先驗。
INFERENCE_M_EXON = 4.0
INFERENCE_M_BUCKET = 4.0
INFERENCE_M_BAND = 8.0
# TTN 結構 band 的 genomic POS 分箱邊界（GRCh38, chr2；來自 research/pipeline_figures.ipynb §5）
# pd.cut(POS, bins=TTN_BAND_EDGES, labels=TTN_BAND_LABELS) → 低 POS 到高 POS
TTN_BAND_EDGES = [178525989, 178535982, 178618584, 178775356, 178807423]
TTN_BAND_LABELS = ["M-band", "A-band", "I-band", "Z-disk"]
# transcript interval 檔（用來判斷變異影響哪些 isoform）
TRANSCRIPT_INTERVAL_PATH = DATA_DIR / "transcript_interval.xlsx"
# TTN RefSeq 轉錄本外顯子模型：讓輸入可以是 HGVS（c./p./g.）而不只是 chr-pos-ref-alt。
# 由 build_ttn_transcript_model.py 從 NCBI RefSeq GFF3 產生（離線換算，不需 VEP）。
TTN_TRANSCRIPT_MODEL_PATH = DATA_DIR / "ttn_transcript_model.json"
# rsID 是唯一需要連網（dbSNP）的解析路徑；預設關閉，需要時再開
HGVS_ALLOW_NETWORK = False

# === 資料庫查表：這顆變異在 HGMD / PubTator 裡登記的疾病是什麼 ===
# 純查表、不做推論；命中就把資料庫原文列在報告的 Curated Database Annotations 區塊。
ENABLE_DATABASE_ANNOTATIONS = True
# HGMD 完整匯出（含 disease / mutype / inheritance / pmid；REF/ALT 靠 acc_num 接 HGMD_KB_PATH）
HGMD_FULL_CSV_PATH = INFERENCE_DATA_DIR / "HGMD_latest.csv"
# PubTator 文獻探勘出的 TTN 變異。這一版多了 ClinVar 標註欄（position /
# chromosome / clinvar_match_status）：唯一命中的 2,142 列直接有權威座標，
# 其餘仍靠 hgvs_resolver 離線換算。
PUBTATOR_VARIANTS_CSV_PATH = DATA_DIR / "pubtator" / "pubtator_ttn_variants_clinvar_annotated.csv"
# PubTator HGVS → 座標的解析快取（第一次建索引後就不用重算）
PUBTATOR_COORD_CACHE_PATH = OUTPUT_DIR / "pubtator_hgvs_coord_cache.json"
# PubMind 文獻探勘出的 TTN 變異（parsed_variants 為正股 allele，比對時取互補）
PUBMIND_JSON_PATH = DATA_DIR / "pubmind_ttn.json"

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