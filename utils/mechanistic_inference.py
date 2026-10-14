#!/usr/bin/env python3
"""
Mechanistic / Analogy-based Inference (Track B)

當文獻完全沒有提到輸入變異時（ClinVar 沒 cite、或 cite 的文章沒提到該變異），
現有的 evidence-extraction 流程無法產生任何有依據的疾病判斷。

本模組實作教授建議的「類比式推理（analogy-based / case-based reasoning）」：

    Step 1  Evo2 給出致病性分數與信心（沿用既有 evo2_result）
    Step 2  判斷變異落在哪個結構 band / 外顯子、是否為 truncating、影響哪些 isoform
    Step 3  從 HGMD 檢索同一外顯子（不足才放寬到同 band）、同一變異類型的已知致病
            變異作為類比案例
    Step 3b 用該外顯子的 isoform 收錄度（exon-level PSI 的序位代理）算出組織的
            likelihood ratio，替類比投票加權——低收錄度的外顯子替代性剪接多、
            心肌 PSI 低，票數應偏向骨骼肌表型
    Step 4  由 LLM 依這些「間接證據」撰寫標示為 inference（非 citation）的推論說明

Step 3b 的依據見 RELATED_WORK.md §3：TTN 領域用 exon-level PSI 而非 band 判斷表型
（Roberts 2015；Circ Genom Precis Med 2025 量到 106 個「骨骼肌高 PSI／心肌低 PSI」
的外顯子且全部落在 I-band）。真 PSI 需要 cardiodb.org/titin 的組織 RNA-seq 註釋，
目前取不到，因此改用 data/ttn_transcript_model.json 的 6 條真實轉錄本算收錄度代理。
量測結果（HGMD leave-one-out / 獨立的 PubTator 標註）：

    組織準確率  72.9% → 74.2%  /  69.5% → 70.3%
    疾病命中率  67.9% → 68.7%  /  55.3% → 55.9%

注意 transcript_interval.xlsx 只有每條 isoform 2–9 段粗區間，且 N2BA 幾乎涵蓋 N2A
的全部範圍，所以「是否落在心肌 isoform 內」這個二元判斷在整條基因上幾乎恆為真
（實測 363 個 meta 外顯子中只有 1 個是心肌專屬、0 個是骨骼肌專屬）。收錄度必須用
per-exon 模型算，不能用那份區間表。

輸出明確標記為 inference（低證據強度），與 evidence-grounded 的結果在報告上分軌呈現。
技術上是 retrieval + prompting，不需 fine-tune。

輸入格式：除了 `variant_info` dict（含 pos/ref/alt）之外，也接受文獻常見的 HGVS
字串，例如 `c.98671-2A>G`、`NM_001267550.2:c.29800T>C`、`p.R31056*`、
`g.178704672A>G`、`2-178527098-C-T`（見 utils/hgvs_resolver.py）：

    mi.infer(extractor, "c.98671-2A>G", clinvar_info=None, evo2_result=None)
"""

import logging
import re
from bisect import bisect_right
from typing import Any, Dict, List, Optional, Set, Union

logger = logging.getLogger(__name__)


# ── Isoform → 組織類別（論文 §1.5.1 / §1.5.3）─────────────────────────────────
# N2BA / N2B 為成人心肌主要 isoform；N2A 為骨骼肌 isoform（永遠缺 N2B 外顯子）；
# Novex-1/2/3 為次要 isoform；Meta 為 meta-transcript（參考，含幾乎所有外顯子）。
TRANSCRIPT_TISSUE_CLASS = {
    "N2BA": "cardiac",
    "N2B": "cardiac",
    "Cronos": "cardiac",
    "N2A": "skeletal",
    "Nvx1": "minor",
    "Nvx2": "minor",
    "Nvx3": "minor",
    "Novex1": "minor",
    "Novex2": "minor",
    "Novex3": "minor",
    "Meta": "reference",
}
# 論文 §1.5.1：Z-disk / A-band / M-band 外顯子幾乎全為 constitutive（所有 isoform 都保留），
# 唯有 I-band 的 spring 區大量替代性剪接。因此只有 I-band 的 isoform 覆蓋差異能區分組織。
_CONSTITUTIVE_BANDS = {"Z-disk", "A-band", "M-band"}


# ── 外顯子 isoform 收錄度（PSI 代理）─────────────────────────────────────────
# RELATED_WORK.md §3：TTN 領域用 exon-level PSI（percent spliced in）而非 band 判斷
# 表型；Circ Genom Precis Med 2025（PMC11999099）量到心肌有 187 個高 PSI 外顯子、
# 另有 106 個「骨骼肌高 PSI／心肌低 PSI」的外顯子且全部落在 I-band。
#
# 真正的 PSI 需要 cardiodb.org/titin 的組織 RNA-seq 註釋（目前取不到）。這裡用
# data/ttn_transcript_model.json 的 6 條真實 RefSeq 轉錄本算「有幾條收錄該外顯子」
# 當序位代理：收錄的 isoform 越少 = 替代性剪接越多 = 心肌 PSI 越低。
#
# 心肌 / 骨骼肌的歸類不靠 accession 命名，而是由資料判定：meta exon 48（2,781 bp）
# 是 N2B 專屬元件，只有心肌 isoform 收錄它。含此外顯子者即心肌側。
_N2B_UNIQUE_PROBE = 48  # 文獻外顯子編號
CARDIAC_TRANSCRIPTS = (
    "NM_003319.4",      # N2B    26,927 aa
    "NM_133437.4",      # 27,119 aa
    "NM_133432.3",      # 27,052 aa
    "NM_001256850.1",   # N2BA   34,351 aa（長心肌 isoform，含 N2B 專屬外顯子）
)
SKELETAL_TRANSCRIPTS = ("NM_133378.4",)   # N2A 33,424 aa（缺 N2B 專屬外顯子）
MINOR_TRANSCRIPTS = ("NM_133379.5",)      # Novex-3 5,605 aa
# 收錄度分桶。HGMD 上的心肌比例：0 → 18%、1-2 → 48%、3-4 → 86%、5-6 → 72%
# （整體基準 69%），所以低收錄度是明確的骨骼肌訊號。
INCLUSION_BUCKETS = ((0, 0, "meta-only"), (1, 2, "low"), (3, 4, "mid"), (5, 6, "high"))

# 真 PSI（優先於上面的代理）。來源：Circ Genom Precis Med 2025 補充表 S2
# （circcvg-2024-004982-t-s01.pdf），已解析成 data/ttn_exon_psi.csv，欄位為
# exon / domain / psi_heart / psi_skeletal / meta_only。
# 該表用的就是 NM_001267550.2 的外顯子編號（排除非編碼的 exon 1），與本模組的
# 文獻編號對齊：323 個外顯子中有 320 個的 band 判定與 TTN_BAND_EDGES 一致，
# 3 個不一致全部落在 band 邊界（exon 29 / 252 / 358）。
# 解析結果可重現論文的兩個關鍵數字：心肌高 PSI 187 個外顯子、
# 「骨骼肌高 PSI 但心肌低 PSI」106 個且全部位於 I-band。
EXON_PSI_FILENAME = "ttn_exon_psi.csv"
PSI_HIGH = 90  # 論文定義：spliced transcripts 中保留率 >= 90% 為 high PSI


# ── 外顯子層級定位（RELATED_WORK.md：Savarese 2020、Roberts 2015）─────────────
# band 只有 4 格、跨度 28 萬 bp，解析度太粗；文獻談 TTN 熱點一律以外顯子為單位。
#
# 編號校正：data/ttn_transcript_model.json 的 NM_001267550.2 只有 363 個外顯子，
# 而文獻與 RefSeq 註釋用的是 364 個外顯子的編號，兩者差一。用三顆已知變異對過：
#   HMERF p.C31712R (c.95134)、p.N31786K (c.95358) → 模型 exon 343 = 文獻 exon 344
#   TMD FINmaj (c.107780_107790)                   → 模型 exon 363 = 文獻 exon 364
# 對外（features、報告、熱點表）一律用文獻編號。
EXON_TRANSCRIPT = "NM_001267550.2"
LITERATURE_EXON_OFFSET = 1

# 文獻已確立的外顯子級熱點，等同 ACMG 的 PM1（突變熱點 / 關鍵功能域）。
# 這兩個外顯子的表型專一性遠高於所在 band 的整體先驗，值得凌駕距離加權投票。
EXON_HOTSPOTS: Dict[int, Dict[str, Any]] = {
    344: {
        "domain": "FN3-119 (A150)",
        "vgroups": ("missense",),
        "disease": "Myopathy, myofibrillar, 9, with early respiratory failure (HMERF)",
        "tissue": "Skeletal",
        "strength": "strong",  # HGMD 17 顆 missense 中 16 顆為 HMERF
        "note": (
            "Exon 344 encodes the FN3-119 (A150) domain, the established HMERF mutation "
            "hotspot: essentially all reported HMERF variants are missense changes in this "
            "single exon, so exon membership is far more informative than the A-band prior."
        ),
    },
    364: {
        "domain": "Mex6 (C-terminal M-band)",
        # 不限類型：HGMD 在這個外顯子的 9 顆變異裡 8 顆是 missense，和文獻偏重
        # FINmaj indel 的印象相反，限定 truncating 會讓規則幾乎永遠不觸發。
        "vgroups": None,
        "disease": "Tibial muscular dystrophy (TMD) / LGMD2J",
        "tissue": "Skeletal",
        "strength": "moderate",  # HGMD 9 顆中 7 顆為骨骼肌（TMD 4、LGMD 2）
        "note": (
            "Exon 364 (Mex6) is the last TTN exon. Variants here cause tibial muscular "
            "dystrophy when heterozygous and LGMD2J when homozygous — the two are allelic, "
            "so they are reported together. Being in the final exon these alleles escape "
            "nonsense-mediated decay, so the mechanism differs from A-band TTNtv."
        ),
    },
}


# ── 變異類型分類 ─────────────────────────────────────────────────────────────
# vgroup：truncating（截斷型 / LoF：nonsense、frameshift、splice-loss）
#         vs missense（非截斷型），對齊教授說的「truncating / missense」
# vclass：更細的類別，供報告顯示與 prompt 說明


def evo2_confidence_bucket(delta_score) -> str:
    """把 Evo2 delta score 轉成論文 §2.2.3 定義的三級信心桶。

    High:   |delta| > 5e-3
    Medium: 1e-3 <= |delta| <= 5e-3
    Low:    otherwise
    """
    try:
        a = abs(float(delta_score))
    except (TypeError, ValueError):
        return "Unknown"
    if a > 5e-3:
        return "High"
    if a >= 1e-3:
        return "Medium"
    return "Low"


def _classify_variant(prot_hgvs: str, nuc_hgvs: str) -> Dict[str, str]:
    """
    根據蛋白質 / 核苷酸 HGVS 判斷變異類型。

    prot_hgvs 例：NP_596869.4:p.Q33396*、p.Ile35947Asn、p.Arg23029ThrfsX9
    nuc_hgvs  例：NM_133378.4:c.100186C>T、c.98671-2A>G（splice）
    """
    p = (prot_hgvs or "").strip()
    d = (nuc_hgvs or "").strip()

    # splice：核苷酸座標帶 ±1 / ±2 內含子位置（c.####+1 / c.####-2 …）
    if re.search(r"\d+[+-][12](?:[ACGTacgt]|del|ins|dup|>)", d):
        return {"vclass": "splice", "vgroup": "truncating"}

    # frameshift
    if re.search(r"fs", p, re.IGNORECASE):
        return {"vclass": "frameshift", "vgroup": "truncating"}

    # nonsense / stop-gain：p.Xxx###Ter、p.X###*、單字母 *、或末端 X
    p_core = re.sub(r"^[A-Za-z0-9_.]+:", "", p)  # 去掉 NP_xxx: 前綴
    if (
        "Ter" in p_core
        or p_core.endswith("*")
        or re.search(r"\d+\*", p_core)
        or re.search(r"[A-Za-z]\d+X$", p_core)
    ):
        return {"vclass": "nonsense", "vgroup": "truncating"}

    # missense：p.Xxx###Yyy（三字母）或 p.X###Y（單字母），且結尾不是終止
    if re.search(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}$", p_core) or re.search(
        r"p\.[A-Za-z]\d+[A-Za-z]$", p_core
    ):
        return {"vclass": "missense", "vgroup": "missense"}

    return {"vclass": "unknown", "vgroup": "unknown"}


class MechanisticInference:
    """類比式推理引擎（Track B）。"""

    def __init__(
        self,
        kb_path: str,
        band_edges: List[int],
        band_labels: List[str],
        top_k: int = 8,
        transcript_interval_path: Optional[str] = None,
        transcript_model_path: Optional[str] = None,
        reference_fasta_path: Optional[str] = None,
        allow_network: bool = False,
        distance_scale: float = 150.0,
        distance_weight_mode: str = "exponential",
        retrieval_mode: str = "exon",
        exon_hotspot_rules: bool = True,
        tissue_prior: Optional[Dict[str, float]] = None,
        rankscore_weighting: float = 0.5,
        vague_penalty: float = 0.6,
        inclusion_tissue_weighting: bool = True,
        inclusion_lr_strength: float = 2.0,
        inclusion_lr_cap: float = 4.0,
        inclusion_min_n: int = 15,
        exon_psi_path: Optional[str] = None,
        predict_mode: str = "analogy",
        m_exon: float = 4.0,
        m_bucket: float = 4.0,
        m_band: float = 8.0,
        cardiodb_path: Optional[str] = None,
        hgmd_csv_path: Optional[str] = None,
    ):
        self.kb_path = str(kb_path)
        self.band_edges = list(band_edges)
        self.band_labels = list(band_labels)
        self.top_k = int(top_k)
        # 規則式投票時「距離 → 權重」的形狀。預設為每 distance_scale bp 減半；
        # 這兩個參數可調，是為了讓 research/test_inference_strategies.ipynb 掃描不同策略。
        self.distance_scale = float(distance_scale)
        self.distance_weight_mode = str(distance_weight_mode)
        # "band"：只用 4 格結構 band 過濾類比（改動前的行為，保留供對照）
        # "exon"：先找同一外顯子的類比，不足才依序放寬到同 band（見 retrieve_analogs）
        # 預設用 exon：在證據等級清理過的 PubTator 跨來源測試上四項指標全正
        # （hit +0.4%、macro +0.5%、tissue +0.7%）。幅度小，但跨來源才是「遇到未研究
        # 變異」這個實際情境的代理指標；HGMD 同源測試的 -0.1% 在雜訊範圍內。
        self.retrieval_mode = str(retrieval_mode)
        self.exon_hotspot_rules = bool(exon_hotspot_rules)
        # 目標族群的組織基準率（例：{"cardiac": 0.79, "skeletal": 0.17}）。
        # 類比投票估的是「HGMD 這個抽樣分佈下」的機率，但 HGMD 收錄偏誤讓骨骼肌
        # 變異群聚得特別緊（HMERF 全擠在 exon 344、TMD 在 364），局部投票會被它們
        # 洗掉。設了這個之後，每票會再乘上 目標基準率 / 類比庫基準率 做校正。
        self.tissue_prior = dict(tissue_prior) if tissue_prior else None
        self._tissue_correction = None
        # 讓類比的票數依 HGMD 自己的標註信心（RANKSCORE）打折。0=不打折，
        # 1=直接乘上 rankscore，中間值是線性內插：(1-s) + s * rankscore。
        # 低分項目的標註明顯較不可靠（LOO 命中率 53.5% vs 高分的 77.5%），
        # 不該和高分項目等值投票，但整個歸零也會丟掉它們僅有的資訊。
        #
        # 預設 0.5：跨來源測試四項指標全正（hit +1.0%、macro +1.7%、tissue +1.2%），
        # 而 HGMD 同源只掉 0.5% hit、macro 反而是所有組態最高。強度拉到 1.0 跨來源
        # 還會更好（hit +1.6%），但同源會掉 1.5%——那個退步即使只看最可靠的答案也
        # 存在，屬於同源過擬合（低分類比共用 HGMD 的標註慣例，對外不迴化），所以
        # 取折衷值而非全開。
        self.rankscore_weighting = float(rankscore_weighting)
        # 疾病名籠統（Cardiomyopathy / Titinopathy / Myopathy …）的類比投票折扣。
        # HGMD 約三成的 PHEN 是這類；它們拿滿票時會把預測也拉成籠統答案，
        # 雖然算「相容」但沒有臨床價值。1.0=不打折。
        #
        # 預設 0.6：兩套 ground truth 上 specific 都 +3 個百分點左右、hit 掉約 2。
        # 這是刻意的取捨——Track B 的任務定義就是「必須給出單一最可能的疾病」，
        # 回答 Titinopathy 雖然算相容卻等於沒回答，hit 有一部分是靠模稜兩可撐的。
        # 折扣再往下拉（0.15）specific 還能多 4 點，但具體而錯誤的診斷臨床風險較高，
        # 所以停在溫和的一檔。
        self.vague_penalty = float(vague_penalty)
        # 外顯子 isoform 收錄度（PSI 代理）→ 組織 likelihood ratio。
        # strength 是 LR 的指數（0 = 關閉、1 = 原始 LR）；cap 夾住極端值；
        # min_n 是桶內最少樣本數，低於就不表態。
        self.inclusion_tissue_weighting = bool(inclusion_tissue_weighting)
        self.inclusion_lr_strength = float(inclusion_lr_strength)
        self.inclusion_lr_cap = float(inclusion_lr_cap)
        self.inclusion_min_n = int(inclusion_min_n)
        # 實測 PSI 表；None 時用 transcript_model_path 同目錄的 ttn_exon_psi.csv
        self.exon_psi_path = str(exon_psi_path) if exon_psi_path else None
        # 第 3 層的決策方式：
        #   "analogy"      距離加權類比投票（原本的做法）
        #   "hierarchical" 階層式收縮 exon → PSI 桶 → band → 全庫，不看距離
        # m_* 是每一層往上一層收縮的強度（m-estimate 的等效樣本數），
        # 值越大代表越不信任該層自己的計數。預設值由 PubTator 掃描選出。
        self.predict_mode = str(predict_mode)
        self.m_exon = float(m_exon)
        self.m_bucket = float(m_bucket)
        self.m_band = float(m_band)
        self._hier_cache = None
        # cardiodb.xlsx：364 個外顯子的 Hg38 座標、結構 band、心肌／骨骼肌實測 PSI。
        # 設定後即取代轉錄本模型成為 exon / band / PSI 的來源（涵蓋率 364 vs 323，
        # 且 band 直接來自該表而非座標切點）。未設定時完全維持原本行為。
        self.cardiodb_path = str(cardiodb_path) if cardiodb_path else None
        # HGMD 完整匯出 CSV（2,193 筆），取代舊的 xlsx 知識庫作為類比／計數來源。
        # 這份沒有 REF/ALT 欄，所以 hgmd_lookup 的 exact 比對會關閉，改由
        # main.py Step 5.5 的 variant_db_lookup 負責（那裡靠 acc_num 接回 REF/ALT）。
        self.hgmd_csv_path = str(hgmd_csv_path) if hgmd_csv_path else None
        self._cardio_rows = None   # list[dict]，依 lo 遞增
        self._cardio_band = None   # dict[exon -> band]
        self.transcript_interval_path = (
            str(transcript_interval_path) if transcript_interval_path else None
        )
        # HGVS → genomic 座標所需（皆可為 None，屆時用 hgvs_resolver 的預設路徑）
        self.transcript_model_path = transcript_model_path
        self.reference_fasta_path = reference_fasta_path
        self.allow_network = bool(allow_network)  # 只影響 rsID 查 dbSNP

        self._kb = None  # list[dict]，延遲載入
        self._transcripts = None  # dict[name -> list[(lo, hi)]]
        self._resolver = None  # HgvsResolver，延遲載入
        self._exons = None  # list[(g_start, g_end, literature_exon_no)]，延遲載入
        self._exon_starts = None  # 供 bisect 用的 g_start 遞增列表
        self._exon_inclusion = None  # dict[文獻外顯子編號 -> 收錄度資訊]
        self._exon_psi = None  # dict[外顯子 -> 實測 PSI]，延遲載入
        self._inclusion_tally = None  # dict[bucket -> Counter(tissue)]，供 LR 估計

    # ── 載入 HGMD 類比庫 ──────────────────────────────────────────────────────

    # HGMD mutype → (vclass, vgroup)。直接用 HGMD 自己的分類，比正則可靠。
    _MUTYPE_MAP = {
        "frameshift": ("frameshift", "truncating"),
        "nonsense": ("nonsense", "truncating"),
        "canonical-splice": ("splice", "truncating"),
        "exonic-splice": ("splice", "truncating"),
        "gross": ("gross", "truncating"),
        "missense": ("missense", "missense"),
        "inframe": ("inframe", "missense"),
        "noncoding": ("noncoding", "unknown"),
    }

    def _load_kb_from_csv(self) -> List[Dict[str, Any]]:
        """解析 HGMD 完整匯出 CSV（含 disease / mutype / rankscore / clinvarID）。

        沒有 REF/ALT 欄，所以 ref/alt 留空，比對一律走同位置（見 hgmd_lookup）。
        """
        import pandas as pd

        from utils.disease_normalizer import (
            canonical_labels,
            is_vague_disease,
            tissue_of_labels,
        )

        df = pd.read_csv(self.hgmd_csv_path, low_memory=False)
        rows: List[Dict[str, Any]] = []
        for r in df.itertuples():
            try:
                pos = int(r.startCoord)
            except (TypeError, ValueError):
                continue
            disease = str(getattr(r, "disease", "") or "").strip()
            vclass, vgroup = self._MUTYPE_MAP.get(
                str(getattr(r, "mutype", "") or "").strip(), ("unknown", "unknown")
            )
            rank = getattr(r, "rankscore", None)
            cvid = getattr(r, "clinvarID", None)
            rows.append({
                "pos": pos,
                "ref": "", "alt": "",
                "dna": str(getattr(r, "hgvs", "") or ""),
                "prot": str(getattr(r, "amino", "") or ""),
                "phen": disease,
                # 引擎內部用小寫組織標籤；tissue_of_labels 回傳 "Cardiac" 這種首字大寫
                "phenotype": tissue_of_labels(canonical_labels(disease)).lower(),
                "rankscore": None if rank is None or pd.isna(rank) else float(rank),
                "rsid": str(getattr(r, "dbsnp", "") or ""),
                "clinvar_id": None if cvid is None or pd.isna(cvid) else str(cvid),
                "phen_vague": is_vague_disease(disease),
                "band": self.assign_band(pos),
                "exon": self.assign_exon(pos),
                "vclass": vclass, "vgroup": vgroup,
            })
        logger.info(f"Loaded {len(rows)} HGMD variants from {self.hgmd_csv_path}")
        return rows

    def _load_kb(self) -> List[Dict[str, Any]]:
        if self._kb is not None:
            return self._kb

        import pandas as pd

        from utils.disease_normalizer import is_vague_disease

        if self.hgmd_csv_path:
            self._kb = self._load_kb_from_csv()
            return self._kb

        df = pd.read_excel(self.kb_path)
        rows: List[Dict[str, Any]] = []
        for _, r in df.iterrows():
            try:
                pos = int(r["POS"])
            except (ValueError, TypeError):
                continue
            prot = str(r.get("PROT") or "")
            dna = str(r.get("DNA") or "")
            cls = _classify_variant(prot, dna)
            rows.append(
                {
                    "pos": pos,
                    "ref": str(r.get("REF") or ""),
                    "alt": str(r.get("ALT") or ""),
                    "dna": dna,
                    "prot": prot,
                    "phen": str(r.get("PHEN") or "").strip(),
                    # HGMD 自己的標註信心（0~1）。實測它強烈預測標註可靠度：
                    # leave-one-out 命中率從 <0.3 的 53.5% 一路到 >0.95 的 77.5%。
                    "rankscore": (
                        float(r["RANKSCORE"])
                        if "RANKSCORE" in df.columns and not pd.isna(r.get("RANKSCORE"))
                        else None
                    ),
                    "phenotype": str(r.get("phenotype") or "").strip().lower(),
                    "rsid": str(r.get("DB") or "").strip(),
                    "clinvar_id": (
                        str(r.get("ClinVar_Variation_ID"))
                        if "ClinVar_Variation_ID" in df.columns
                        and not pd.isna(r.get("ClinVar_Variation_ID"))
                        else None
                    ),
                    "phen_vague": is_vague_disease(r.get("PHEN")),
                    "band": self.assign_band(pos),
                    "exon": self.assign_exon(pos),
                    "vclass": cls["vclass"],
                    "vgroup": cls["vgroup"],
                }
            )
        self._kb = rows
        logger.info(f"Loaded HGMD analogy KB: {len(rows)} variants from {self.kb_path}")
        return rows

    # ── HGVS 輸入 ────────────────────────────────────────────────────────────

    @property
    def resolver(self):
        """HGVS 解析器（延遲建立；會用 HGMD 的 c.HGVS 當作精確索引）。"""
        if self._resolver is None:
            from utils.hgvs_resolver import HgvsResolver, build_hgmd_hgvs_index

            try:
                index = build_hgmd_hgvs_index(self._load_kb())
            except Exception as e:  # noqa: BLE001
                logger.warning(f"無法由 HGMD 建立 HGVS 索引：{e}")
                index = {}
            self._resolver = HgvsResolver(
                model_path=self.transcript_model_path,
                fasta_path=self.reference_fasta_path,
                hgmd_index=index,
                allow_network=self.allow_network,
            )
        return self._resolver

    def resolve_variant(self, query: str, hgvs_p: Optional[str] = None) -> Dict[str, Any]:
        """
        把 HGVS 字串解析成 pipeline 慣用的 variant_info dict。

        query 可為：c.98671-2A>G / NM_001267550.2:c.29800T>C / p.R31056* /
        g.178704672A>G / 2-178527098-C-T / rs794729424（需 allow_network）。
        解析不出來會丟 HgvsResolveError。
        """
        return self.resolver.resolve(str(query), hgvs_p=hgvs_p).to_variant_info()

    def normalize_variant_input(
        self, variant: Union[str, Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        統一 infer() / build_features() 的輸入：接受 HGVS 字串或 variant_info dict。

        dict 若已有 pos 就直接用；只給 hgvs / hgvs_c / hgvs_p / variant_id 時，
        會用 HGVS 解析器補上 pos / ref / alt / 變異類型。
        """
        if isinstance(variant, str):
            return self.resolve_variant(variant)

        info = dict(variant or {})
        if info.get("pos") not in (None, ""):
            return info

        prot = info.get("hgvs_p") or info.get("hgvs_protein") or info.get("protein_change")
        candidates = [
            info.get("hgvs_c"), info.get("hgvs"), info.get("hgvs_nucleotide"),
            prot, info.get("variant_id"), info.get("rsid"),
        ]
        flat: List[str] = []
        for c in candidates:
            if isinstance(c, (list, tuple)):
                flat.extend(str(x) for x in c if x)
            elif c not in (None, ""):
                flat.append(str(c))
        resolved = self.resolver.resolve_best(flat)
        if resolved is None:
            raise ValueError(
                "variant_info 沒有 pos，也無法從 HGVS 解析出座標："
                f"{[c for c in flat][:4]}"
            )
        merged = resolved.to_variant_info()
        # 呼叫端原本給的欄位（例如 disease 標註）保留，但座標以解析結果為準
        for k, v in info.items():
            merged.setdefault(k, v)
        return merged

    # ── Step 2: band / 類型 / isoform ────────────────────────────────────────

    def assign_band(self, pos: int) -> str:
        """把 genomic POS 分箱到結構 band。

        有 cardiodb 時以該表登記的 band 為準（外顯子層級）；落在內含子或
        表外的位置退回 band_edges 的座標切點。
        """
        if self.cardiodb_path:
            self._load_cardio()
            exon = self.assign_exon(pos)
            if exon is not None:
                band = (self._cardio_band or {}).get(int(exon))
                if band:
                    return band
        # band_edges 為遞增邊界；label i 對應區間 [edges[i], edges[i+1]]
        idx = bisect_right(self.band_edges, pos) - 1
        if idx < 0 or idx >= len(self.band_labels):
            return "outside"
        return self.band_labels[idx]

    def _load_cardio(self) -> List[Dict[str, Any]]:
        """讀 cardiodb.xlsx：每個外顯子的 Hg38 區間、結構 band、兩個組織的 PSI。

        回傳依 lo 遞增排序的 list[dict]；未設定路徑或讀取失敗時回傳空 list，
        呼叫端一律退回轉錄本模型的原本行為。
        """
        if self._cardio_rows is not None:
            return self._cardio_rows

        self._cardio_rows = []
        self._cardio_band = {}
        if not self.cardiodb_path:
            return self._cardio_rows
        try:
            import os

            import pandas as pd

            if not os.path.exists(self.cardiodb_path):
                logger.warning(f"找不到 cardiodb，改用轉錄本模型：{self.cardiodb_path}")
                return self._cardio_rows

            # 第 0 列是 isoform 名稱標頭，不是資料。欄名含空白，itertuples 會把它們
            # 改寫成 _1.._16 而容易對錯位置，所以用 iterrows + 欄名存取。
            df = pd.read_excel(self.cardiodb_path).iloc[1:]
            rows: List[Dict[str, Any]] = []
            for _, r in df.iterrows():
                try:
                    ex = int(r["Exon number"])
                    a, b = int(r["Hg38 start"]), int(r["Hg38 end"])
                except (TypeError, ValueError, KeyError):
                    continue
                # cardiodb 對 36 個外顯子沒有骨骼肌 PSI（exon 148/150/172-205，全在
                # I-band 且心肌 PSI 也都 <90%）。缺值視為 0，它們因此落入 both-low，
                # 語意上正確：這些是兩個組織都幾乎不收錄的替代性剪接外顯子。
                ph, ps = r.get("PSI in heart"), r.get("PSI in skeletal")
                rows.append({
                    "exon": ex,
                    "lo": min(a, b), "hi": max(a, b),
                    "band": str(r.get("Band") or "").strip(),
                    "psi_heart": 0.0 if pd.isna(ph) else float(ph),
                    "psi_skeletal": 0.0 if pd.isna(ps) else float(ps),
                })
            rows.sort(key=lambda x: x["lo"])
            self._cardio_rows = rows
            self._cardio_band = {r["exon"]: r["band"] for r in rows if r["band"]}
            logger.info(f"Loaded {len(rows)} exons from cardiodb ({self.cardiodb_path})")
        except Exception as e:  # noqa: BLE001
            logger.warning(f"讀取 cardiodb 失敗，改用轉錄本模型：{e}")
            self._cardio_rows = []
            self._cardio_band = {}
        return self._cardio_rows

    def _load_exons(self) -> List[tuple]:
        """載入 meta transcript 的外顯子區間，回傳 (g_start, g_end, 文獻外顯子編號)。

        依 g_start 遞增排序（TTN 在負鏈，所以外顯子編號是遞減的），供 bisect 查找。
        """
        if self._exons is not None:
            return self._exons

        # cardiodb 優先：涵蓋 364 個外顯子且座標直接來自該表
        cardio = self._load_cardio()
        if cardio:
            spans = [(r["lo"], r["hi"], r["exon"]) for r in cardio]
            self._exons = spans
            self._exon_starts = [s[0] for s in spans]
            return self._exons

        self._exons = []
        self._exon_starts = []
        if not self.transcript_model_path:
            return self._exons
        try:
            import json
            import os

            if not os.path.exists(str(self.transcript_model_path)):
                logger.warning(f"找不到轉錄本模型，exon 檢索停用：{self.transcript_model_path}")
                return self._exons
            with open(str(self.transcript_model_path), encoding="utf-8") as fh:
                model = json.load(fh)
            tx = (model.get("transcripts") or {}).get(EXON_TRANSCRIPT)
            if not tx:
                logger.warning(f"轉錄本模型缺少 {EXON_TRANSCRIPT}，exon 檢索停用")
                return self._exons

            # exons 依轉錄本順序排列，索引 i 即模型編號 i+1；文獻編號再加偏移。
            spans = [
                (int(e["g_start"]), int(e["g_end"]), i + 1 + LITERATURE_EXON_OFFSET)
                for i, e in enumerate(tx["exons"])
            ]
            spans.sort(key=lambda s: s[0])
            self._exons = spans
            self._exon_starts = [s[0] for s in spans]
            logger.info(f"Loaded {len(spans)} exons from {EXON_TRANSCRIPT} for exon-level retrieval")
        except Exception as e:  # noqa: BLE001
            logger.warning(f"載入外顯子模型失敗，exon 檢索停用：{e}")
            self._exons = []
            self._exon_starts = []
        return self._exons

    def assign_exon(self, pos: int) -> Optional[int]:
        """回傳位置所在的外顯子編號（文獻編號）；落在內含子或基因外則為 None。"""
        exons = self._load_exons()
        if not exons:
            return None
        try:
            p = int(pos)
        except (TypeError, ValueError):
            return None
        # 找最後一個 g_start <= p 的外顯子；外顯子不重疊，所以只需檢查這一個
        idx = bisect_right(self._exon_starts, p) - 1
        if idx < 0:
            return None
        g_start, g_end, exon_no = exons[idx]
        return exon_no if g_start <= p <= g_end else None

    def exon_hotspot(self, exon: Optional[int], vgroup: str) -> Optional[Dict[str, Any]]:
        """若變異落在文獻確立的外顯子熱點、且變異類型相符，回傳該熱點的註記。"""
        if not self.exon_hotspot_rules or exon is None:
            return None
        hs = EXON_HOTSPOTS.get(int(exon))
        if not hs:
            return None
        # 類型不符時不套用：例如 exon 344 的截斷變異走的是一般 A-band TTNtv 機制，
        # 不是 HMERF（HMERF 是 FN3-119 域的 missense 專屬表型；HGMD 在該外顯子唯一的
        # truncating 變異確實註記為 DCM）。vgroups=None 表示該熱點不分類型。
        vgroups = hs.get("vgroups")
        if vgroups and vgroup != "unknown" and vgroup not in vgroups:
            return None
        return dict(hs, exon=int(exon))

    # ── 外顯子 isoform 收錄度（PSI 代理）──────────────────────────────────────

    def _load_exon_inclusion(self) -> Dict[int, Dict[str, Any]]:
        """算出每個 meta 外顯子被幾條真實 isoform 收錄，當 exon-level PSI 的序位代理。

        回傳 {文獻外顯子編號: {"n_inc", "cardiac", "skeletal", "bucket"}}。
        取不到轉錄本模型時回傳空 dict，呼叫端一律退回原本的 band 行為。
        """
        if self._exon_inclusion is not None:
            return self._exon_inclusion

        # cardiodb 優先：364 個外顯子全部有實測 PSI，不需要 isoform 收錄度代理，
        # 也就沒有 psi:absent 這個桶。
        cardio = self._load_cardio()
        if cardio:
            self._exon_inclusion = {
                r["exon"]: {
                    "psi_heart": r["psi_heart"],
                    "psi_skeletal": r["psi_skeletal"],
                    "band": r["band"],
                    "meta_only": 0,
                    "n_inc": None,
                    "bucket": "psi:" + self._psi_bucket(r["psi_heart"], r["psi_skeletal"], 0),
                    "source": "cardiodb",
                }
                for r in cardio
            }
            return self._exon_inclusion

        self._exon_inclusion = {}
        exons = self._load_exons()  # [(g_start, g_end, 文獻編號)]，已依 g_start 排序
        if not exons or not self.transcript_model_path:
            return self._exon_inclusion

        # 真 PSI 優先。補充表只涵蓋 323 個外顯子（排除非編碼與不在可報告範圍者），
        # 其餘外顯子仍走下面的轉錄本收錄度代理，兩種桶名分開命名（psi: / proxy:），
        # LR 各自估計、互不混用。
        psi_table = self._load_exon_psi()
        try:
            import json
            import os

            if not os.path.exists(str(self.transcript_model_path)):
                return self._exon_inclusion
            with open(str(self.transcript_model_path), encoding="utf-8") as fh:
                tx = (json.load(fh).get("transcripts") or {})

            others = [t for t in tx if t != EXON_TRANSCRIPT]
            spans = {
                t: sorted((int(e["g_start"]), int(e["g_end"])) for e in tx[t]["exons"])
                for t in others
            }

            def covered(t: str, lo: int, hi: int) -> bool:
                # 重疊達該 meta 外顯子一半以上才算收錄（容忍邊界微調）
                for a, b in spans[t]:
                    if a > hi:
                        break
                    ov = min(hi, b) - max(lo, a) + 1
                    if ov > 0 and ov >= 0.5 * (hi - lo + 1):
                        return True
                return False

            for lo, hi, exon_no in exons:
                pr = psi_table.get(int(exon_no))
                if pr:
                    self._exon_inclusion[int(exon_no)] = {
                        "psi_heart": pr["psi_heart"],
                        "psi_skeletal": pr["psi_skeletal"],
                        "meta_only": pr["meta_only"],
                        "domain": pr["domain"],
                        "n_inc": None,
                        "bucket": "psi:" + self._psi_bucket(
                            pr["psi_heart"], pr["psi_skeletal"], pr["meta_only"]
                        ),
                        "source": "psi",
                    }
                    continue
                card = sum(1 for t in CARDIAC_TRANSCRIPTS if t in spans and covered(t, lo, hi))
                skel = sum(1 for t in SKELETAL_TRANSCRIPTS if t in spans and covered(t, lo, hi))
                minor = sum(1 for t in MINOR_TRANSCRIPTS if t in spans and covered(t, lo, hi))
                n_inc = card + skel + minor
                self._exon_inclusion[int(exon_no)] = {
                    "n_inc": n_inc,
                    "cardiac": card,
                    "skeletal": skel,
                    # 補充表刻意排除的外顯子（exon 45-46/147/149/164/172-201，
                    # 皆為 I-band 的替代性剪接區）自成一類：分開成三個代理桶會讓每桶
                    # 都掉到 min_n 以下而全部棄權，反而丟掉訊號。
                    "bucket": "psi:absent",
                    "source": "proxy",
                }
            n_psi = sum(1 for v in self._exon_inclusion.values() if v.get("source") == "psi")
            logger.info(
                f"Exon usage resolved for {len(self._exon_inclusion)} exons "
                f"({n_psi} from measured PSI, {len(self._exon_inclusion) - n_psi} from "
                f"transcript-inclusion proxy)"
            )
        except Exception as e:  # noqa: BLE001
            logger.warning(f"外顯子收錄度計算失敗，PSI 代理停用：{e}")
            self._exon_inclusion = {}
        return self._exon_inclusion

    def _load_exon_psi(self) -> Dict[int, Dict[str, Any]]:
        """讀 data/ttn_exon_psi.csv（論文補充表 S2 解析結果）。取不到就回空 dict。"""
        if self._exon_psi is not None:
            return self._exon_psi
        self._exon_psi = {}
        try:
            import os

            path = self.exon_psi_path
            if not path:
                base = os.path.dirname(str(self.transcript_model_path or ""))
                path = os.path.join(base, EXON_PSI_FILENAME) if base else None
            if not path or not os.path.exists(str(path)):
                logger.info("找不到 ttn_exon_psi.csv，改用轉錄本收錄度代理")
                return self._exon_psi

            import pandas as pd

            df = pd.read_csv(str(path))
            for r in df.itertuples():
                self._exon_psi[int(r.exon)] = {
                    "psi_heart": float(r.psi_heart),
                    "psi_skeletal": float(r.psi_skeletal),
                    "meta_only": int(getattr(r, "meta_only", 0) or 0),
                    "domain": str(getattr(r, "domain", "") or ""),
                }
            logger.info(f"Loaded measured exon PSI for {len(self._exon_psi)} exons from {path}")
        except Exception as e:  # noqa: BLE001
            logger.warning(f"讀取外顯子 PSI 失敗，改用代理：{e}")
            self._exon_psi = {}
        return self._exon_psi

    @staticmethod
    def _psi_bucket(psi_heart: float, psi_skeletal: float, meta_only: int) -> str:
        """依論文自己的分類把外顯子分桶。

        skeletal-high 就是論文點名的那 106 個外顯子（骨骼肌高 PSI、心肌低 PSI，
        且全部落在 I-band）——心肌／骨骼肌的判別力集中在這一類。
        """
        if meta_only:
            return "meta-only"
        h, s = psi_heart >= PSI_HIGH, psi_skeletal >= PSI_HIGH
        if h and s:
            return "both-high"
        if s:
            return "skeletal-high"
        if h:
            return "cardiac-high"
        return "both-low"

    def exon_inclusion(self, pos: int) -> Optional[Dict[str, Any]]:
        """回傳該位置所在外顯子的 isoform 收錄度；落在內含子或無模型時為 None。"""
        exon = self.assign_exon(pos)
        if exon is None:
            return None
        info = self._load_exon_inclusion().get(int(exon))
        return dict(info, exon=int(exon)) if info else None

    def _load_inclusion_tally(self):
        """統計每個收錄度分桶的組織分佈，並保留各外顯子的分項以便 leave-one-out。"""
        if self._inclusion_tally is not None:
            return self._inclusion_tally
        from collections import Counter, defaultdict

        by_bucket: Dict[str, Counter] = defaultdict(Counter)
        by_exon: Dict[int, Counter] = defaultdict(Counter)
        exon_bucket: Dict[int, str] = {}
        incl = self._load_exon_inclusion()
        for r in self._load_kb():
            exon = self.assign_exon(r["pos"])
            info = incl.get(int(exon)) if exon is not None else None
            tissue = r.get("phenotype")
            if not info or tissue not in ("cardiac", "skeletal", "both"):
                continue
            b = info["bucket"]
            by_bucket[b][tissue] += 1
            by_exon[int(exon)][tissue] += 1
            exon_bucket[int(exon)] = b
        self._inclusion_tally = (dict(by_bucket), dict(by_exon), exon_bucket)
        return self._inclusion_tally

    def inclusion_tissue_lr(self, exon: Optional[int]) -> Dict[str, float]:
        """該外顯子收錄度桶對「組織」的 likelihood ratio。

        形式沿用 RELATED_WORK.md §5 建議的校準方式（PM5 pLR 8.4–16.3、同源位置
        pLR 13.0）：回報 LR 而非 High/Low 標籤。

        為避免 RELATED_WORK.md §8 的 circularity，估計時會扣掉**同一外顯子**的所有
        HGMD 變異——查詢變異自己與它的外顯子鄰居都不參與自己的先驗。
        """
        neutral = {"cardiac": 1.0, "skeletal": 1.0, "both": 1.0}
        if exon is None or not self.inclusion_tissue_weighting:
            return neutral
        incl = self._load_exon_inclusion().get(int(exon))
        if not incl:
            return neutral
        by_bucket, by_exon, _ = self._load_inclusion_tally()
        bucket = incl["bucket"]

        own = by_exon.get(int(exon), {})
        in_b = {t: by_bucket.get(bucket, {}).get(t, 0) - own.get(t, 0)
                for t in ("cardiac", "skeletal", "both")}
        overall = {
            t: sum(b.get(t, 0) for b in by_bucket.values()) - own.get(t, 0)
            for t in ("cardiac", "skeletal", "both")
        }
        n_b, n_all = sum(in_b.values()), sum(overall.values())
        # 樣本太少就不表態（避免用 8 顆變異推出 LR=∞）
        if n_b < self.inclusion_min_n or n_all <= 0:
            return neutral

        out = {}
        alpha = self.inclusion_lr_strength
        for t in ("cardiac", "skeletal", "both"):
            # Laplace 平滑，桶內與整體都加 1
            p_b = (in_b[t] + 1.0) / (n_b + 3.0)
            p_all = (overall[t] + 1.0) / (n_all + 3.0)
            lr = (p_b / p_all) if p_all > 0 else 1.0
            lr = lr ** alpha
            out[t] = min(max(lr, 1.0 / self.inclusion_lr_cap), self.inclusion_lr_cap)
        return out

    def classify_input_variant(
        self, variant_info: Dict, clinvar_info: Optional[Dict]
    ) -> Dict[str, str]:
        """
        由 ClinVar 回傳的 HGVS（或輸入自帶的 HGVS）判斷輸入變異的類型，
        全程離線、不呼叫 VEP。若完全沒有 HGVS（真正 novel），回傳 unknown。
        """
        # 1) HGVS 解析器已經算出類型（含離線 CDS 轉譯的 missense/nonsense）→ 直接用
        hv = (variant_info or {}).get("hgvs_vclass")
        if hv and hv != "unknown":
            return {"vclass": hv, "vgroup": (variant_info or {}).get("hgvs_vgroup") or "unknown"}

        prot_candidates: List[str] = []
        nuc_candidates: List[str] = []
        for key in ("hgvs_p",):
            v = (variant_info or {}).get(key)
            if v:
                prot_candidates.append(str(v))
        for key in ("hgvs_c", "hgvs"):
            v = (variant_info or {}).get(key)
            if v:
                nuc_candidates.append(str(v))
        if clinvar_info:
            for key in ("hgvs_protein", "protein_change"):
                v = clinvar_info.get(key)
                if isinstance(v, list):
                    prot_candidates.extend(str(x) for x in v)
                elif v:
                    prot_candidates.append(str(v))
            for key in ("hgvs_nucleotide", "hgvs"):
                v = clinvar_info.get(key)
                if isinstance(v, list):
                    nuc_candidates.extend(str(x) for x in v)
                elif v:
                    nuc_candidates.append(str(v))

        # 依序嘗試每個候選，取第一個非 unknown 的結果
        best = {"vclass": "unknown", "vgroup": "unknown"}
        for prot in prot_candidates or [""]:
            for nuc in nuc_candidates or [""]:
                cls = _classify_variant(prot, nuc)
                if cls["vgroup"] != "unknown":
                    return cls
        # 只有核苷酸（例如純 splice）也要能判
        for nuc in nuc_candidates:
            cls = _classify_variant("", nuc)
            if cls["vgroup"] != "unknown":
                return cls
        return best

    def _load_transcripts(self) -> Dict[str, List]:
        if self._transcripts is not None:
            return self._transcripts
        self._transcripts = {}
        if not self.transcript_interval_path:
            return self._transcripts
        try:
            import os

            import pandas as pd

            if not os.path.exists(self.transcript_interval_path):
                return self._transcripts
            df = pd.read_excel(self.transcript_interval_path)
            for _, row in df.iterrows():
                name = str(row.get("transcript") or "").strip()
                if not name:
                    continue
                intervals = []
                for col in df.columns:
                    if not str(col).startswith("interval"):
                        continue
                    val = row[col]
                    if pd.isna(val):
                        continue
                    parts = str(val).strip().split("-")
                    if len(parts) != 2:
                        continue
                    try:
                        a, b = int(parts[0]), int(parts[1])
                        intervals.append((min(a, b), max(a, b)))
                    except ValueError:
                        continue
                if intervals:
                    self._transcripts[name] = intervals
        except Exception as e:
            logger.warning(f"Failed to load transcript intervals: {e}")
        return self._transcripts

    def affected_isoforms(self, pos: int) -> List[str]:
        """回傳位置落在其外顯子區間內的 transcript 名稱清單。"""
        out = []
        for name, intervals in self._load_transcripts().items():
            if any(lo <= pos <= hi for lo, hi in intervals):
                out.append(name)
        return out

    def isoform_tissue_hint(self, band: str, isoforms: List[str]) -> Dict[str, str]:
        """
        根據論文 §1.5.1 / §1.5.3 的 isoform-組織對應，推斷變異落點對組織的暗示。

        - constitutive 區帶（Z-disk / A-band / M-band）：所有 isoform 都保留該外顯子，
          isoform 覆蓋無區分力 → 組織由 band 關聯與類比案例決定。
        - I-band：用心肌 isoform（N2BA/N2B）與骨骼肌 isoform（N2A）的覆蓋差異區分：
            同時覆蓋 → 心肌+骨骼肌共用 → 可能為 mixed (Both)
            只覆蓋心肌 → cardiac-predominant
            只覆蓋骨骼肌 → skeletal-predominant
        """
        cardiac = [i for i in isoforms if TRANSCRIPT_TISSUE_CLASS.get(i) == "cardiac"]
        skeletal = [i for i in isoforms if TRANSCRIPT_TISSUE_CLASS.get(i) == "skeletal"]

        if band in _CONSTITUTIVE_BANDS:
            return {
                "class": "constitutive",
                "text": (
                    f"{band} exons are constitutively expressed across all isoforms "
                    f"(cardiac {cardiac or '—'}, skeletal {skeletal or '—'}), so isoform "
                    f"coverage is NOT discriminative here; tissue is driven by band "
                    f"association and analogous variants."
                ),
            }

        if cardiac and skeletal:
            return {
                "class": "shared",
                "text": (
                    f"Exon is retained in BOTH cardiac ({', '.join(cardiac)}) and skeletal "
                    f"({', '.join(skeletal)}) isoforms → a mixed (Both) phenotype is possible."
                ),
            }
        if cardiac:
            return {
                "class": "cardiac",
                "text": (
                    f"Exon is retained mainly in cardiac isoforms ({', '.join(cardiac)}) and "
                    f"not the skeletal N2A isoform → cardiac-predominant effect expected."
                ),
            }
        if skeletal:
            return {
                "class": "skeletal",
                "text": (
                    f"Exon is retained in the skeletal isoform ({', '.join(skeletal)}) → "
                    f"skeletal-predominant effect expected."
                ),
            }
        return {
            "class": "uncertain",
            "text": "Isoform coverage could not be resolved for this position.",
        }

    # ── Step 3: 類比檢索 ─────────────────────────────────────────────────────

    def retrieve_analogs(
        self,
        pos: int,
        band: str,
        vgroup: str,
        top_k: Optional[int] = None,
        exon: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        從 HGMD 檢索已知致病變異當類比案例，依 genomic 距離由近到遠排序，取 top_k。

        retrieval_mode="band"（預設）：同 band + 同變異類型；同類型無結果才放寬成 band-only。

        retrieval_mode="exon"：改成分層檢索，優先用解析度更高的定位條件，湊不滿 top_k
        才往下一層放寬（RELATED_WORK.md：TTN 的表型專一性是外顯子級的，band 跨度 28 萬
        bp 太粗）。每筆類比會標上 match_tier，讓下游可以依定位強度加權：
            exon+type > exon > band+type > band
        """
        k = top_k or self.top_k
        kb = self._load_kb()

        if self.retrieval_mode != "exon":
            def match(row):
                if row["band"] != band:
                    return False
                if vgroup != "unknown" and row["vgroup"] != vgroup:
                    return False
                return True

            candidates = [r for r in kb if match(r)]
            used_band_only = False
            if not candidates and vgroup != "unknown":
                candidates = [r for r in kb if r["band"] == band]
                used_band_only = True

            candidates.sort(key=lambda r: abs(r["pos"] - pos))
            results = []
            for r in candidates[:k]:
                item = dict(r)
                item["distance"] = abs(r["pos"] - pos)
                item["band_only_match"] = used_band_only
                item["match_tier"] = "band" if used_band_only else "band+type"
                results.append(item)
            return results

        typed = vgroup != "unknown"
        # 層級順序是量出來的，不是猜的：把「同 exon 但不同變異類型」排在「同 band 同
        # 變異類型」之前，HGMD LOO 命中率會掉 2.1%（該層命中率只有 47.6%，對照組
        # 67.9%）。截斷 vs missense 是 TTN 的機制分水嶺（單倍劑量不足 vs 其他），
        # 比共處同一外顯子更能決定表型，所以類型比對優先於外顯子共位。
        tiers = [
            ("exon+type", lambda r: exon is not None and r["exon"] == exon and (not typed or r["vgroup"] == vgroup)),
            ("band+type", lambda r: r["band"] == band and (not typed or r["vgroup"] == vgroup)),
            ("exon", lambda r: exon is not None and r["exon"] == exon),
            ("band", lambda r: r["band"] == band),
        ]

        results: List[Dict[str, Any]] = []
        seen: Set[tuple] = set()
        for tier_name, pred in tiers:
            if len(results) >= k:
                break
            hits = [r for r in kb if (r["pos"], r["ref"], r["alt"]) not in seen and pred(r)]
            hits.sort(key=lambda r: abs(r["pos"] - pos))
            for r in hits[: k - len(results)]:
                item = dict(r)
                item["distance"] = abs(r["pos"] - pos)
                item["match_tier"] = tier_name
                # 沿用舊欄位語意：只靠 band 湊出來、沒對上變異類型的類比
                item["band_only_match"] = tier_name == "band"
                results.append(item)
                seen.add((r["pos"], r["ref"], r["alt"]))
        return results

    def hgmd_lookup(
        self, pos: int, ref: Optional[str] = None, alt: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        查 HGMD 是否收錄「同一位置」的已知致病變異。

        - exact：pos + ref + alt 完全相同 → 等同「這顆變異本身就是已知致病變異」，
          雖然 ClinVar 沒 cite 文獻，但 HGMD 已有註記，是最強證據（近似 ground truth）。
        - same_position：同 pos 但 ref/alt 不同（同一 codon 的其他置換）→ 極強證據。
        """
        kb = self._load_kb()
        try:
            p = int(pos)
        except (TypeError, ValueError):
            return {"exact": None, "same_position": []}
        same_pos = [dict(r) for r in kb if r["pos"] == p]
        exact = None
        # HGMD 完整匯出 CSV 沒有 REF/ALT 欄，庫內與查詢的 ref/alt 都是空字串，
        # 直接比字串會讓同位置的任何變異都被當成完全命中（實測 PubTator 有 49.2%
        # 中招，等於把答案讀回來）。這條路一律只回報同位置，exact 交給 main.py
        # Step 5.5 的 variant_db_lookup 處理，那裡靠 acc_num 接回 REF/ALT。
        if self.hgmd_csv_path:
            return {"exact": None, "same_position": same_pos}
        if ref is not None and alt is not None:
            rr, aa = str(ref).upper(), str(alt).upper()
            for r in same_pos:
                if str(r.get("ref", "")).upper() == rr and str(r.get("alt", "")).upper() == aa:
                    exact = r
                    break
        return {"exact": exact, "same_position": same_pos}

    def band_phenotype_prior(self, band: str) -> Dict[str, Any]:
        """該 band 內 HGMD 變異的 cardiac/skeletal/both 分佈（統計先驗）。"""
        kb = self._load_kb()
        counts = {"cardiac": 0, "skeletal": 0, "both": 0, "other": 0}
        for r in kb:
            if r["band"] != band:
                continue
            ph = r["phenotype"]
            if ph in counts:
                counts[ph] += 1
            else:
                counts["other"] += 1
        total = sum(counts.values())
        prop = {
            k: (v / total if total else 0.0) for k, v in counts.items()
        }
        return {"counts": counts, "proportions": prop, "total": total}

    # ── Step 2-3 彙整成 features ─────────────────────────────────────────────

    def build_features(
        self,
        variant_info: Union[str, Dict],
        clinvar_info: Optional[Dict] = None,
        evo2_result: Optional[Dict] = None,
    ) -> Dict[str, Any]:
        variant_info = self.normalize_variant_input(variant_info)
        pos = int(variant_info["pos"])
        ref = variant_info.get("ref")
        alt = variant_info.get("alt")
        band = self.assign_band(pos)
        cls = self.classify_input_variant(variant_info, clinvar_info)

        # HGMD 直接命中（同位置的已知致病變異）——最強證據，也可用來補變異類型
        hgmd_hit = self.hgmd_lookup(pos, ref, alt)
        exact = hgmd_hit["exact"]
        same_pos = hgmd_hit["same_position"]
        # ClinVar 無 HGVS 判不出類型時，改用 HGMD 命中列的 prot/dna（離線，不呼叫 VEP）
        if cls["vgroup"] == "unknown":
            src = exact or (same_pos[0] if same_pos else None)
            if src:
                cls = _classify_variant(src.get("prot", ""), src.get("dna", ""))

        exon = self.assign_exon(pos)
        isoforms = self.affected_isoforms(pos)
        iso_hint = self.isoform_tissue_hint(band, isoforms)
        # 外顯子 isoform 收錄度（PSI 代理）與它帶來的組織 likelihood ratio
        incl = self.exon_inclusion(pos)
        tissue_lr = self.inclusion_tissue_lr(incl.get("exon") if incl else None)
        analogs = self.retrieve_analogs(pos, band, cls["vgroup"], exon=exon)
        # 把「輸入變異本身」從類比清單移除（它會另外以 hgmd_exact_match 呈現）
        if exact:
            analogs = [
                a for a in analogs
                if not (
                    a.get("pos") == pos
                    and str(a.get("ref", "")).upper() == str(ref or "").upper()
                    and str(a.get("alt", "")).upper() == str(alt or "").upper()
                )
            ]
        prior = self.band_phenotype_prior(band)

        evo2 = {}
        if evo2_result and evo2_result.get("success"):
            evo2 = {
                "prediction": evo2_result.get("prediction"),
                "delta_score": evo2_result.get("delta_score"),
                "confidence": evo2_result.get("confidence"),
                # 論文 §2.2.3 定義的三級信心桶（predictor 本身未輸出，這裡補算）
                "confidence_level": evo2_confidence_bucket(evo2_result.get("delta_score")),
                "threshold": evo2_result.get("threshold"),
            }

        return {
            "variant_id": variant_info.get("variant_id"),
            "pos": pos,
            "ref": variant_info.get("ref"),
            "alt": variant_info.get("alt"),
            # HGVS 輸入時的來源資訊（供報告 / 除錯顯示）
            "hgvs_c": variant_info.get("hgvs_c"),
            "hgvs_p": variant_info.get("hgvs_p"),
            "transcript": variant_info.get("transcript"),
            "hgvs_resolved_by": variant_info.get("hgvs_resolved_by"),
            "hgvs_confidence": variant_info.get("hgvs_confidence"),
            "hgvs_warnings": variant_info.get("hgvs_warnings") or [],
            "band": band,
            "exon": exon,
            "exon_hotspot": self.exon_hotspot(exon, cls["vgroup"]),
            "vclass": cls["vclass"],
            "vgroup": cls["vgroup"],
            "isoforms": isoforms,
            "isoform_tissue_hint": iso_hint,
            "exon_inclusion": incl,
            "inclusion_tissue_lr": tissue_lr,
            "evo2": evo2,
            "analogs": analogs,
            "band_prior": prior,
            "hgmd_exact_match": exact,
            "hgmd_same_position": same_pos,
        }

    # ── Step 4: LLM 推論 prompt ──────────────────────────────────────────────

    @staticmethod
    def build_prompt(features: Dict[str, Any], is_medgemma: bool = True) -> str:
        v = features
        evo2 = v.get("evo2") or {}
        evo2_line = (
            f"{evo2.get('prediction', 'N/A')} "
            f"(delta_score={evo2.get('delta_score')}, "
            f"confidence={evo2.get('confidence_level', 'N/A')})"
            if evo2
            else "Not available"
        )

        # 類比案例表
        analog_lines = []
        for i, a in enumerate(v.get("analogs", []), 1):
            analog_lines.append(
                f"{i}. {a['ref']}>{a['alt']} @chr2:{a['pos']} | {a['prot'] or a['dna']} "
                f"| type={a['vclass']} | disease={a['phen']} "
                f"| tissue={a['phenotype']} | dist={a['distance']}bp"
                + (f" | exon={a['exon']}" if a.get("exon") else "")
                + (f" | match={a['match_tier']}" if a.get("match_tier") else "")
            )
        analog_block = "\n".join(analog_lines) if analog_lines else "(no analogous variants found)"

        prior = v.get("band_prior", {})
        prior_props = prior.get("proportions", {})
        prior_line = (
            f"cardiac={prior_props.get('cardiac', 0):.0%}, "
            f"skeletal={prior_props.get('skeletal', 0):.0%}, "
            f"both={prior_props.get('both', 0):.0%} "
            f"(n={prior.get('total', 0)} known pathogenic variants in this band)"
        )

        isoforms = ", ".join(v.get("isoforms", [])) or "unknown"
        iso_hint = v.get("isoform_tissue_hint", {}) or {}
        iso_hint_text = iso_hint.get("text", "Not available")

        # 外顯子熱點：文獻等級的定位證據，比 band 先驗強得多，要讓 LLM 優先採信
        hs = v.get("exon_hotspot")
        hotspot_line = (
            f"\n*** EXON HOTSPOT — exon {hs['exon']} ({hs['domain']}) ***\n"
            f"{hs['note']}\n"
            f"Unless the HGMD lookup below contradicts it, your inferred_disease should be "
            f"'{hs['disease']}' and inferred_tissue '{hs['tissue']}'."
            if hs
            else ""
        )

        # HGVS 輸入時，把 c./p. 與轉錄本一起帶進 prompt（LLM 判讀 splice/截斷位置有用）
        hgvs_parts = [
            p for p in (
                f"{v.get('transcript')}:{v.get('hgvs_c')}" if v.get("hgvs_c") else None,
                v.get("hgvs_p"),
            ) if p
        ]
        hgvs_line = f"\nHGVS: {'  '.join(hgvs_parts)}" if hgvs_parts else ""

        # HGMD 直接命中（同位置）——最強證據，放在最上面要求 LLM 以此為主
        exact = v.get("hgmd_exact_match")
        same_pos = v.get("hgmd_same_position") or []
        if exact:
            hgmd_note = (
                f"DIRECT HGMD MATCH — this EXACT variant is already catalogued in HGMD as pathogenic for "
                f"'{exact.get('phen')}' (tissue={exact.get('phenotype')}; {exact.get('prot') or exact.get('dna')}). "
                f"This is near-ground-truth. Your inferred_disease and inferred_tissue MUST equal this HGMD "
                f"annotation. It only reached this analogy step because ClinVar had no cited article."
            )
        elif same_pos:
            diseases = "; ".join(sorted({str(r.get('phen', '')) for r in same_pos if r.get('phen')}))
            hgmd_note = (
                f"SAME-POSITION HGMD VARIANT(S) at chr2:{v.get('pos')} (a different substitution at the identical "
                f"codon): {diseases}. This is very strong evidence — weight it ABOVE any more-distant analog."
            )
        else:
            hgmd_note = (
                "No exact or same-position HGMD match. Rely on the NEAREST analogs below (a 0 bp / same-codon "
                "analog is near-decisive; weight closer variants far more than distant ones)."
            )

        system_message = f"""You are a clinical genetics expert performing MECHANISTIC, ANALOGY-BASED reasoning for a TTN (titin) variant that has NO direct literature evidence.

IMPORTANT: There is NO published article describing this specific variant. You must NOT invent citations. Your output is an INFERENCE (a hypothesis), grounded ONLY in:
  (a) the Evo2 sequence-level pathogenicity prediction,
  (b) the variant's structural band and functional consequence,
  (c) which titin ISOFORMS retain the affected exon,
  (d) analogous KNOWN pathogenic TTN variants in the SAME band and of the SAME type,
  (e) the statistical phenotype prior of that band (HGMD).

## TTN MECHANISTIC BACKGROUND (use this to reason)
- Truncating variants (TTNtv: nonsense, frameshift, canonical splice-site) are ~an order of magnitude
  more likely to be pathogenic than missense variants, ESPECIALLY in the A-band. A-band TTNtv are the
  single most frequent genetic cause of dilated cardiomyopathy (DCM). The mechanism is haploinsufficiency
  of full-length titin (truncated transcript undergoes nonsense-mediated decay) → incomplete penetrance
  and a "second-hit" model, so a truncating call is a risk allele, not a deterministic diagnosis.
- Missense variants are very common in the population and are overwhelmingly benign or of uncertain
  significance; only assign a disease to a missense variant when the analogs strongly support it.

## ISOFORM → DISEASE LOGIC (titin isoforms, thesis §1.5.1 / §1.5.3)
- Cardiac isoforms: N2BA (long, compliant) and N2B (short, stiff) → dominant in HEART. Variants affecting
  these point to cardiac disease (DCM, and less commonly RCM/ARVC/LVNC).
- Skeletal isoform: N2A (longest, retains PEVK + proximal Ig, always lacks the N2B exon) → SKELETAL muscle.
  The C-terminus (M-band / near-C-terminal I-band) causes tibial muscular dystrophy (TMD) and LGMD2J;
  other skeletal titinopathies include HMERF and centronuclear myopathy.
- Z-disk / A-band / M-band exons are CONSTITUTIVELY expressed (present in ALL isoforms), so isoform
  coverage does not discriminate tissue there — rely on band association + analogs. Alternative splicing
  is confined to the I-band spring region, so isoform coverage is only discriminative for I-band variants.
- An exon retained ONLY in cardiac isoforms → cardiac-predominant; an exon SHARED by cardiac and skeletal
  isoforms → a mixed (Both) phenotype is possible.
- Caveat: splice-site variants can get an UNEXPECTEDLY LOW Evo2 score (effect is via mis-splicing) —
  do not down-weight a splice variant just because Evo2 is low.

## STRONGEST EVIDENCE — HGMD SAME-POSITION LOOKUP
{hgmd_note}

## TARGET VARIANT (no ClinVar-cited literature)
Variant: {v.get('variant_id')}  (chr2:{v.get('pos')} {v.get('ref')}>{v.get('alt')}){hgvs_line}
Structural band: {v.get('band')}
Exon (NM_001267550.2): {v.get('exon') or 'intronic / not in a coding exon'}{hotspot_line}
Variant type: {v.get('vclass')} ({v.get('vgroup')})
Affected isoforms: {isoforms}
Isoform-based tissue implication: {iso_hint_text}
Evo2 pathogenicity: {evo2_line}

## ANALOGOUS KNOWN PATHOGENIC VARIANTS (HGMD, same band + same type, nearest first)
NOTE: proximity matters — the closer the analog (smaller dist), the more it should weigh; the
nearest 1-2 analogs dominate. Do NOT let several distant analogs outvote a much closer one.
{analog_block}

## BAND-LEVEL PHENOTYPE PRIOR ({v.get('band')})
{prior_line}

## YOUR TASK
You MUST commit to a single best-guess disease for this variant — this is the whole point of the analysis.
Do NOT answer "Uncertain" unless there are literally ZERO analogous variants above. When the analogs list
several diseases, pick the MOST LIKELY one by (1) frequency among the nearest analogs, (2) the band prior,
(3) the isoform-based tissue implication, and (4) the Evo2 call (respecting the splice caveat). Name a
specific disease (e.g. "Dilated Cardiomyopathy (DCM)", "Tibial Muscular Dystrophy (TMD)"), not a vague
category, and list up to two alternative diagnoses. Always state clearly that this is an inference without
direct literature support and that penetrance may be incomplete.

## OUTPUT FORMAT
Respond with ONLY a JSON object:
```json
{{
  "inferred_disease": "<your single MOST LIKELY disease, full name — REQUIRED, do not say 'Uncertain' if any analog exists>",
  "inferred_tissue": "<Cardiac|Skeletal|Both|Uncertain>",
  "confidence": "<High|Moderate|Low>",
  "differential_diagnoses": ["<up to 2 alternative diseases, full names>", "..."],
  "rationale": "<2-4 sentences. MUST mention: the band, the variant type, the isoform/tissue implication, that it is mechanistically similar to specific analog variant(s), your predicted disease, and that there is NO direct literature evidence. Example: 'This is a truncating variant in the A-band, whose exon is retained in the cardiac N2BA/N2B isoforms; it is mechanistically similar to known DCM-causing TTNtv such as p.Q18753*, so it is predicted to be associated with dilated cardiomyopathy. No direct literature describes this variant and penetrance may be incomplete.'>",
  "analog_variants_used": ["<e.g. p.Q33396* (Tibial muscular dystrophy)>", "..."]
}}
```"""

        user_message = "Produce the analogy-based inference JSON. You MUST name a specific most-likely disease. Output ONLY the JSON object."

        if is_medgemma:
            return (
                f"<start_of_turn>user\n{system_message}\n\n{user_message}"
                f"<end_of_turn>\n<start_of_turn>model\n```json"
            )
        return (
            f"<|begin_of_text|><|start_header_id|>system<|end_header_id|>\n"
            f"{system_message}<|eot_id|>"
            f"<|start_header_id|>user<|end_header_id|>\n"
            f"{user_message}<|eot_id|>"
            f"<|start_header_id|>assistant<|end_header_id|>\n```json"
        )

    # ── 主入口 ───────────────────────────────────────────────────────────────

    def infer(
        self,
        extractor,
        variant_info: Union[str, Dict],
        clinvar_info: Optional[Dict] = None,
        evo2_result: Optional[Dict] = None,
    ) -> Dict[str, Any]:
        """
        執行完整 Track B 推論。

        extractor: 已載入 LLM 的 LocalClinicalExtractor（重用其 vLLM）；
                   若為 None，回傳只含 features 的規則式 fallback。
        variant_info: variant_info dict（含 pos/ref/alt），或直接給 HGVS 字串
                   （c./p./g. 或 `2-178527098-C-T`）。
        回傳 dict：features + inference（若有 LLM）+ mode 標記。
        """
        features = self.build_features(variant_info, clinvar_info, evo2_result)

        result: Dict[str, Any] = {
            "mode": "inference",
            "evidence_strength": "inferred_no_literature",
            "features": features,
        }

        if extractor is None:
            result["inference"] = self._rule_based_fallback(features)
            logger.info("Mechanistic inference: no LLM available -> rule-based fallback")
            return result

        is_medgemma = getattr(extractor, "is_medgemma", True)
        prompt = self.build_prompt(features, is_medgemma=is_medgemma)
        try:
            raw = self._generate(extractor, prompt, is_medgemma)
        except Exception as e:
            logger.error(f"Mechanistic LLM inference failed: {e}", exc_info=True)
            result["inference"] = self._rule_based_fallback(features)
            result["inference"]["source"] += "+after_llm_error"
            return result

        result["inference"] = self._finalize_inference(extractor, features, prompt, raw)
        return result

    def infer_batch(
        self,
        extractor,
        variants: List[Union[str, Dict]],
        clinvar_infos: Optional[List[Optional[Dict]]] = None,
        evo2_results: Optional[List[Optional[Dict]]] = None,
    ) -> List[Dict[str, Any]]:
        """
        一次推論多顆變異：所有 prompt 合併成單一次 vLLM generate 呼叫。

        逐顆呼叫時 GPU 幾乎在跑單一序列（約 4 tok/s），批次化可以把吞吐拉高
        一個量級，適合大批驗證。回傳順序與輸入一致。
        """
        n = len(variants)
        cis = clinvar_infos or [None] * n
        ers = evo2_results or [None] * n
        features = [self.build_features(v, cis[i], ers[i]) for i, v in enumerate(variants)]

        if extractor is None:
            out = []
            for f in features:
                inf = self._rule_based_fallback(f)
                out.append({
                    "mode": "inference",
                    "evidence_strength": "inferred_no_literature",
                    "features": f,
                    "inference": inf,
                })
            return out

        is_medgemma = getattr(extractor, "is_medgemma", True)
        prompts = [self.build_prompt(f, is_medgemma=is_medgemma) for f in features]
        try:
            raws = self._generate_many(extractor, prompts, is_medgemma)
        except Exception as e:
            logger.error(f"Mechanistic batch inference failed: {e}", exc_info=True)
            raws = [None] * n

        results = []
        for f, prompt, raw in zip(features, prompts, raws):
            if raw is None:
                inf = self._rule_based_fallback(f)
                inf["source"] += "+after_llm_error"
            else:
                inf = self._finalize_inference(extractor, f, prompt, raw)
            results.append({
                "mode": "inference",
                "evidence_strength": "inferred_no_literature",
                "features": f,
                "inference": inf,
            })
        return results

    def _finalize_inference(
        self, extractor, features: Dict[str, Any], prompt: str, raw: str
    ) -> Dict[str, Any]:
        """把 LLM 原始輸出解析、正規化，並套用 HGMD 命中覆寫與 Uncertain 保底。"""
        try:
            # 記錄 mechanistic prompt/response 以利除錯（重用 extractor 的 log 工具）
            try:
                extractor._log_prompt_response(
                    "mechanistic_inference", prompt, raw, features.get("variant_id")
                )
            except Exception:
                pass
            data = extractor._robust_json_parse(raw)
            data["inferred_tissue"] = self._normalize_tissue(data.get("inferred_tissue"))
            data["source"] = "llm_analogy"
            data["raw_llm_output"] = raw

            # HGMD 直接命中同一變異 → 以其註記為準，覆寫 LLM 的猜測（近似 ground truth）
            data = self._apply_hgmd_exact(data, features)

            # 保底：只要有類比案例，就不允許最終輸出停在 Uncertain。
            # LLM 面對異質類比清單常會 hedge 成 Uncertain（見 I-band 混合案例），
            # 這裡用類比多數決把疾病/組織補上，確保一定給出可解釋的猜測。
            data = self._backfill_if_uncertain(data, features)

            logger.info(
                f"Mechanistic inference -> disease={data.get('inferred_disease')} "
                f"tissue={data.get('inferred_tissue')} conf={data.get('confidence')} "
                f"source={data.get('source')}"
            )
            return data
        except Exception as e:
            logger.error(f"Mechanistic inference post-processing failed: {e}", exc_info=True)
            data = self._rule_based_fallback(features)
            data["source"] += "+after_llm_error"
            return data

    # ── LLM 生成（推論專用的 sampling / schema）──────────────────────────────

    # 推論輸出的 JSON schema。LocalClinicalExtractor 的 sampling_params 綁的是
    # 「逐篇文獻擷取」的 schema（disease / tissue_affected / …），直接沿用會讓
    # guided decoding 逼模型吐錯欄位、inferred_disease 永遠拿不到值，因此這裡
    # 另外建一組 sampling params。
    INFERENCE_JSON_SCHEMA: Dict[str, Any] = {
        "type": "object",
        "properties": {
            "inferred_disease": {"type": "string"},
            "inferred_tissue": {
                "type": "string",
                "enum": ["Cardiac", "Skeletal", "Both", "Uncertain"],
            },
            "confidence": {"type": "string", "enum": ["High", "Moderate", "Low"]},
            "differential_diagnoses": {"type": "array", "items": {"type": "string"}},
            "rationale": {"type": "string"},
            "analog_variants_used": {"type": "array", "items": {"type": "string"}},
        },
        "required": ["inferred_disease", "inferred_tissue", "confidence", "rationale"],
    }

    def _generate(self, extractor, prompt: str, is_medgemma: bool = True) -> str:
        return self._generate_many(extractor, [prompt], is_medgemma)[0]

    def _generate_many(
        self, extractor, prompts: List[str], is_medgemma: bool = True
    ) -> List[str]:
        """一次送出多個 prompt；取不到 vLLM 物件時退回 extractor 的逐筆方法。"""
        llm = getattr(extractor, "llm", None)
        params = None
        if llm is not None:
            params = getattr(extractor, "_mechanistic_sampling_params", None)
            if params is None:
                params = self._build_sampling_params(is_medgemma)
                try:
                    extractor._mechanistic_sampling_params = params
                except Exception:  # noqa: BLE001
                    pass
        if llm is None or params is None:
            return [extractor._generate_vllm(p) for p in prompts]
        outs = llm.generate(prompts, params)
        return [o.outputs[0].text for o in outs]

    def _build_sampling_params(self, is_medgemma: bool):
        try:
            from vllm import SamplingParams
        except ImportError:
            return None
        kwargs: Dict[str, Any] = {
            "n": 1,
            "temperature": 0.2,
            "top_p": 0.95,
            "max_tokens": 1536,
            "stop": ["<end_of_turn>"] if is_medgemma else ["<|eot_id|>"],
        }
        try:
            from vllm.sampling_params import StructuredOutputsParams

            kwargs["structured_outputs"] = StructuredOutputsParams(
                json=self.INFERENCE_JSON_SCHEMA
            )
        except Exception as e:  # noqa: BLE001
            logger.warning(f"推論的 guided JSON 無法啟用（改用純 sampling）：{e}")
        return SamplingParams(**kwargs)

    # ── HGMD 直接命中：以 DB 註記覆寫 LLM 猜測 ───────────────────────────────

    def _apply_hgmd_exact(
        self, data: Dict[str, Any], features: Dict[str, Any]
    ) -> Dict[str, Any]:
        """若 HGMD 收錄了同一顆變異（pos+ref+alt 相同），其疾病/組織近似 ground truth，
        直接覆寫 LLM 的預測，避免被較遠的異質類比案例帶偏。"""
        exact = features.get("hgmd_exact_match")
        if not exact:
            return data
        disease = str(exact.get("phen") or "").strip()
        if not disease:
            return data
        tissue = self._normalize_tissue(exact.get("phenotype"))
        prev = str(data.get("inferred_disease") or "").strip()

        data["inferred_disease"] = disease
        data["inferred_tissue"] = tissue
        data["confidence"] = "High"
        data["source"] = f"{data.get('source', 'llm_analogy')}+hgmd_exact_match"
        note = (
            f" [Note: this exact variant is directly catalogued in HGMD as {disease} ({tissue}); "
            f"the prediction was anchored to that HGMD annotation"
            + (
                f", overriding the model's initial guess ('{prev}')."
                if prev and prev.lower() != disease.lower()
                else "."
            )
            + "]"
        )
        data["rationale"] = str(data.get("rationale") or "").rstrip() + note

        # differential：LLM 若沒給，用最近的相異類比疾病補上
        if not data.get("differential_diagnoses"):
            diffs: List[str] = []
            for a in sorted(features.get("analogs", []), key=lambda x: x.get("distance", 1e18)):
                d = a.get("phen")
                if d and d != disease and d not in diffs:
                    diffs.append(d)
                if len(diffs) >= 2:
                    break
            data["differential_diagnoses"] = diffs
        return data

    # ── 保底：LLM 回 Uncertain 時用類比多數決補上 ─────────────────────────────

    def _backfill_if_uncertain(
        self, data: Dict[str, Any], features: Dict[str, Any]
    ) -> Dict[str, Any]:
        """若 LLM 給的 disease 是 Uncertain/空，但確實有類比案例，改用規則式多數決補上。"""
        analogs = features.get("analogs", [])
        disease = str(data.get("inferred_disease") or "").strip()
        is_uncertain = disease.lower() in ("", "uncertain", "not specified", "n/a", "none", "unknown")

        if not is_uncertain or not analogs:
            return data

        fb = self._rule_based_fallback(features)
        if str(fb.get("inferred_disease", "")).lower() in ("", "uncertain"):
            return data  # 類比案例也沒有可用疾病，維持原狀

        data["inferred_disease"] = fb["inferred_disease"]
        # tissue 若也 Uncertain 一併補
        if self._normalize_tissue(data.get("inferred_tissue")) == "Uncertain":
            data["inferred_tissue"] = fb["inferred_tissue"]
        # differential 若空也補
        if not data.get("differential_diagnoses"):
            data["differential_diagnoses"] = fb.get("differential_diagnoses", [])
        if not data.get("analog_variants_used"):
            data["analog_variants_used"] = fb.get("analog_variants_used", [])
        # rationale 附註來源
        note = (
            " [Note: the LLM was undecided among the heterogeneous analogs; the predicted "
            "disease was set to the most frequent disease among the nearest known variants.]"
        )
        data["rationale"] = str(data.get("rationale") or fb.get("rationale", "")).rstrip() + note
        data["source"] = f"{data.get('source', 'llm_analogy')}+analog_majority"
        if str(data.get("confidence") or "").strip().lower() not in ("low", "moderate", "high"):
            data["confidence"] = "Low"
        return data

    # ── Fallback（無 LLM 時）────────────────────────────────────────────────

    @staticmethod
    def _normalize_tissue(t: Optional[str]) -> str:
        if not t:
            return "Uncertain"
        s = str(t).strip().lower()
        if "both" in s:
            return "Both"
        if "cardiac" in s or "heart" in s:
            return "Cardiac"
        if "skeletal" in s or "muscle" in s:
            return "Skeletal"
        return "Uncertain"

    def _distance_weight(self, dist: Any, scale: Optional[float] = None) -> float:
        """把 genomic 距離轉成投票權重：越近權重越高。

        mode（由 distance_weight_mode 決定）：
          exponential  每 distance_scale bp 減半（預設；同位置權重=1）
          inverse      1 / (1 + d / distance_scale)，尾巴比指數衰減厚
          uniform      所有類比等權，等於純多數決、完全不看距離
          nearest      只有最近的那一組距離有投票權（勝者全拿）
        """
        try:
            d = max(float(dist), 0.0)
        except (TypeError, ValueError):
            return 1.0

        mode = self.distance_weight_mode
        if mode == "uniform":
            return 1.0
        s = float(scale if scale is not None else self.distance_scale)
        if s <= 0:
            return 1.0 if d == 0 else 0.0
        if mode == "inverse":
            return 1.0 / (1.0 + d / s)
        return 0.5 ** (d / s)

    # 定位層級加權：同一外顯子的類比落在同一個功能域，比單純同 band 可信得多。
    # 純距離看不出這件事——相隔 200bp 的兩顆變異可能跨了內含子、分屬不同結構域。
    # band 模式下所有類比都是同一層級，等於乘上同一常數，不影響排序。
    TIER_WEIGHTS = {"exon+type": 1.0, "band+type": 0.6, "exon": 0.4, "band": 0.2}

    def _tissue_correction_factors(self) -> Dict[str, float]:
        """目標基準率 / 類比庫基準率，用來把投票從 HGMD 的抽樣分佈校正回目標族群。"""
        if self._tissue_correction is not None:
            return self._tissue_correction
        pool: Dict[str, float] = {}
        for r in self._load_kb():
            ph = r.get("phenotype")
            if ph:
                pool[ph] = pool.get(ph, 0.0) + 1
        total = sum(pool.values()) or 1.0
        self._tissue_correction = {
            t: (self.tissue_prior.get(t, 0.0) / (n / total)) if n else 1.0
            for t, n in pool.items()
        }
        return self._tissue_correction

    def _analog_weight(
        self, analog: Dict[str, Any], tissue_lr: Optional[Dict[str, float]] = None
    ) -> float:
        """類比案例的投票權重 = 距離權重 × 定位層級權重 × 組織基準率校正 × 收錄度 LR。"""
        tier = analog.get("match_tier") or "band+type"
        w = self._distance_weight(analog.get("distance", 0)) * self.TIER_WEIGHTS.get(tier, 0.4)
        if self.tissue_prior:
            w *= self._tissue_correction_factors().get(analog.get("phenotype"), 1.0)
        if tissue_lr:
            # 查詢變異所在外顯子的收錄度，對「該類比的組織」打分：
            # 低收錄度外顯子會壓低心肌類比、抬高骨骼肌類比。
            w *= tissue_lr.get(analog.get("phenotype"), 1.0)
        if self.rankscore_weighting > 0:
            rs = analog.get("rankscore")
            # 缺 RANKSCORE 的列（約 5%）給中位數 0.75，不獎不罰
            rs = float(rs) if rs is not None else 0.75
            s = self.rankscore_weighting
            w *= (1.0 - s) + s * rs
        if self.vague_penalty != 1.0 and analog.get("phen_vague"):
            w *= self.vague_penalty
        return w

    # ── 階層式收縮預測（predict_mode="hierarchical"）──────────────────────

    def _hier_tallies(self):
        """統計每個外顯子 / PSI 桶 / band / 全庫的疾病與組織計數。

        以 self._kb 的物件識別當快取鍵：評估時 run() 會整包換掉 _kb 做
        leave-one-out，換掉後必須重算，否則會拿到含答案的統計。
        """
        kb = self._load_kb()
        if self._hier_cache is not None and self._hier_cache[0] is kb:
            return self._hier_cache[1]

        from collections import Counter, defaultdict

        incl = self._load_exon_inclusion()
        lv = {k: (defaultdict(Counter), defaultdict(Counter))
              for k in ("exon", "bucket", "band")}
        g_d, g_t = Counter(), Counter()
        for r in kb:
            d, t = r.get("phen"), r.get("phenotype")
            if not d or t not in ("cardiac", "skeletal", "both"):
                continue
            ex = r.get("exon")
            keys = {"exon": ex,
                    "bucket": (incl.get(ex) or {}).get("bucket") if ex else None,
                    "band": r.get("band")}
            for name, k in keys.items():
                if k is not None:
                    lv[name][0][k][d] += 1
                    lv[name][1][k][t] += 1
            g_d[d] += 1
            g_t[t] += 1
        out = (lv, g_d, g_t)
        self._hier_cache = (kb, out)
        return out

    def hierarchical_scores(self, exon, band):
        """回傳 (疾病機率 dict, 組織機率 dict)，由全庫逐層收縮到該外顯子。"""
        from collections import Counter

        lv, g_d, g_t = self._hier_tallies()
        incl = self._load_exon_inclusion()
        bucket = (incl.get(exon) or {}).get("bucket") if exon else None

        def chain(idx, glob):
            n_g = sum(glob.values()) or 1
            steps = [(lv["band"][idx].get(band, Counter()), self.m_band),
                     (lv["bucket"][idx].get(bucket, Counter()), self.m_bucket),
                     (lv["exon"][idx].get(exon, Counter()), self.m_exon)]
            out = {}
            for lab in glob:
                p = glob[lab] / n_g
                for cnt, m in steps:          # 由粗到細，每層往上一層收縮
                    p = (cnt.get(lab, 0) + m * p) / (sum(cnt.values()) + m)
                out[lab] = p
            return out

        return chain(0, g_d), chain(1, g_t)

    def _hierarchical_predict(self, features: Dict[str, Any]) -> Dict[str, Any]:
        exon, band = features.get("exon"), features.get("band")
        d_scores, t_scores = self.hierarchical_scores(exon, band)
        if not d_scores:
            return {"inferred_disease": "Uncertain", "inferred_tissue": "Uncertain",
                    "confidence": "Low", "differential_diagnoses": [],
                    "rationale": "No annotated variants available for this region.",
                    "analog_variants_used": [], "source": "hierarchical_no_data"}

        ranked = sorted(d_scores.items(), key=lambda kv: -kv[1])
        top, p_top = ranked[0]
        tissue = max(t_scores, key=t_scores.get) if t_scores else "uncertain"
        incl = (self._load_exon_inclusion().get(exon) or {}) if exon else {}
        psi_txt = (f"cardiac PSI {incl['psi_heart']:.0f}%, skeletal PSI "
                   f"{incl['psi_skeletal']:.0f}%"
                   if incl.get("psi_heart") is not None else "PSI unavailable")
        return {
            "inferred_disease": top,
            "inferred_tissue": self._normalize_tissue(tissue),
            "confidence": "High" if p_top >= 0.6 else ("Moderate" if p_top >= 0.4 else "Low"),
            "differential_diagnoses": [d for d, _ in ranked[1:3]],
            "rationale": (
                f"Exon {exon} ({band}; {psi_txt}) is estimated to be associated with "
                f"{top} (p = {p_top:.2f}) by hierarchical shrinkage over the exon, its "
                f"exon-usage class, and the structural band. No direct literature "
                f"describes this variant and penetrance may be incomplete."
            ),
            "analog_variants_used": [],
            "source": "hierarchical",
        }

    def _rule_based_fallback(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """無 LLM 時（或 LLM 失敗時）的保守推論：
        1) HGMD 同一變異直接命中 → 採用其疾病/組織（近似 ground truth）。
        2) 否則用「距離加權」的類比多數決，最近鄰主導。"""
        analogs = features.get("analogs", [])
        prior = features.get("band_prior", {}).get("proportions", {})
        band = features.get("band")
        vgroup = features.get("vgroup")
        exact = features.get("hgmd_exact_match")

        # ── 1) HGMD 直接命中同一變異 ──────────────────────────────────────────
        if exact:
            disease = str(exact.get("phen") or "").strip() or "Uncertain"
            tissue = self._normalize_tissue(exact.get("phenotype"))
            differential = []
            for a in sorted(analogs, key=lambda x: x.get("distance", 1e18)):
                d = a.get("phen")
                if d and d != disease and d not in differential:
                    differential.append(d)
                if len(differential) >= 2:
                    break
            rationale = (
                f"This EXACT variant (chr2:{features.get('pos')} {features.get('ref')}>{features.get('alt')}, "
                f"{exact.get('prot') or exact.get('dna')}) is catalogued in HGMD as pathogenic for "
                f"{disease} ({tissue}). It lies in the {band} and is a {features.get('vclass')} variant. "
                f"ClinVar has no cited article for it, so it is surfaced via the HGMD knowledge base rather "
                f"than as a de novo inference; penetrance may still be incomplete."
            )
            return {
                "inferred_disease": disease,
                "inferred_tissue": tissue,
                "confidence": "High",
                "differential_diagnoses": differential,
                "rationale": rationale,
                "analog_variants_used": [f"{exact.get('prot') or exact.get('dna')} ({disease}, 0bp exact)"],
                "source": "hgmd_exact_match",
            }

        # ── 2) 文獻確立的外顯子熱點（ACMG PM1 等級的定位證據）──────────────────
        hotspot = features.get("exon_hotspot")
        if hotspot:
            differential = []
            for a in sorted(analogs, key=lambda x: x.get("distance", 1e18)):
                d = a.get("phen")
                if d and d != hotspot["disease"] and d not in differential:
                    differential.append(d)
                if len(differential) >= 2:
                    break
            nearest = sorted(analogs, key=lambda x: x.get("distance", 1e18))[:3]
            return {
                "inferred_disease": hotspot["disease"],
                "inferred_tissue": self._normalize_tissue(hotspot["tissue"]),
                "confidence": "High" if hotspot.get("strength") == "strong" else "Moderate",
                "differential_diagnoses": differential,
                "rationale": (
                    f"This {features.get('vclass')} variant falls in exon {hotspot['exon']} "
                    f"({hotspot['domain']}), an established TTN mutation hotspot. "
                    f"{hotspot['note']} The variant is therefore predicted to be associated with "
                    f"{hotspot['disease']}. No direct literature describes this specific variant; "
                    f"this is a mechanistic inference and penetrance may be incomplete."
                ),
                "analog_variants_used": [
                    f"{a['prot'] or a['dna']} ({a['phen']}, {a.get('distance', '?')}bp)"
                    for a in nearest
                ],
                "source": f"exon_hotspot_{hotspot['exon']}",
            }

        # ── 3) 依 predict_mode 決定：階層式收縮 或 距離加權投票 ────────────────
        if self.predict_mode == "hierarchical":
            return self._hierarchical_predict(features)

        if self.distance_weight_mode == "nearest" and analogs:
            # 勝者全拿：只讓距離最近的那一組（可能並列）投票
            nearest_d = min(a.get("distance", 0) or 0 for a in analogs)
            analogs = [a for a in analogs if (a.get("distance", 0) or 0) == nearest_d]

        # 查詢變異所在外顯子的收錄度 LR：低收錄度（替代性剪接多、心肌 PSI 低）
        # 會把票數推向骨骼肌。詳見 inclusion_tissue_lr()。
        tissue_lr = features.get("inclusion_tissue_lr") or None

        tissue_w: Dict[str, float] = {}
        disease_w: Dict[str, float] = {}
        for a in analogs:
            w = self._analog_weight(a, tissue_lr)
            ph = a.get("phenotype")
            if ph:
                tissue_w[ph] = tissue_w.get(ph, 0.0) + w
            if a.get("phen"):
                disease_w[a["phen"]] = disease_w.get(a["phen"], 0.0) + w

        ranked_diseases = sorted(disease_w.items(), key=lambda kv: -kv[1])
        top_disease = ranked_diseases[0][0] if ranked_diseases else "Uncertain"
        differential = [d for d, _ in ranked_diseases[1:3]]

        iso_hint = features.get("isoform_tissue_hint", {}) or {}
        iso_class = iso_hint.get("class")

        # tissue：只有「單一組織」的 isoform 覆蓋才有決定性；shared/constitutive 太弱，
        # 不用來壓過就近的類比案例。否則讓 tissue 與 top_disease 一致（距離加權）。
        if iso_class in ("cardiac", "skeletal"):
            top_tissue = iso_class
        else:
            dt_w: Dict[str, float] = {}
            for a in analogs:
                if a.get("phen") == top_disease and a.get("phenotype"):
                    dt_w[a["phenotype"]] = dt_w.get(a["phenotype"], 0.0) + self._analog_weight(a, tissue_lr)
            if dt_w:
                top_tissue = max(dt_w, key=dt_w.get)
            elif tissue_w:
                top_tissue = max(tissue_w, key=tissue_w.get)
            elif prior:
                cs = {k: prior.get(k, 0) for k in ("cardiac", "skeletal")}
                top_tissue = max(cs, key=cs.get) if any(cs.values()) else "uncertain"
            else:
                top_tissue = "uncertain"

        nearest = sorted(analogs, key=lambda x: x.get("distance", 1e18))[:3]
        used = [f"{a['prot'] or a['dna']} ({a['phen']}, {a.get('distance', '?')}bp)" for a in nearest]
        rationale = (
            f"This is a {features.get('vclass')} variant in the {band}. "
            f"{iso_hint.get('text', '')} "
            f"By distance-weighted analogy to the NEAREST known pathogenic {vgroup} variants "
            f"({'; '.join(used) if used else 'none found'}), it is predicted to be "
            f"associated with {top_disease}. No direct literature evidence exists for "
            f"this variant, and penetrance may be incomplete; this is a mechanistic inference only."
        )
        return {
            "inferred_disease": top_disease,
            "inferred_tissue": self._normalize_tissue(top_tissue),
            "confidence": "Low",
            "differential_diagnoses": differential,
            "rationale": rationale,
            "analog_variants_used": used,
            "source": "analogy_voting",
        }
