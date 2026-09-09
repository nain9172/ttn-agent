#!/usr/bin/env python3
"""
TTN 疾病名稱正規化與比對。

HGMD 的 PHEN（"Cardiomyopathy; dilated"）、PubTator 抽取的 disease
（"Dilated Cardiomyopathy (DCM)"）、LLM 推論輸出（"Dilated Cardiomyopathy"）
寫法都不同，直接字串比對會嚴重低估正確率。這裡把它們映射到同一組
canonical label，並用「疾病階層」判斷兩者是相同、相容（一般↔特定）、
只有組織相同、或完全不同。

用法：
    from utils.disease_normalizer import compare_disease, canonical_labels
    compare_disease("Dilated Cardiomyopathy (DCM)", "Cardiomyopathy; dilated")
    # -> {'level': 'exact', ...}
"""

from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Optional, Set

# ── canonical label 的關鍵字（順序重要：特定疾病要在一般疾病之前）──────────────
# label -> regex（case-insensitive，會找出字串中「所有」命中的 label）
_PATTERNS: List[tuple] = [
    # 心肌 —— 特定
    ("DCM", r"dilated\s+cardiomyopath|cardiomyopathy[;,]\s*dilated|\bDCM\b|\bDCMP\b|\bIDC\b"),
    ("HCM", r"hypertroph\w*\s+cardiomyopath|cardiomyopathy[;,]\s*hypertrophic|\bHCM\b"),
    ("ARVC", r"arrhythmogenic\s+right\s+ventricular|arrhythmogenic\s+cardiomyopath|\bARVC\b|\bARVD\b|\bACM\b"),
    ("RCM", r"restrictive\s+cardiomyopath|\bRCM\b"),
    ("LVNC", r"non-?compaction|noncompaction|\bLVNC\b|\bLVNCC\b"),
    ("PPCM", r"peripartum"),
    ("MYOCARDITIS", r"myocarditis"),
    # 心臟 —— 其他
    ("ARRHYTHMIA", r"atrial\s+fibrillation|long\s*qt|\bLQT\w*\b|brugada|\bBrS\b|\bCPVT\b|"
                   r"ventricular\s+tachycardia|sick\s+sinus|\bSSS\b|arrhythmi|conduction\s+(defect|disease)|\bAF\b"),
    ("SUDDEN_DEATH", r"sudden\s+(unexplained|unexpected|cardiac|arrhythmic|infant)|\bSIDS\b|"
                     r"\bSUDS?\b|\bSCD\b|\bSADS\b|sudden\s+death"),
    ("CHD", r"congenital\s+heart|septal\s+defect|\bASD\b|\bVSD\b|ebstein|tetralogy|patent\s+ductus"),
    ("HEART_FAILURE", r"heart\s+failure|cardiac\s+failure"),
    # 骨骼肌 —— 特定
    ("EOMFC", r"early[-\s]?onset.{0,40}(fatal\s+)?cardiomyopath|\bEOMFC\b|salih"),
    ("TMD", r"tibial\s+muscular\s+dystroph|\bTMD\b|udd\s+myopath"),
    ("LGMD", r"limb[-\s]?girdle|\bLGMD\w*\b"),
    ("HMERF", r"early\s+respiratory\s+failure|\bHMERF\b"),
    ("CNM", r"centronuclear|\bCNM\b|internal\s+nuclei"),
    # 不收 \bMmD\b：pattern 是 case-insensitive，會把 Moyamoya Disease (MMD) 這種
    # 和 TTN 無關的疾病誤判成 multi-minicore（資料裡有 32 筆）。真正寫 MmD 的地方
    # 一定同時把 minicore 拼出來（8/8 筆都是），靠下面的 minicore 就抓得到。
    ("CORE_MYOPATHY", r"multi-?minicore|minicore|core\s+myopathy|core\s+disease"),
    ("CFTD", r"fibre\s+type\s+disproportion|fiber\s+type\s+disproportion|\bCFTD\b"),
    ("ARTHROGRYPOSIS", r"arthrogryposis|fetal\s+akinesia|hydrops\s+fetalis|\bAMC\b"),
    ("RIMMED_VACUOLAR", r"rimmed\s+vacuolar"),
    ("DISTAL_MYOPATHY", r"distal\s+titinopath|myopathy[;,]\s*distal|distal\s+myopath"),
    ("CONGENITAL_MYOPATHY", r"congenital\s+myopath|myopathy[;,]\s*congenital|congenital\s+titinopath|"
                            r"childhood[-\s]onset\s+myopath|congenital\s+muscular"),
    ("MUSCULAR_DYSTROPHY", r"muscular\s+dystroph|\bMD\b|dystroph"),
    # 一般
    # 不收單獨的 "CM"：資料裡它同時是 congenital myopathy 與 cardiomyopathy 的縮寫
    # （"Congenital Myopathy (CM)" vs "Cardiomyopathy (CM)"），而後者一定也把字
    # 拼出來，靠 cardiomyopath 就抓得到，收 \bCM\b 只會把骨骼肌病誤判成心肌病。
    ("CARDIOMYOPATHY", r"cardiomyopath|cardiac\s+hypertrophy"),
    ("CARDIOSKELETAL", r"cardioskeletal|neuromuscular/cardioskeletal"),
    # 「myopathy」不可命中 cardiomyopathy（否則每個心肌病都會被歸成骨骼肌病）
    ("MYOPATHY", r"(?<!cardio)(?<!cardio-)(?<!cardio )myopath|myofibrillar|"
                 r"muscle\s+weakness|muscle\s+disease|myositis|rhabdomyolysis"),
    ("TITINOPATHY", r"titinopath"),
    ("NEUROMUSCULAR", r"neuromuscular"),
    ("CARDIAC_OTHER", r"cardiovascular\s+disease|cardiac\s+disease|\bcardiac\b|\bheart\b"),
    # 與 TTN 機制無關（癌症、其他系統疾病）→ 評分時單獨列出
    ("CANCER", r"cancer|carcinom|adenocarcinom|tumou?r|neoplas|melanom|leuk[a]?emi|lymphom|sarcom|glioma|\bPCa\b|\bPDAC\b"),
]

_NOT_SPECIFIED_RE = re.compile(
    r"^\s*(not\s+(la\s+)?specified|not\s+available|others?|none|n/?a|unknown|uncertain|-{0,2})\s*$",
    re.IGNORECASE,
)

# label -> 上位（較一般）label 集合，用來判斷「相容」
_PARENTS: Dict[str, Set[str]] = {
    "DCM": {"CARDIOMYOPATHY"},
    "HCM": {"CARDIOMYOPATHY"},
    "ARVC": {"CARDIOMYOPATHY"},
    "RCM": {"CARDIOMYOPATHY"},
    "LVNC": {"CARDIOMYOPATHY"},
    "PPCM": {"CARDIOMYOPATHY"},
    "MYOCARDITIS": {"CARDIAC_OTHER"},
    "CARDIOMYOPATHY": {"CARDIAC_OTHER"},
    # TITINOPATHY 是最沒有資訊量的上位詞：任何 TTN 相關疾病都算它的下位
    "CARDIAC_OTHER": {"TITINOPATHY"},
    "NEUROMUSCULAR": {"TITINOPATHY"},
    "ARRHYTHMIA": {"CARDIAC_OTHER"},
    "SUDDEN_DEATH": {"CARDIAC_OTHER"},
    "CHD": {"CARDIAC_OTHER"},
    "HEART_FAILURE": {"CARDIAC_OTHER"},
    "TMD": {"DISTAL_MYOPATHY", "MUSCULAR_DYSTROPHY"},
    "LGMD": {"MUSCULAR_DYSTROPHY"},
    "HMERF": {"MYOPATHY"},
    "CNM": {"CONGENITAL_MYOPATHY"},
    "CORE_MYOPATHY": {"CONGENITAL_MYOPATHY"},
    "CFTD": {"CONGENITAL_MYOPATHY"},
    "ARTHROGRYPOSIS": {"CONGENITAL_MYOPATHY"},
    "RIMMED_VACUOLAR": {"MYOPATHY"},
    "DISTAL_MYOPATHY": {"MYOPATHY"},
    "CONGENITAL_MYOPATHY": {"MYOPATHY"},
    "MUSCULAR_DYSTROPHY": {"MYOPATHY"},
    # 早發型肌病合併致命心肌病：骨骼肌 + 心肌兩軌
    "EOMFC": {"CONGENITAL_MYOPATHY", "CARDIOMYOPATHY"},
    "CARDIOSKELETAL": {"MYOPATHY", "CARDIAC_OTHER"},
    "MYOPATHY": {"NEUROMUSCULAR"},
}

# 每個 label 對應的組織（Cardiac / Skeletal / Both / Other）
_TISSUE: Dict[str, str] = {
    "DCM": "Cardiac", "HCM": "Cardiac", "ARVC": "Cardiac", "RCM": "Cardiac",
    "LVNC": "Cardiac", "PPCM": "Cardiac", "MYOCARDITIS": "Cardiac",
    "CARDIOMYOPATHY": "Cardiac", "ARRHYTHMIA": "Cardiac", "SUDDEN_DEATH": "Cardiac",
    "CHD": "Cardiac", "HEART_FAILURE": "Cardiac", "CARDIAC_OTHER": "Cardiac",
    "TMD": "Skeletal", "LGMD": "Skeletal", "HMERF": "Skeletal", "CNM": "Skeletal",
    "CORE_MYOPATHY": "Skeletal", "CFTD": "Skeletal", "ARTHROGRYPOSIS": "Skeletal",
    "RIMMED_VACUOLAR": "Skeletal", "DISTAL_MYOPATHY": "Skeletal",
    "CONGENITAL_MYOPATHY": "Skeletal", "MUSCULAR_DYSTROPHY": "Skeletal",
    "MYOPATHY": "Skeletal", "NEUROMUSCULAR": "Skeletal",
    "EOMFC": "Both", "CARDIOSKELETAL": "Both",
    "TITINOPATHY": "Both",
    "CANCER": "Other",
}

# 「一般到幾乎沒有資訊」的 label。預測只給這種答案時雖然能算 compatible
# （"Titinopathy" 對任何 TTN 變異都成立），但幾乎沒有臨床價值，
# 因此在比對結果標記出來，讓報表能把它和真正具體的命中分開看。
_VAGUE = {"TITINOPATHY", "NEUROMUSCULAR", "CARDIAC_OTHER", "MYOPATHY", "CARDIOMYOPATHY"}

_TTN_RELEVANT = set(_TISSUE) - {"CANCER"}


def is_unspecified(text: Optional[str]) -> bool:
    return not text or bool(_NOT_SPECIFIED_RE.match(str(text)))


def is_vague_disease(text: Optional[str]) -> bool:
    """疾病名是否只落在「幾乎沒有資訊」的上位詞（Cardiomyopathy / Titinopathy …）。

    HGMD 有約三成的 PHEN 屬於這類，拿它們當類比案例投票會把預測拉向同樣籠統的
    答案，所以檢索端也需要判斷這件事，不是只有評分時要用。
    """
    labels = canonical_labels(text)
    return bool(labels) and all(x in _VAGUE for x in labels)


def canonical_labels(text: Optional[str]) -> List[str]:
    """把自由文字疾病名映射成 canonical label 清單（可能多個，例如合併表型）。"""
    if is_unspecified(text):
        return []
    s = str(text)
    out: List[str] = []
    for label, pattern in _PATTERNS:
        if re.search(pattern, s, re.IGNORECASE) and label not in out:
            out.append(label)
    # 只保留最具體的 label：凡是其他命中 label 的上位詞都丟掉。
    # "Tibial Muscular Dystrophy (TMD)" 會同時命中 TMD 與 MUSCULAR_DYSTROPHY，
    # 若把上位詞留著，它就會和 LGMD2J（LGMD + MUSCULAR_DYSTROPHY）共用一個 label
    # 而被判成 exact——兩個其實不同的表型會被算成命中。
    return [x for x in out if not _is_ancestor_of_any(x, out)]


def _is_ancestor_of_any(label: str, labels: Iterable[str]) -> bool:
    return any(other != label and label in ancestors(other) for other in labels)


def ancestors(label: str, _seen: Optional[Set[str]] = None) -> Set[str]:
    """label 的所有上位 label（含間接）。"""
    seen = _seen or set()
    for p in _PARENTS.get(label, set()):
        if p in seen:
            continue
        seen.add(p)
        ancestors(p, seen)
    return seen


def tissue_of_labels(labels: Iterable[str]) -> str:
    """一組 label 的組織歸屬（Cardiac / Skeletal / Both / Other / Uncertain）。"""
    tissues = {_TISSUE.get(x) for x in labels if _TISSUE.get(x)}
    tissues.discard(None)
    if "Both" in tissues or {"Cardiac", "Skeletal"} <= tissues:
        return "Both"
    for t in ("Cardiac", "Skeletal", "Other"):
        if t in tissues:
            return t
    return "Uncertain"


def normalize_tissue(text: Optional[str]) -> str:
    """把 tissue 字串（或疾病名）正規化成 Cardiac / Skeletal / Both / Uncertain。"""
    if is_unspecified(text):
        return "Uncertain"
    s = str(text).strip().lower()
    if "both" in s:
        return "Both"
    if s in ("cardiac", "heart", "cardiac muscle"):
        return "Cardiac"
    if s in ("skeletal", "skeletal muscle", "muscle"):
        return "Skeletal"
    if s in ("others", "other"):
        return "Other"
    return tissue_of_labels(canonical_labels(text))


def is_ttn_relevant(labels: Iterable[str]) -> bool:
    """是否屬於 TTN 機制可解釋的表型（排除癌症等旁支關聯）。"""
    ls = list(labels)
    return bool(ls) and any(x in _TTN_RELEVANT for x in ls)


def _pair_level(pred: str, truth: str) -> Optional[str]:
    if pred == truth:
        return "exact"
    if truth in ancestors(pred) or pred in ancestors(truth):
        return "compatible"
    return None


def compare_disease(
    predicted: Optional[str],
    truth: Optional[str],
    differentials: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """
    比對推論疾病與 ground truth。

    level：
      exact       canonical label 完全相同（DCM vs "Cardiomyopathy; dilated"）
      compatible  一般 ↔ 特定（"Cardiomyopathy" vs DCM、Titinopathy vs TMD）
      differential 主答案沒中，但列在 differential_diagnoses 裡
      same_tissue 疾病不同但組織相同（DCM vs HCM）
      miss        不同組織 / 完全不同
      n/a         ground truth 或預測沒有可比對的疾病
    """
    pred_labels = canonical_labels(predicted)
    truth_labels = canonical_labels(truth)
    result: Dict[str, Any] = {
        "level": "n/a",
        "pred_labels": pred_labels,
        "truth_labels": truth_labels,
        "pred_tissue": tissue_of_labels(pred_labels),
        "truth_tissue": tissue_of_labels(truth_labels),
        "pred_vague": bool(pred_labels) and all(x in _VAGUE for x in pred_labels),
    }
    if not truth_labels or not pred_labels:
        return result

    best = None
    for p in pred_labels:
        for t in truth_labels:
            lvl = _pair_level(p, t)
            if lvl == "exact":
                best = "exact"
                break
            if lvl == "compatible" and best != "exact":
                best = "compatible"
        if best == "exact":
            break

    if best is None and differentials:
        for d in differentials:
            for dl in canonical_labels(d):
                for t in truth_labels:
                    if _pair_level(dl, t):
                        best = "differential"
                        break
                if best:
                    break
            if best:
                break

    if best is None:
        pt, tt = result["pred_tissue"], result["truth_tissue"]
        if pt != "Uncertain" and (pt == tt or "Both" in (pt, tt)):
            best = "same_tissue"
        else:
            best = "miss"
    result["level"] = best
    return result


def tissue_match(pred: str, truth: str) -> str:
    """組織比對：Both 視為與任一單一組織相容。"""
    if pred in ("Uncertain", "Other") or truth in ("Uncertain", "Other"):
        return "n/a"
    if pred == truth:
        return "exact"
    if "Both" in (pred, truth):
        return "partial"
    return "miss"


__all__ = [
    "canonical_labels",
    "compare_disease",
    "is_ttn_relevant",
    "is_unspecified",
    "normalize_tissue",
    "tissue_match",
    "tissue_of_labels",
]
