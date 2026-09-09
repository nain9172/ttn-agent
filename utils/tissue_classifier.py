#!/usr/bin/env python3
"""
組織（心肌 / 骨骼肌）判定器。

為什麼要獨立出來：benchmark_inference.py 量到組織判定是整個 Track B 的閘門——
組織判對時疾病命中率 69.8%，判錯時只剩 10.0%。但現行邏輯把組織當成疾病投票的
副產品（先選出票數最高的疾病，再看支持它的類比是什麼組織），所以疾病投票一旦
被 HGMD 的收錄偏誤帶偏，組織跟著錯，疾病也就不可能對。

這裡改成獨立預測：組織有自己的一組特徵（isoform 覆蓋、外顯子、變異類型、局部
類比的組織分佈、band、基因上的相對位置），先判組織，再用組織去約束疾病選擇。

特徵抽取（tissue_features）刻意寫成 build_features 輸出的純函式，訓練與線上推論
共用同一份程式碼，避免 train/serve skew。

模型是 5~15 個特徵的 logistic regression / gradient boosting，不是 fine-tune LLM；
定位上等同 RELATED_WORK.md 裡 PM5 校準論文做的「把局部證據轉成校準過的機率」。
"""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# TTN 在 GRCh38 chr2 的區間，用來把座標正規化成 0~1
TTN_START, TTN_END = 178525989, 178807423
N_EXONS = 364  # 文獻編號的外顯子總數

# 類比組織投票的距離半衰期。和 MechanisticInference 的預設一致，但這裡刻意寫死：
# 特徵定義必須和訓練當下完全一樣，不能被呼叫端的參數改掉。
ANALOG_HALF_LIFE_BP = 150.0

CARDIAC_ISOFORMS = ("N2BA", "N2B", "Cronos")
SKELETAL_ISOFORMS = ("N2A",)

FEATURE_NAMES = [
    "is_truncating",
    "is_missense",
    "band_Z_disk",
    "band_I_band",
    "band_A_band",
    "band_M_band",
    "iso_cardiac",
    "iso_skeletal",
    "iso_cardiac_only",
    "iso_skeletal_only",
    "in_exon",
    "exon_norm",
    "pos_norm",
    "n_analogs",
    "log_nearest_bp",
    "analog_cardiac_frac",
    "analog_cardiac_frac_weighted",
    "nearest_is_cardiac",
]

_BAND_KEYS = {
    "Z-disk": "band_Z_disk",
    "I-band": "band_I_band",
    "A-band": "band_A_band",
    "M-band": "band_M_band",
}


def tissue_features(features: Dict[str, Any]) -> Dict[str, float]:
    """把 MechanisticInference.build_features() 的輸出轉成組織分類器的數值特徵。"""
    out = {name: 0.0 for name in FEATURE_NAMES}

    vgroup = features.get("vgroup")
    out["is_truncating"] = 1.0 if vgroup == "truncating" else 0.0
    out["is_missense"] = 1.0 if vgroup == "missense" else 0.0

    band_key = _BAND_KEYS.get(features.get("band"))
    if band_key:
        out[band_key] = 1.0

    isoforms = features.get("isoforms") or []
    has_cardiac = any(i in isoforms for i in CARDIAC_ISOFORMS)
    has_skeletal = any(i in isoforms for i in SKELETAL_ISOFORMS)
    out["iso_cardiac"] = float(has_cardiac)
    out["iso_skeletal"] = float(has_skeletal)
    # 只被單一組織的 isoform 保留，才是有鑑別力的訊號；兩邊都保留等於沒說
    out["iso_cardiac_only"] = float(has_cardiac and not has_skeletal)
    out["iso_skeletal_only"] = float(has_skeletal and not has_cardiac)

    exon = features.get("exon")
    out["in_exon"] = 1.0 if exon else 0.0
    out["exon_norm"] = (float(exon) / N_EXONS) if exon else 0.0

    try:
        pos = float(features.get("pos"))
        out["pos_norm"] = (pos - TTN_START) / (TTN_END - TTN_START)
    except (TypeError, ValueError):
        out["pos_norm"] = 0.0

    analogs: List[Dict[str, Any]] = features.get("analogs") or []
    out["n_analogs"] = float(len(analogs))
    if analogs:
        nearest = min(analogs, key=lambda a: a.get("distance", 10**9))
        d = float(nearest.get("distance") or 0.0)
        out["log_nearest_bp"] = math.log10(d + 1.0)
        out["nearest_is_cardiac"] = 1.0 if nearest.get("phenotype") == "cardiac" else 0.0

        # 只在明確標了 cardiac / skeletal 的類比之間算比例，both / unsure 不投票
        sided = [a for a in analogs if a.get("phenotype") in ("cardiac", "skeletal")]
        if sided:
            out["analog_cardiac_frac"] = sum(
                1.0 for a in sided if a["phenotype"] == "cardiac"
            ) / len(sided)
            wsum = card = 0.0
            for a in sided:
                w = 0.5 ** (float(a.get("distance") or 0.0) / ANALOG_HALF_LIFE_BP)
                wsum += w
                if a["phenotype"] == "cardiac":
                    card += w
            out["analog_cardiac_frac_weighted"] = (card / wsum) if wsum else 0.0
    else:
        out["log_nearest_bp"] = math.log10(1e6)

    return out


class TissueClassifier:
    """訓練好的組織判定器。predict_proba 回傳「是心肌」的機率。"""

    def __init__(self, model, feature_names: List[str], metadata: Optional[Dict] = None):
        self.model = model
        self.feature_names = list(feature_names)
        self.metadata = metadata or {}

    # ── 推論 ─────────────────────────────────────────────────────────────────

    def _vector(self, features: Dict[str, Any]) -> List[float]:
        f = tissue_features(features)
        return [f[name] for name in self.feature_names]

    def predict_cardiac_proba(self, features: Dict[str, Any]) -> float:
        import numpy as np

        x = np.asarray([self._vector(features)], dtype=float)
        return float(self.model.predict_proba(x)[0, 1])

    def predict_tissue(
        self, features: Dict[str, Any], uncertain_margin: float = 0.0
    ) -> Dict[str, Any]:
        """
        回傳 {"tissue": Cardiac|Skeletal|Uncertain, "cardiac_proba": p}。

        uncertain_margin > 0 時，機率落在 0.5 ± margin 的案例回 Uncertain，
        讓下游可以選擇不硬猜（預設 0，維持一定給答案的行為）。
        """
        p = self.predict_cardiac_proba(features)
        if uncertain_margin > 0 and abs(p - 0.5) < uncertain_margin:
            tissue = "Uncertain"
        else:
            tissue = "Cardiac" if p >= 0.5 else "Skeletal"
        return {"tissue": tissue, "cardiac_proba": p}

    # ── 存取 ─────────────────────────────────────────────────────────────────

    def save(self, path: Path) -> None:
        import joblib

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(
            {"model": self.model, "feature_names": self.feature_names,
             "metadata": self.metadata},
            path,
        )
        path.with_suffix(".json").write_text(
            json.dumps(self.metadata, ensure_ascii=False, indent=2), encoding="utf-8"
        )

    @classmethod
    def load(cls, path: Path) -> Optional["TissueClassifier"]:
        """載入失敗一律回 None，讓呼叫端安靜退回原本的規則式組織判定。"""
        path = Path(path)
        if not path.exists():
            return None
        try:
            import joblib

            blob = joblib.load(path)
            return cls(blob["model"], blob["feature_names"], blob.get("metadata"))
        except Exception as e:  # noqa: BLE001
            logger.warning(f"載入組織判定器失敗，改用規則式組織判定：{e}")
            return None
