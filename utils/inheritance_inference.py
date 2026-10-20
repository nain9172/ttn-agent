"""推論 TTN 變異的遺傳模式：體染色體顯性（AD）或隱性（AR）。

TTN 同一個基因兩種遺傳模式都有，而且主要由「變異類型 × 所在外顯子 × 疾病」決定：
  - 心臟高收錄外顯子（尤其 A-band）的截斷變異（TTNtv）：單一 allele 造成單倍劑量不足 → AD 擴張型心肌病變
  - 心臟不收錄（低 PSI、meta-only）外顯子的截斷變異：心臟不受影響，要兩個 allele 都壞才致病 → AR 骨骼肌病變
  - Mex1–Mex4（exon 359–362）的截斷：Salih myopathy 等 AR 早發肌病變
  - Mex6（exon 364）：異合子 TMD（AD）、同合子／複合異合子 LGMD2J（AR），同一顆變異兩種都有
  - exon 344（A150 FN3 domain）missense：HMERF（AD）
  - 其他 missense：多半以複合異合子（配另一顆截斷）出現在隱性肌病變，或證據不足

輸出的是「推論」，報告會標成 inference，不是資料庫紀錄。證據分四層：
  1. 這顆變異的紀錄：HGMD 精確命中的 comments（homozygous / compound het）與疾病、文獻擷取的 inheritance
  2. 疾病對應：預測或紀錄的疾病本身的遺傳模式（OMIM）
  3. 機轉規則：上面列的變異類型 × 外顯子規則
  4. 類比：同外顯子／同 band、同類型的 HGMD 已知變異的遺傳模式
有 LLM 時把這些證據交給 LLM 做最後判斷；沒有時用規則分數。validate() 對 HGMD 做 leave-one-out。

注意 HGMD 的 expected_inheritance 欄位全部都是 AD（基因層級的預設），不能拿來當變異層級的答案；
這裡的 HGMD 標籤改由 comments 的 zygosity 與疾病推得。
"""

from __future__ import annotations

import json
import os
import logging
import re
from collections import Counter
from typing import Any, Dict, Iterable, List, Optional

logger = logging.getLogger(__name__)

AD, AR, BOTH, UNCERTAIN = "Autosomal dominant", "Autosomal recessive", "AD or AR (allelic)", "Uncertain"
MODES = (AD, AR, BOTH, UNCERTAIN)
SHORT = {AD: "AD", AR: "AR", BOTH: "AD/AR", UNCERTAIN: "?"}

# 疾病 → 遺傳模式（OMIM 對 TTN 表型的登記；比對時不分大小寫，先比 AR 再比 AD）
_AR_DISEASE = re.compile(
    r"limb.?girdle|LGMD|salih|centronuclear|\bCNM\b|congenital myopathy|myopathy, congenital|arthrogryposis|"
    r"fetal akinesia|multi.?minicore|core myopathy|early.?onset myopathy|EOMFC|autosomal recessive|"
    r"emery.?dreifuss|prenatal|recessive",
    re.IGNORECASE,
)
_AD_DISEASE = re.compile(
    r"dilated|\bDCM\b|hypertrophic|\bHCM\b|arrhythmogenic|ARVC|non.?compaction|LVNC|atrial fibrillation|"
    r"peripartum|tibial muscular dystrophy|\bTMD\b|udd|early respiratory failure|HMERF|"
    r"cardiomyopathy|sudden (?:unexpected|cardiac)? ?death|left ventricular dysfunction|heart failure",
    re.IGNORECASE,
)
# TMD 與 LGMD2J 是同一組變異的異合子／同合子表現
_ALLELIC_DISEASE = re.compile(r"(?=.*(?:tibial|\bTMD\b|udd))(?=.*(?:LGMD|limb.?girdle))", re.IGNORECASE)

_AR_ZYGOSITY = re.compile(r"homozyg|compound het|comp\.? ?het|compound-het|biallelic|in trans", re.IGNORECASE)
_HET = re.compile(r"heterozyg|\bhet\b", re.IGNORECASE)

MEX_AR_EXONS = range(359, 363)  # Mex1–Mex4：截斷 → Salih myopathy 等 AR
MEX6_EXON = 364
HMERF_EXON = 344
PSI_HIGH = 90.0
MECH_LOW_PSI_WEIGHT = 1.5


def disease_inheritance(disease: Any) -> Optional[str]:
    """疾病名稱對應的遺傳模式；對不到回 None。"""
    text = str(disease or "")
    if not text.strip():
        return None
    if _ALLELIC_DISEASE.search(text):
        return BOTH
    if _AR_DISEASE.search(text):
        return AR
    if _AD_DISEASE.search(text):
        return AD
    return None


def record_label(disease: Any, comment: Any) -> Optional[str]:
    """一筆 HGMD 紀錄的遺傳模式：comments 寫同合子／複合異合子 → AR，否則看疾病。"""
    if _AR_ZYGOSITY.search(str(comment or "")):
        return AR
    return disease_inheritance(disease)


def normalize_mode(text: Any) -> str:
    t = str(text or "").strip().lower()
    if not t or t in ("not specified", "n/a", "none", "unknown", "uncertain"):
        return UNCERTAIN
    has_ad = "dominant" in t or re.search(r"\bad\b", t) is not None
    has_ar = "recessive" in t or re.search(r"\bar\b", t) is not None or "compound" in t or "homozyg" in t
    if has_ad and has_ar:
        return BOTH
    if has_ad:
        return AD
    if has_ar:
        return AR
    return UNCERTAIN


class InheritanceInference:
    INHERITANCE_JSON_SCHEMA: Dict[str, Any] = {
        "type": "object",
        "properties": {
            "inheritance": {"type": "string", "enum": list(MODES)},
            "confidence": {"type": "string", "enum": ["High", "Moderate", "Low"]},
            "rationale": {"type": "string"},
        },
        "required": ["inheritance", "confidence", "rationale"],
    }

    def __init__(self, mechanistic, hgmd_csv_path=None, n_analogs: int = 8):
        """mechanistic: utils.mechanistic_inference.MechanisticInference（共用外顯子／band／PSI 查表）。"""
        self.mi = mechanistic
        self.hgmd_csv_path = hgmd_csv_path or getattr(mechanistic, "hgmd_csv_path", None)
        self.n_analogs = n_analogs
        self._labelled: Optional[List[Dict[str, Any]]] = None

    # ── HGMD 標籤表 ─────────────────────────────────────────────────────────

    def labelled_hgmd(self) -> List[Dict[str, Any]]:
        """HGMD 每筆紀錄加上由 comments／疾病推得的遺傳模式（推不出來的不收）。"""
        if self._labelled is not None:
            return self._labelled
        self._labelled = []
        if not self.hgmd_csv_path or not os.path.exists(str(self.hgmd_csv_path)):
            logger.warning(f"HGMD CSV not found, inheritance uses rules only: {self.hgmd_csv_path}")
            return self._labelled
        import pandas as pd

        from utils.mechanistic_inference import _classify_variant, _text

        mutype_map = getattr(self.mi, "_MUTYPE_MAP", {})
        df = pd.read_csv(self.hgmd_csv_path, low_memory=False)
        for r in df.itertuples():
            try:
                pos = int(r.startCoord)
            except (TypeError, ValueError):
                continue
            disease = str(getattr(r, "disease", "") or "")
            comment = getattr(r, "comments", "")
            comment = "" if pd.isna(comment) else str(comment)
            label = record_label(disease, comment)
            if label is None:
                continue
            vclass, vgroup = mutype_map.get(str(getattr(r, "mutype", "") or "").strip(), (None, None))
            if vgroup is None:
                c = _classify_variant(str(getattr(r, "amino", "") or ""), str(getattr(r, "hgvs", "") or ""))
                vclass, vgroup = c["vclass"], c["vgroup"]
            self._labelled.append({
                "pos": pos,
                "exon": self.exon_for(pos, self.mi.assign_exon(pos)),
                "band": self.mi.assign_band(pos),
                "vclass": vclass, "vgroup": vgroup,
                "disease": disease,
                "label": label,
                "zygosity_ar": bool(_AR_ZYGOSITY.search(comment)),
                "prot": _text(getattr(r, "amino", "")),
                "dna": _text(getattr(r, "hgvs", "")),
                "acc": str(getattr(r, "acc_num", "") or ""),
            })
        logger.info(f"Inheritance: {len(self._labelled)} HGMD records with an inferable mode of inheritance")
        return self._labelled

    def _psi(self, exon: Optional[int]) -> Optional[Dict[str, Any]]:
        """實測 PSI：論文補充表（323 個外顯子）優先，缺的（exon 171–204 等）改用 cardiodb。"""
        if exon is None:
            return None
        try:
            hit = self.mi._load_exon_psi().get(int(exon))
            if hit:
                return hit
            for r in self.mi._load_cardio():
                if r["exon"] == int(exon):
                    return {"psi_heart": float(r["psi_heart"]), "psi_skeletal": float(r["psi_skeletal"]),
                            "meta_only": 0}
        except Exception:  # noqa: BLE001
            pass
        return None

    def exon_for(self, pos: int, exon: Optional[int] = None, max_dist: int = 20) -> Optional[int]:
        """變異所在外顯子；剪接位點等內含子變異取 max_dist bp 內最近的外顯子。"""
        if exon is not None:
            return exon
        best = None
        for lo, hi, no in self.mi._load_exons():
            d = 0 if lo <= pos <= hi else min(abs(pos - lo), abs(pos - hi))
            if d <= max_dist and (best is None or d < best[0]):
                best = (d, no)
        return best[1] if best else None

    # ── 證據 ────────────────────────────────────────────────────────────────

    def mechanism_rule(self, exon: Optional[int], band: str, vgroup: str, vclass: str = "") -> Dict[str, Any]:
        """變異類型 × 外顯子的機轉規則（TTN 文獻共識）。回傳 {mode, weight, reason}。"""
        psi = self._psi(exon)
        heart = psi["psi_heart"] if psi else None
        ex = f"exon {exon}" if exon else "an intron / unmapped position"
        if exon == MEX6_EXON:
            return {"mode": BOTH, "weight": 2.0,
                    "reason": "Mex6 (exon 364): heterozygous variants cause tibial muscular dystrophy (AD); "
                              "homozygous / compound-heterozygous ones cause LGMD2J (AR)"}
        if vgroup == "truncating":
            if exon in MEX_AR_EXONS:
                return {"mode": BOTH, "weight": 1.0,
                        "reason": f"truncating variant in {ex} (Mex1–Mex4, M-band): such variants are reported both "
                                  "heterozygous in DCM (AD) and biallelic in recessive early-onset myopathies "
                                  "(e.g. Salih myopathy)"}
            if psi and (psi.get("meta_only") or (heart is not None and heart < PSI_HIGH)):
                return {"mode": AR, "weight": MECH_LOW_PSI_WEIGHT,
                        "reason": f"truncating variant in {ex}, which is spliced out of most cardiac transcripts "
                                  f"(heart PSI {heart:.0f}%): heterozygous carriers are not expected to develop "
                                  "DCM, so disease needs a second allele (recessive myopathy)"}
            if band == "A-band":
                return {"mode": AD, "weight": 2.0,
                        "reason": f"truncating variant in {ex} (A-band, constitutively expressed in heart): "
                                  "one allele causes haploinsufficiency, the established AD mechanism of TTN DCM"}
            if heart is not None:
                return {"mode": AD, "weight": 1.2,
                        "reason": f"truncating variant in {ex} ({band}, heart PSI {heart:.0f}%): highly expressed "
                                  "in cardiac isoforms, so a single allele can cause AD cardiomyopathy "
                                  "(weaker enrichment outside the A-band)"}
            return {"mode": UNCERTAIN, "weight": 0.0,
                    "reason": f"truncating variant in {ex}; exon usage unknown"}
        if vgroup == "missense":
            if exon == HMERF_EXON:
                return {"mode": AD, "weight": 2.0,
                        "reason": "missense variant in exon 344 (A150 fibronectin-3 domain): the HMERF hotspot, "
                                  "inherited as AD"}
            return {"mode": AR, "weight": 0.5,
                    "reason": "TTN missense variants are rarely disease-causing alone; when pathogenic they are "
                              "usually found in trans with a truncating variant (recessive myopathy)"}
        return {"mode": UNCERTAIN, "weight": 0.0, "reason": "variant type unknown"}

    def analogs(self, pos: int, exon: Optional[int], band: str, vgroup: str,
                exclude: Iterable[str] = ()) -> List[Dict[str, Any]]:
        """同外顯子+同類型 → 同 band+同類型 → 同外顯子，依距離取最近的已標註 HGMD 變異。"""
        ex = set(exclude)
        rows = [r for r in self.labelled_hgmd() if r["acc"] not in ex and r["pos"] != pos]
        typed = vgroup not in ("unknown", None)
        tiers = [
            ("exon+type", lambda r: exon is not None and r["exon"] == exon and (not typed or r["vgroup"] == vgroup)),
            ("band+type", lambda r: r["band"] == band and (not typed or r["vgroup"] == vgroup)),
            ("exon", lambda r: exon is not None and r["exon"] == exon),
        ]
        out: List[Dict[str, Any]] = []
        seen = set()
        for tier, pred in tiers:
            hits = sorted((r for r in rows if r["acc"] not in seen and pred(r)), key=lambda r: abs(r["pos"] - pos))
            for r in hits[: self.n_analogs - len(out)]:
                out.append(dict(r, distance=abs(r["pos"] - pos), tier=tier))
                seen.add(r["acc"])
            if len(out) >= self.n_analogs:
                break
        return out

    def collect(
        self,
        features: Dict[str, Any],
        db_annotations: Optional[Any] = None,
        pubmed_results: Optional[List[Dict]] = None,
        predicted_disease: Optional[str] = None,
        exclude_acc: Iterable[str] = (),
    ) -> Dict[str, Any]:
        pos, band = int(features["pos"]), features.get("band") or ""
        exon = self.exon_for(pos, features.get("exon"))
        vgroup, vclass = features.get("vgroup") or "unknown", features.get("vclass") or "unknown"

        recorded = []
        data = {}
        if db_annotations is not None:
            data = db_annotations.to_dict() if hasattr(db_annotations, "to_dict") else dict(db_annotations)
        for r in data.get("hgmd_exact") or []:
            mode = record_label(r.get("disease"), r.get("comments"))
            if mode:
                recorded.append({"source": "HGMD", "disease": r.get("disease"), "mode": mode,
                                 "note": "homozygous / compound heterozygous" if _AR_ZYGOSITY.search(
                                     str(r.get("comments") or "")) else "from the recorded disease"})
        for r in data.get("pubtator_exact") or []:
            mode = disease_inheritance(r.get("disease"))
            if mode:
                recorded.append({"source": "PubTator", "disease": r.get("disease"), "mode": mode,
                                 "note": "from the recorded disease"})

        literature = []
        for a in pubmed_results or []:
            info = a.get("clinical_info") or {}
            mode = normalize_mode(info.get("inheritance"))
            if mode != UNCERTAIN:
                literature.append({"pmid": a.get("pmid"), "mode": mode, "text": info.get("inheritance")})

        disease_mode = disease_inheritance(predicted_disease) if predicted_disease else None
        return {
            "pos": pos, "exon": exon, "band": band, "vgroup": vgroup, "vclass": vclass,
            "psi": self._psi(exon),
            "recorded": recorded,
            "literature": literature,
            "predicted_disease": predicted_disease,
            "disease_mode": disease_mode,
            "mechanism": self.mechanism_rule(exon, band, vgroup, vclass),
            "analogs": self.analogs(pos, exon, band, vgroup, exclude_acc),
        }

    # ── 規則式判斷 ──────────────────────────────────────────────────────────

    @staticmethod
    def rule_predict(ev: Dict[str, Any]) -> Dict[str, Any]:
        score = Counter()
        reasons: List[str] = []

        def add(mode: str, w: float) -> None:
            if mode == BOTH:
                score[AD] += w / 2
                score[AR] += w / 2
                score[BOTH] += w
            elif mode in (AD, AR):
                score[mode] += w

        for r in ev["recorded"]:
            add(r["mode"], 4.0)
        if ev["recorded"]:
            reasons.append("database record for this variant: " + "; ".join(
                f"{r['source']} {r['disease']} → {SHORT[r['mode']]} ({r['note']})" for r in ev["recorded"]))
        for r in ev["literature"]:
            add(r["mode"], 2.0)
        if ev["literature"]:
            reasons.append("inheritance extracted from articles about this variant: " + ", ".join(
                f"PMID {r['pmid']} {SHORT[r['mode']]}" for r in ev["literature"]))
        if ev["disease_mode"]:
            add(ev["disease_mode"], 1.5)
            reasons.append(f"predicted disease {ev['predicted_disease']} is typically {SHORT[ev['disease_mode']]}")
        mech = ev["mechanism"]
        if mech["mode"] != UNCERTAIN:
            add(mech["mode"], mech["weight"])
            reasons.append(mech["reason"])
        votes = Counter()
        for a in ev["analogs"]:
            w = {"exon+type": 1.0, "band+type": 0.5, "exon": 0.4}[a["tier"]] / (1 + a["distance"] / 20000)
            votes[a["label"]] += w
        total = sum(votes.values())
        if total:
            for mode, w in votes.items():
                add(mode, 2.0 * w / total)
            reasons.append("nearest known HGMD variants (same exon / band and type): " + ", ".join(
                f"{SHORT[m]} {sum(1 for a in ev['analogs'] if a['label'] == m)}" for m in (AD, AR, BOTH)
                if any(a["label"] == m for a in ev["analogs"])))

        ad, ar = score[AD], score[AR]
        if ad + ar == 0:
            return {"inheritance": UNCERTAIN, "confidence": "Low", "reasons": reasons or ["no usable evidence"],
                    "scores": {"AD": 0.0, "AR": 0.0}}
        share = max(ad, ar) / (ad + ar)
        if score[BOTH] >= 2.0 and share < 0.75:
            mode = BOTH
        elif share < 0.6:
            mode = BOTH if min(ad, ar) >= 1.5 else UNCERTAIN
        else:
            mode = AD if ad > ar else AR
        strength = ad + ar
        if ev["recorded"] and mode != UNCERTAIN:
            conf = "High"
        elif share >= 0.8 and strength >= 3:
            conf = "Moderate"
        else:
            conf = "Low"
        return {"inheritance": mode, "confidence": conf, "reasons": reasons,
                "scores": {"AD": round(ad, 2), "AR": round(ar, 2)}}

    # ── LLM ─────────────────────────────────────────────────────────────────

    def build_prompt(self, ev: Dict[str, Any], rule: Dict[str, Any], is_medgemma: bool = True) -> str:
        psi = ev.get("psi")
        psi_line = (f"heart PSI {psi['psi_heart']:.0f}%, skeletal PSI {psi['psi_skeletal']:.0f}%"
                    + (", meta-transcript only" if psi.get("meta_only") else "")) if psi else "not available"
        recorded = "\n".join(f"- {r['source']}: {r['disease']} ({r['note']})" for r in ev["recorded"]) or "- none"
        lit = "\n".join(f"- PMID {r['pmid']}: {r['text']}" for r in ev["literature"]) or "- none"
        analogs = "\n".join(
            f"- {a['prot'] or a['dna']} | {a['vclass']} | exon {a['exon']} | {a['disease']} | "
            f"{'homozygous/compound het' if a['zygosity_ar'] else 'zygosity not stated'} | "
            f"mode {SHORT[a['label']]} | {a['distance']} bp away"
            for a in ev["analogs"]
        ) or "- none"
        system = f"""You are a clinical geneticist specialising in titin (TTN). Decide the most likely mode of
inheritance for a single TTN variant: autosomal dominant (one allele is enough) or autosomal recessive
(biallelic: homozygous or compound heterozygous), or both (allelic disorders such as TMD/LGMD2J).

Key TTN principles:
- Truncating variants in exons constitutively expressed in the heart (especially the A-band) cause
  dilated cardiomyopathy by haploinsufficiency: autosomal dominant with incomplete penetrance.
- Truncating variants in exons spliced out of cardiac isoforms (low heart PSI, meta-only exons) do not
  cause DCM in carriers; they cause recessive skeletal/congenital myopathies when biallelic.
- Truncations in Mex1-Mex4 (exons 359-362) cause recessive early-onset myopathies (e.g. Salih myopathy).
- Mex6 (exon 364) variants: tibial muscular dystrophy when heterozygous (AD), LGMD2J when homozygous (AR).
- Missense variants in exon 344 (A150 FN3 domain) cause HMERF (AD). Other missense variants are rarely
  pathogenic alone and are mostly seen in trans with a truncating variant (AR) or are of uncertain effect.

## VARIANT
chr2:{ev['pos']} | type: {ev['vclass']} ({ev['vgroup']}) | band: {ev['band']} | exon: {ev['exon']} | {psi_line}
Predicted / recorded disease: {ev.get('predicted_disease') or 'not available'}

## RECORDS FOR THIS EXACT VARIANT
{recorded}

## INHERITANCE STATED IN ARTICLES ABOUT THIS VARIANT
{lit}

## NEAREST KNOWN HGMD VARIANTS (same exon / band and type)
{analogs}

## RULE-BASED ESTIMATE
{rule['inheritance']} ({rule['confidence']}); reasons: {'; '.join(rule['reasons'])}

## OUTPUT
Respond with ONLY a JSON object:
{{"inheritance": "<{'|'.join(MODES)}>", "confidence": "<High|Moderate|Low>",
  "rationale": "<2-3 sentences: the variant type and exon, the key evidence, and that this is an inference>"}}"""
        user = "Give the inheritance JSON for this variant. Output ONLY the JSON object."
        if is_medgemma:
            return f"<start_of_turn>user\n{system}\n\n{user}<end_of_turn>\n<start_of_turn>model\n```json"
        return (f"<|begin_of_text|><|start_header_id|>system<|end_header_id|>\n{system}<|eot_id|>"
                f"<|start_header_id|>user<|end_header_id|>\n{user}<|eot_id|>"
                f"<|start_header_id|>assistant<|end_header_id|>\n```json")

    def _sampling_params(self, is_medgemma: bool):
        try:
            from vllm import SamplingParams
        except ImportError:
            return None
        kwargs: Dict[str, Any] = {"n": 1, "temperature": 0.0, "max_tokens": 512,
                                  "stop": ["<end_of_turn>"] if is_medgemma else ["<|eot_id|>"]}
        try:
            from vllm.sampling_params import StructuredOutputsParams

            kwargs["structured_outputs"] = StructuredOutputsParams(json=self.INHERITANCE_JSON_SCHEMA)
        except Exception:  # noqa: BLE001
            pass
        return SamplingParams(**kwargs)

    def generate(self, extractor, prompts: List[str]) -> List[str]:
        is_medgemma = getattr(extractor, "is_medgemma", True)
        llm = getattr(extractor, "llm", None)
        params = self._sampling_params(is_medgemma) if llm is not None else None
        if llm is None or params is None:
            return [extractor._generate_vllm(p) for p in prompts]
        return [o.outputs[0].text for o in llm.generate(prompts, params)]

    @staticmethod
    def parse(raw: str) -> Optional[Dict[str, Any]]:
        text = (raw or "").strip().strip("`")
        m = re.search(r"\{.*\}", text, re.DOTALL)
        if not m:
            return None
        try:
            data = json.loads(m.group(0))
        except json.JSONDecodeError:
            return None
        mode = data.get("inheritance")
        if mode not in MODES:
            mode = normalize_mode(mode)
        return {"inheritance": mode, "confidence": data.get("confidence") or "Low",
                "rationale": str(data.get("rationale") or "").strip()}

    def infer(
        self,
        extractor,
        features: Dict[str, Any],
        db_annotations: Optional[Any] = None,
        pubmed_results: Optional[List[Dict]] = None,
        predicted_disease: Optional[str] = None,
    ) -> Dict[str, Any]:
        ev = self.collect(features, db_annotations, pubmed_results, predicted_disease)
        rule = self.rule_predict(ev)
        result = {
            "inheritance": rule["inheritance"], "confidence": rule["confidence"],
            "rationale": "Rule-based: " + "; ".join(rule["reasons"]) + ".",
            "source": "rule-based", "rule": rule, "evidence": _public_evidence(ev),
        }
        if extractor is None:
            return result
        is_medgemma = getattr(extractor, "is_medgemma", True)
        try:
            raw = self.generate(extractor, [self.build_prompt(ev, rule, is_medgemma)])[0]
        except Exception as e:  # noqa: BLE001
            logger.warning(f"Inheritance LLM call failed, keeping the rule-based estimate: {e}")
            result["source"] = "rule-based (LLM error)"
            return result
        parsed = self.parse(raw)
        if not parsed:
            result["source"] = "rule-based (LLM output unparsable)"
            return result
        result.update(parsed, source="LLM (" + str(getattr(extractor, "model_name", "") or "local") + ")")
        return result


def _public_evidence(ev: Dict[str, Any]) -> Dict[str, Any]:
    """報告用的精簡證據（不含整列 HGMD 資料）。"""
    return {
        "recorded": ev["recorded"],
        "literature": ev["literature"],
        "predicted_disease": ev.get("predicted_disease"),
        "disease_mode": ev.get("disease_mode"),
        "mechanism": ev["mechanism"]["reason"] if ev["mechanism"]["mode"] != UNCERTAIN else "",
        "analogs": [
            {"variant": a["prot"] or a["dna"], "disease": a["disease"], "mode": SHORT[a["label"]],
             "exon": a["exon"], "tier": a["tier"]}
            for a in ev["analogs"][:5]
        ],
    }
