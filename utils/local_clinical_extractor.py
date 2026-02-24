#!/usr/bin/env python3
"""
Improved Local Clinical Information Extractor (v2)

Key improvements over v1:
1. Multi-window context extraction — captures ALL alias occurrences, not just the first
2. Better prompt with few-shot examples and permissive disease extraction
3. Regex-based fallback when LLM fails or returns "Not specified"
4. Structured context with labeled sections (Abstract / Results / Tables)
5. Smarter aggregation with more phenotype categories
"""

import logging
import json
import re
import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from utils.variant_utils import get_variant_aliases

try:
    from config import ENABLE_EASY_PROMPT
except ImportError:
    ENABLE_EASY_PROMPT = False

logger = logging.getLogger(__name__)

# ── Disease keyword → canonical disease name ─────────────────────────────────
DISEASE_PATTERNS: List[Tuple[str, str, str]] = [
    # pattern (case-insensitive), canonical name, tissue
    (r'\bDCM\b|dilated cardiomyopathy',                'DCM (Dilated Cardiomyopathy)',   'Cardiac'),
    (r'\bHCM\b|hypertrophic cardiomyopathy',           'HCM (Hypertrophic Cardiomyopathy)', 'Cardiac'),
    (r'\bARVC\b|\bACM\b|arrhythmogenic.*cardiomyopathy', 'ARVC/ACM',                       'Cardiac'),
    (r'\bLVNC\b|left ventricular non-compaction',      'LVNC',                            'Cardiac'),
    (r'\bheart failure\b',                             'Heart Failure',                   'Cardiac'),
    (r'\batrial fibrillation\b|\bAF\b',                'Atrial Fibrillation',             'Cardiac'),
    (r'\bcardiomyopathy\b',                            'Cardiomyopathy (unspecified)',     'Cardiac'),
    (r'\bTMD\b|tibial muscular dystrophy',             'TMD (Tibial Muscular Dystrophy)', 'Skeletal'),
    (r'\bLGMD\b|limb.girdle muscular dystrophy',       'LGMD',                            'Skeletal'),
    (r'\bmuscular dystrophy\b',                        'Muscular Dystrophy',              'Skeletal'),
    (r'\bmyopathy\b',                                  'Myopathy',                        'Skeletal'),
    (r'\bUDD\b|udd myopathy',                          'UDD Myopathy',                   'Skeletal'),
    (r'\bcardioskeletal myopathy\b',                   'Cardioskeletal Myopathy',         'Both'),
]

INHERITANCE_PATTERNS = [
    (r'autosomal dominant|AD\b',    'Autosomal Dominant'),
    (r'autosomal recessive|AR\b',   'Autosomal Recessive'),
    (r'X-linked',                   'X-linked'),
    (r'de novo',                    'De Novo'),
    (r'\bsporadic\b',               'Sporadic'),
    (r'heterozygous',               'Autosomal Dominant'),  # weak signal
    (r'homozygous|compound het',    'Autosomal Recessive'),  # weak signal
]

AGE_PATTERNS = [
    (r'congenital|neonatal|infant|newborn|birth',                          'Congenital'),
    (r'pediatric|child|childhood|juvenile|\b[1-9] year|\b1[0-7] year',    'Pediatric'),
    (r'adult|\b[2-9]\d year|\b1[89] year',                                'Adult'),
]


class LocalClinicalExtractor:
    def __init__(
        self,
        backend: str = "vllm",
        model_name: str = "meta-llama/Llama-3.2-3B",
        tensor_parallel_size: int = 2,
        max_model_len: int = 32768,
        max_context_length: int = 24000,
    ):
        self.backend = backend
        self.model_name = model_name
        self.tensor_parallel_size = tensor_parallel_size
        self.max_model_len = max_model_len
        self.max_context_length = max_context_length

        self.is_medgemma = "medgemma" in model_name.lower() or "gemma" in model_name.lower()
        self.is_llama = "llama" in model_name.lower()

        logger.info(f"Initializing {backend} backend with model: {model_name}")

        if backend == "vllm":
            self._init_vllm()
        elif backend == "ollama":
            self._init_ollama()
        else:
            raise ValueError(f"Unsupported backend: {backend}")

    # ── Initialisation ────────────────────────────────────────────────────────

    def _init_vllm(self):
        from vllm import LLM, SamplingParams

        gpu_memory_utilization = 0.85
        stop_tokens = (
            ["<end_of_turn>", "\n\n\n"] if self.is_medgemma else ["\n\n\n", "```\n\n"]
        )

        self.llm = LLM(
            model=self.model_name,
            tensor_parallel_size=self.tensor_parallel_size,
            trust_remote_code=True,
            gpu_memory_utilization=gpu_memory_utilization,
            max_model_len=self.max_model_len,
            enforce_eager=True,
        )

        temperature = 0.2 if self.is_medgemma else 0.3  # lower = more deterministic

        self.sampling_params = SamplingParams(
            temperature=temperature,
            top_p=0.9,
            max_tokens=8192,
            stop=stop_tokens,
        )

    def _init_ollama(self):
        import ollama
        self.ollama_client = ollama

    # ── Logging ───────────────────────────────────────────────────────────────

    def _log_prompt_response(self, pmid: str, prompt: str, response: str, variant_id: str = None):
        log_base_dir = os.path.join(os.path.dirname(__file__), "..", "llm_logs")
        log_dir = os.path.join(log_base_dir, variant_id) if variant_id else log_base_dir
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, f"{pmid}.txt")

        with open(log_file, "w", encoding="utf-8") as f:
            f.write("=" * 80 + "\n")
            f.write(f"Timestamp: {datetime.now().isoformat()}\n")
            if variant_id:
                f.write(f"Variant ID: {variant_id}\n")
            f.write(f"PMID: {pmid}\n")
            f.write(f"Model: {self.model_name}\n")
            f.write(f"Backend: {self.backend}\n")
            f.write("=" * 80 + "\n\n")
            f.write("=" * 40 + " PROMPT " + "=" * 40 + "\n")
            f.write(prompt)
            f.write("\n\n")
            f.write("=" * 40 + " RESPONSE " + "=" * 38 + "\n")
            f.write(response)
            f.write("\n" + "=" * 80 + "\n")

        logger.info(f"Logged prompt/response to: {log_file}")

    # ── Main API ──────────────────────────────────────────────────────────────

    def batch_extract(
        self,
        articles: List[Dict],
        variant_info: Dict,
        clinvar_info: Optional[Dict] = None,
    ) -> List[Dict]:
        variant_id = (
            f"{variant_info['chrom']}-{variant_info['pos']}"
            f"-{variant_info['ref']}-{variant_info['alt']}"
        )
        logger.info(f"Extracting from {len(articles)} articles. Target: {variant_id}")

        enriched = []
        for article in articles:
            if not any([article.get("title"), article.get("abstract"), article.get("full_text")]):
                continue
            info = self._extract_single(article, variant_id, clinvar_info)
            article["clinical_info"] = info

            if info.get("disease") not in ("Not specified", None):
                logger.info(
                    f"  [+] PMID {article.get('pmid')}: "
                    f"{info.get('disease')} | {info.get('tissue_affected')} | "
                    f"source={info.get('extraction_source','llm')}"
                )
            else:
                logger.info(f"  [-] PMID {article.get('pmid')}: No relevant info found")

            enriched.append(article)
        return enriched

    def _extract_single(
        self, article: Dict, variant_id: str, clinvar_info: Optional[Dict]
    ) -> Dict:
        raw_text = (
            article.get("text_for_llm")
            or article.get("full_text")
            or article.get("abstract")
            or ""
        )
        aliases = get_variant_aliases(variant_id, clinvar_info)

        # Build structured context with labelled sections
        context = self._prepare_context(raw_text, aliases, article)

        if len(context.strip()) < 50:
            return self._get_empty_result("Text too short")

        # ── Step 1: Try LLM extraction ────────────────────────────────────────
        llm_result = None
        if not ENABLE_EASY_PROMPT:
            prompt = self._build_prompt(context, aliases, variant_id)
            pmid = article.get("pmid", "unknown")
            try:
                if self.backend == "vllm":
                    raw = self._generate_vllm(prompt)
                else:
                    raw = self._generate_ollama(prompt)

                self._log_prompt_response(pmid, prompt, raw, variant_id)
                data = self._robust_json_parse(raw)
                data["raw_llm_output"] = raw
                data["extraction_source"] = "llm"
                llm_result = data
            except Exception as e:
                logger.error(f"LLM extraction error for PMID {pmid}: {e}")

        # ── Step 2: Regex fallback (always run if LLM gave "Not specified") ──
        if llm_result is None or llm_result.get("disease") in ("Not specified", None, ""):
            regex_result = self._regex_extract(raw_text, aliases, article)
            if regex_result.get("disease") not in ("Not specified", None):
                # Merge: prefer LLM values where non-empty, fill gaps with regex
                if llm_result:
                    merged = llm_result.copy()
                    for key in ("disease", "tissue_affected", "age_onset", "inheritance"):
                        if merged.get(key) in ("Not specified", None, ""):
                            merged[key] = regex_result.get(key, "Not specified")
                    merged["extraction_source"] = "llm+regex"
                    return merged
                else:
                    return regex_result

        return llm_result if llm_result else self._get_empty_result("All methods failed")

    # ── Context preparation ───────────────────────────────────────────────────

    def _prepare_context(
        self, text: str, aliases: List[str], article: Dict
    ) -> str:
        """
        Build a structured context string:
        - [ABSTRACT] section from article.abstract (always included, up to 2 000 chars)
        - [RELEVANT SECTIONS] — windows around each alias occurrence in full text
        Caps total length at self.max_context_length.
        """
        ABSTRACT_CAP = 2_000
        WINDOW = 1_500  # chars around each hit

        parts: List[str] = []
        budget = self.max_context_length

        # Abstract first (most reliable summary)
        abstract = (article.get("abstract") or "").strip()
        if abstract:
            snip = abstract[:ABSTRACT_CAP]
            parts.append(f"[ABSTRACT]\n{snip}")
            budget -= len(snip)

        # Title
        title = (article.get("title") or "").strip()
        if title:
            parts.insert(0, f"[TITLE]\n{title}\n")
            budget -= len(title)

        if budget <= 0 or not text:
            return "\n\n".join(parts)

        # Find all alias hit positions in the full text
        text_lower = text.lower()
        hit_positions: List[int] = []
        for alias in aliases:
            alias_lower = alias.lower()
            idx = 0
            while True:
                pos = text_lower.find(alias_lower, idx)
                if pos == -1:
                    break
                hit_positions.append(pos)
                idx = pos + 1

        hit_positions = sorted(set(hit_positions))

        if not hit_positions:
            # No alias found — take beginning of text (likely abstract-only paper)
            snippet = text[:min(budget, 4_000)]
            parts.append(f"[FULL TEXT START]\n{snippet}")
            return "\n\n".join(parts)

        # Merge overlapping windows and collect unique text chunks
        windows: List[Tuple[int, int]] = []
        for pos in hit_positions:
            start = max(0, pos - WINDOW)
            end = min(len(text), pos + WINDOW)
            if windows and start <= windows[-1][1]:
                windows[-1] = (windows[-1][0], max(windows[-1][1], end))
            else:
                windows.append((start, end))

        chunks: List[str] = []
        for start, end in windows:
            chunk = text[start:end].strip()
            if chunk:
                prefix = "..." if start > 0 else ""
                suffix = "..." if end < len(text) else ""
                chunks.append(f"{prefix}{chunk}{suffix}")

        combined = "\n\n---\n\n".join(chunks)
        if len(combined) > budget:
            combined = combined[:budget]

        parts.append(f"[RELEVANT SECTIONS (around variant mentions)]\n{combined}")
        return "\n\n".join(parts)

    # ── Prompt building ───────────────────────────────────────────────────────

    def _build_prompt(self, context: str, aliases: List[str], variant_id: str) -> str:
        alias_list = "\n".join(f"  - {a}" for a in aliases[:20])  # cap list length

        system_message = f"""You are a clinical genetics expert extracting structured data from biomedical literature about TTN gene variants.

## TARGET VARIANT
Variant ID: {variant_id}
Known aliases (search for ANY of these in the text):
{alias_list}

## YOUR TASK
Read the article text below and answer: does this article provide clinical information (disease, phenotype, patient data) for the target variant listed above?

### IMPORTANT RULES
1. **The variant may appear only in a table row** — that is sufficient. Extract the disease studied in that table/cohort.
2. **Be permissive**: if the study enrolls patients with DCM and lists our variant in a supplementary table, the disease IS DCM.
3. **Do NOT require the variant to be the focus of the paper** — a single row in a patient table counts.
4. If the variant is not mentioned at all, output disease = "Not specified".
5. For tissue_affected, use ONLY: "Cardiac", "Skeletal", "Both", or "Not specified".

### FEW-SHOT EXAMPLES

Example 1 — variant in a cohort table:
Text: "Table 2. Pathogenic TTN variants in our DCM cohort. ... c.12345G>A (p.Gly4115Ser) | 2 patients | het | adult onset"
Output:
```json
{{
  "reasoning": "Variant p.Gly4115Ser appears in a DCM cohort table with 2 patients.",
  "disease": "DCM (Dilated Cardiomyopathy)",
  "tissue_affected": "Cardiac",
  "age_onset": "Adult",
  "inheritance": "Autosomal Dominant",
  "patient_count": 2,
  "evidence_sentence": "c.12345G>A (p.Gly4115Ser) | 2 patients | het | adult onset"
}}
```

Example 2 — variant in HCM paper:
Text: "Among 150 HCM patients screened, we identified TTN variant I35947N in one proband. Echocardiography showed septal hypertrophy."
Output:
```json
{{
  "reasoning": "Variant I35947N found in one HCM patient. Echocardiography confirms cardiac phenotype.",
  "disease": "HCM (Hypertrophic Cardiomyopathy)",
  "tissue_affected": "Cardiac",
  "age_onset": "Not specified",
  "inheritance": "Not specified",
  "patient_count": 1,
  "evidence_sentence": "we identified TTN variant I35947N in one proband"
}}
```

Example 3 — variant not in article:
Text: "We studied 30 LGMD patients and performed TTN exome sequencing. No rare variants were found in exons 1-100."
Output:
```json
{{
  "reasoning": "The target variant is not mentioned in this article.",
  "disease": "Not specified",
  "tissue_affected": "Not specified",
  "age_onset": "Not specified",
  "inheritance": "Not specified",
  "patient_count": 0,
  "evidence_sentence": "Not specified"
}}
```

## OUTPUT FORMAT
Respond with ONLY a JSON object (no extra text):
```json
{{
  "reasoning": "<1-2 sentences on what you found>",
  "disease": "<disease name or 'Not specified'>",
  "tissue_affected": "<Cardiac|Skeletal|Both|Not specified>",
  "age_onset": "<Congenital|Pediatric|Adult|Not specified>",
  "inheritance": "<Autosomal Dominant|Autosomal Recessive|X-linked|De Novo|Sporadic|Not specified>",
  "patient_count": <integer>,
  "evidence_sentence": "<direct quote or 'Not specified'>"
}}
```"""

        user_message = f"""## ARTICLE TEXT

{context}

---

Extract clinical information for the target variant. Output ONLY the JSON object."""

        if self.is_medgemma:
            return (
                f"<start_of_turn>user\n{system_message}\n\n{user_message}"
                f"<end_of_turn>\n<start_of_turn>model\n```json"
            )
        else:
            return (
                f"<|begin_of_text|><|start_header_id|>system<|end_header_id|>\n"
                f"{system_message}<|eot_id|>"
                f"<|start_header_id|>user<|end_header_id|>\n"
                f"{user_message}<|eot_id|>"
                f"<|start_header_id|>assistant<|end_header_id|>\n```json"
            )

    # ── Regex fallback ────────────────────────────────────────────────────────

    def _regex_extract(
        self, text: str, aliases: List[str], article: Dict
    ) -> Dict:
        """
        Fast regex-based extraction as fallback.
        Works best when the article is clearly about a specific disease.
        """
        combined_text = " ".join(
            filter(None, [article.get("title", ""), article.get("abstract", ""), text[:5000]])
        )

        # Check if variant is even mentioned
        variant_mentioned = any(
            alias.lower() in combined_text.lower() for alias in aliases
        )

        # If variant not in full text, still check title+abstract for context
        search_scope = combined_text if variant_mentioned else " ".join(
            filter(None, [article.get("title", ""), article.get("abstract", "")])
        )

        found_disease = "Not specified"
        found_tissue = "Not specified"

        for pattern, disease, tissue in DISEASE_PATTERNS:
            if re.search(pattern, search_scope, re.IGNORECASE):
                found_disease = disease
                found_tissue = tissue
                break  # take most specific match first

        # If variant not explicitly mentioned but abstract strongly implies it
        if not variant_mentioned and found_disease != "Not specified":
            # Only trust this if a co-author-level connection (pmid from clinvar, etc.)
            # Mark as low-confidence
            found_disease = f"{found_disease} (inferred from abstract)"

        found_inheritance = "Not specified"
        for pattern, inheritance in INHERITANCE_PATTERNS:
            if re.search(pattern, search_scope, re.IGNORECASE):
                found_inheritance = inheritance
                break

        found_age = "Not specified"
        for pattern, age in AGE_PATTERNS:
            if re.search(pattern, search_scope, re.IGNORECASE):
                found_age = age
                break

        # Try to count patients
        patient_count = 0
        pat_match = re.search(
            r"(\d+)\s*(?:patients?|cases?|individuals?|probands?|carriers?)",
            search_scope,
            re.IGNORECASE,
        )
        if pat_match:
            try:
                patient_count = int(pat_match.group(1))
                patient_count = min(patient_count, 10_000)  # sanity cap
            except ValueError:
                pass

        return {
            "reasoning": f"Regex fallback. Variant mentioned: {variant_mentioned}.",
            "disease": found_disease,
            "tissue_affected": found_tissue,
            "age_onset": found_age,
            "inheritance": found_inheritance,
            "patient_count": patient_count,
            "evidence_sentence": "Not specified",
            "extraction_source": "regex",
        }

    # ── LLM generation ────────────────────────────────────────────────────────

    def _generate_vllm(self, prompt: str) -> str:
        outputs = self.llm.generate([prompt], self.sampling_params)
        return outputs[0].outputs[0].text

    def _generate_ollama(self, prompt: str) -> str:
        return self.ollama_client.generate(
            model=self.model_name, prompt=prompt
        )["response"]

    # ── JSON parsing ──────────────────────────────────────────────────────────

    def _robust_json_parse(self, text: str) -> Dict:
        clean = text.strip()
        for prefix in ("```json", "```"):
            if clean.startswith(prefix):
                clean = clean[len(prefix):]
        if clean.endswith("```"):
            clean = clean[:-3]
        clean = clean.strip()

        start = clean.find("{")
        end = clean.rfind("}")
        if start != -1 and end != -1:
            json_str = clean[start : end + 1]
            for attempt in (json_str, re.sub(r",\s*([}\]])", r"\1", json_str)):
                try:
                    return json.loads(attempt)
                except json.JSONDecodeError:
                    pass

        logger.warning(f"JSON parse failed. Raw output (first 400): {text[:400]}")
        return self._get_empty_result("JSON Parsing Failed")

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _get_empty_result(self, reason: str) -> Dict:
        return {
            "reasoning": reason,
            "disease": "Not specified",
            "tissue_affected": "Not specified",
            "age_onset": "Not specified",
            "inheritance": "Not specified",
            "evidence_sentence": "N/A",
            "patient_count": 0,
            "extraction_source": "none",
        }

    # ── Statistics aggregation ────────────────────────────────────────────────

    def aggregate_stats(self, articles: List[Dict]) -> Dict:
        stats = {
            "total_articles": len(articles),
            "articles_with_extraction": 0,
            "total_patients": 0,
            "total_families": 0,
            "disease_distribution": {},
            "inheritance_distribution": {},
            "age_onset_distribution": {},
            "tissue_affected": {"cardiac": 0, "skeletal_muscle": 0, "both": 0},
            "cardiac_phenotypes": {},
            "skeletal_phenotypes": {},
            "extraction_source_counts": {"llm": 0, "regex": 0, "llm+regex": 0, "none": 0},
        }

        for a in articles:
            info = a.get("clinical_info", {})
            disease = info.get("disease", "Not specified")
            source = info.get("extraction_source", "none")
            stats["extraction_source_counts"][source] = (
                stats["extraction_source_counts"].get(source, 0) + 1
            )

            skip_conditions = ("Not specified", None, "JSON Parsing Failed", "All methods failed")
            if disease in skip_conditions or (
                isinstance(disease, str) and "(inferred from abstract)" in disease
            ):
                continue

            stats["articles_with_extraction"] += 1
            try:
                stats["total_patients"] += int(info.get("patient_count", 0))
            except (ValueError, TypeError):
                pass

            for key, dist_key in (
                ("disease", "disease_distribution"),
                ("inheritance", "inheritance_distribution"),
                ("age_onset", "age_onset_distribution"),
            ):
                val = info.get(key, "Unknown") or "Unknown"
                if val not in ("Not specified", None):
                    stats[dist_key][val] = stats[dist_key].get(val, 0) + 1

            tissue = (info.get("tissue_affected") or "").lower()
            if "cardiac" in tissue:
                stats["tissue_affected"]["cardiac"] += 1
                d = info.get("disease", "Unknown")
                stats["cardiac_phenotypes"][d] = stats["cardiac_phenotypes"].get(d, 0) + 1
            if "skeletal" in tissue:
                stats["tissue_affected"]["skeletal_muscle"] += 1
                d = info.get("disease", "Unknown")
                stats["skeletal_phenotypes"][d] = stats["skeletal_phenotypes"].get(d, 0) + 1
            if "both" in tissue:
                stats["tissue_affected"]["both"] += 1

        return stats