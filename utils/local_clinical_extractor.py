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
from collections import Counter
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple
from utils.variant_utils import get_variant_aliases

try:
    from config import ENABLE_EASY_PROMPT
except ImportError:
    ENABLE_EASY_PROMPT = False

try:
    from config import (
        LLM_SELF_CONSISTENCY_N,
        LLM_SAMPLING_TEMPERATURE,
        LLM_USE_GUIDED_JSON,
        LLM_VERIFY_EVIDENCE,
        LLM_EVIDENCE_NGRAM_OVERLAP_THRESHOLD,
    )
except ImportError:
    LLM_SELF_CONSISTENCY_N = 1
    LLM_SAMPLING_TEMPERATURE = 0.2
    LLM_USE_GUIDED_JSON = False
    LLM_VERIFY_EVIDENCE = False
    LLM_EVIDENCE_NGRAM_OVERLAP_THRESHOLD = 0.5

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

        # B1: self-consistency — n>1 才有意義（要 sample 多個再投票）
        self.n_samples = max(1, int(LLM_SELF_CONSISTENCY_N))
        # n=1 時保留低溫；n>1 時拉高溫度才有 diversity
        if self.n_samples == 1:
            temperature = 0.2 if self.is_medgemma else 0.3
        else:
            temperature = float(LLM_SAMPLING_TEMPERATURE)

        sampling_kwargs: Dict[str, Any] = dict(
            n=self.n_samples,
            temperature=temperature,
            top_p=0.9,
            max_tokens=8192,
            stop=stop_tokens,
        )

        # B3: guided JSON decoding — 強制輸出符合 schema
        self.use_guided_json = bool(LLM_USE_GUIDED_JSON)
        if self.use_guided_json:
            try:
                from vllm.sampling_params import StructuredOutputsParams

                schema = self._build_json_schema()
                sampling_kwargs["structured_outputs"] = StructuredOutputsParams(
                    json=schema
                )
                logger.info("Guided JSON decoding enabled (vLLM StructuredOutputsParams)")
            except Exception as e:
                logger.warning(
                    f"Could not enable guided JSON decoding "
                    f"(vLLM={getattr(__import__('vllm'), '__version__', '?')}): {e}. "
                    "Falling back to plain sampling."
                )
                self.use_guided_json = False

        self.sampling_params = SamplingParams(**sampling_kwargs)
        logger.info(
            f"Sampling: n={self.n_samples} temperature={temperature} "
            f"guided_json={self.use_guided_json}"
        )

    @staticmethod
    def _build_json_schema() -> Dict[str, Any]:
        """JSON schema for the per-article clinical extraction output (B3)."""
        return {
            "type": "object",
            "properties": {
                "reasoning": {"type": "string"},
                "disease": {"type": "string"},
                "tissue_affected": {
                    "type": "string",
                    "enum": ["Cardiac", "Skeletal", "Both", "Not specified"],
                },
                "age_onset": {
                    "type": "string",
                    "enum": ["Congenital", "Pediatric", "Adult", "Not specified"],
                },
                "inheritance": {
                    "type": "string",
                    "enum": [
                        "Autosomal Dominant",
                        "Autosomal Recessive",
                        "X-linked",
                        "De Novo",
                        "Sporadic",
                        "Not specified",
                    ],
                },
                "patient_count": {"type": "integer", "minimum": 0},
                "evidence_sentence": {"type": "string"},
            },
            "required": [
                "reasoning",
                "disease",
                "tissue_affected",
                "age_onset",
                "inheritance",
                "patient_count",
                "evidence_sentence",
            ],
            "additionalProperties": False,
        }

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

        aliases = get_variant_aliases(variant_id, clinvar_info)

        # Filter out completely empty articles early
        valid_articles = [
            a for a in articles
            if any([a.get("title"), a.get("abstract"), a.get("full_text")])
        ]

        # ── Phase 1: Build contexts and prompts for all articles ──────────────
        # Each entry: (article, raw_text, context, prompt_or_None)
        pipeline: List[tuple] = []
        for article in valid_articles:
            raw_text = (
                article.get("text_for_llm")
                or article.get("full_text")
                or article.get("abstract")
                or ""
            )
            context = self._prepare_context(raw_text, aliases, article)
            if len(context.strip()) < 50:
                pipeline.append((article, raw_text, context, None))
                continue

            if not ENABLE_EASY_PROMPT and self.backend == "vllm":
                prompt = self._build_prompt(context, aliases, variant_id)
            else:
                prompt = None
            pipeline.append((article, raw_text, context, prompt))

        # ── Phase 2: Batch vLLM inference (one call for all articles) ─────────
        # Only used when backend == "vllm" and ENABLE_EASY_PROMPT is False
        prompt_indices: List[int] = []
        prompts_to_run: List[str] = []
        for i, (_, _, _, prompt) in enumerate(pipeline):
            if prompt is not None:
                prompt_indices.append(i)
                prompts_to_run.append(prompt)

        # batch_outputs[i] 是 list[str]（n 個 sample），對應 prompt_indices[i]
        batch_outputs: List[List[str]] = []
        if prompts_to_run:
            logger.info(
                f"  Running vLLM batch inference on {len(prompts_to_run)} prompts "
                f"(n_samples={self.n_samples}, guided_json={self.use_guided_json})..."
            )
            try:
                raw_outputs = self.llm.generate(prompts_to_run, self.sampling_params)
                batch_outputs = [
                    [c.text for c in o.outputs] for o in raw_outputs
                ]
                logger.info("  vLLM batch inference complete.")
            except Exception as e:
                logger.error(f"vLLM batch inference error: {e}")
                batch_outputs = [[""] for _ in prompts_to_run]

        # Map batch output index → pipeline index
        output_map: Dict[int, List[str]] = {
            pipeline_idx: batch_outputs[out_idx]
            for out_idx, pipeline_idx in enumerate(prompt_indices)
        }

        # ── Phase 3: Parse results, apply regex fallback, log ─────────────────
        enriched: List[Dict] = []
        for i, (article, raw_text, context, prompt) in enumerate(pipeline):
            pmid = article.get("pmid", "unknown")

            if prompt is None:
                # Text too short or easy-prompt mode: skip LLM
                llm_result = None
            else:
                raw_samples = output_map.get(i, [""])
                # Log every sample joined by a separator so all attempts are visible
                self._log_prompt_response(
                    pmid,
                    prompt,
                    "\n\n--- SAMPLE BOUNDARY ---\n\n".join(raw_samples),
                    variant_id,
                )
                # B1: parse all n samples then majority-vote
                parsed_samples: List[Dict] = []
                for raw in raw_samples:
                    try:
                        parsed_samples.append(self._robust_json_parse(raw))
                    except Exception as e:
                        logger.debug(f"Parse error for PMID {pmid} sample: {e}")
                if parsed_samples:
                    voted = self._majority_vote(parsed_samples, raw_text, article)
                    voted["raw_llm_output"] = raw_samples[0]
                    voted["extraction_source"] = (
                        f"llm_vote(n={len(parsed_samples)}/{len(raw_samples)})"
                    )
                    llm_result = voted
                else:
                    logger.error(f"All {len(raw_samples)} samples failed JSON parse for PMID {pmid}")
                    llm_result = None

            # Ollama fallback (non-batched path, used only when backend != vllm)
            if self.backend != "vllm" and not ENABLE_EASY_PROMPT and prompt is not None:
                try:
                    raw = self._generate_ollama(prompt)
                    self._log_prompt_response(pmid, prompt, raw, variant_id)
                    data = self._robust_json_parse(raw)
                    data["raw_llm_output"] = raw
                    data["extraction_source"] = "llm"
                    llm_result = data
                except Exception as e:
                    logger.error(f"LLM extraction error for PMID {pmid}: {e}")
                    llm_result = None

            # Regex fallback if LLM gave "Not specified" or failed
            if llm_result is None or llm_result.get("disease") in ("Not specified", None, ""):
                regex_result = self._regex_extract(raw_text, aliases, article)
                if regex_result.get("disease") not in ("Not specified", None):
                    if llm_result:
                        merged = llm_result.copy()
                        for key in ("disease", "tissue_affected", "age_onset", "inheritance"):
                            if merged.get(key) in ("Not specified", None, ""):
                                merged[key] = regex_result.get(key, "Not specified")
                        merged["extraction_source"] = (
                            f"{merged.get('extraction_source', 'llm')}+regex"
                        )
                        info = merged
                    else:
                        info = regex_result
                else:
                    info = llm_result if llm_result else self._get_empty_result("All methods failed")
            else:
                info = llm_result

            # B4: evidence grounding — 確認 evidence_sentence 真的出自原文
            if (
                LLM_VERIFY_EVIDENCE
                and info.get("disease") not in ("Not specified", None, "")
                and not (info.get("extraction_source", "").startswith("regex"))
            ):
                grounded, score = self._verify_evidence_grounding(
                    info.get("evidence_sentence", ""), raw_text, article
                )
                info["evidence_grounded"] = grounded
                info["evidence_overlap_score"] = round(score, 3)
                if not grounded:
                    logger.warning(
                        f"  PMID {pmid}: evidence NOT grounded (overlap={score:.2f}). "
                        f"Downgrading disease to 'Not specified'."
                    )
                    info["reasoning"] = (
                        f"[evidence ungrounded, overlap={score:.2f}] "
                        f"{info.get('reasoning', '')}"
                    )
                    info["disease"] = "Not specified"
                    info["tissue_affected"] = "Not specified"
                    info["extraction_source"] = (
                        f"{info.get('extraction_source', 'llm')}+ungrounded"
                    )

            # Handle text-too-short case
            if prompt is None and not ENABLE_EASY_PROMPT and len(context.strip()) < 50:
                info = self._get_empty_result("Text too short")

            article["clinical_info"] = info

            if info.get("disease") not in ("Not specified", None):
                logger.info(
                    f"  [+] PMID {pmid}: "
                    f"{info.get('disease')} | {info.get('tissue_affected')} | "
                    f"source={info.get('extraction_source', 'llm')}"
                )
            else:
                logger.info(f"  [-] PMID {pmid}: No relevant info found")

            enriched.append(article)
        return enriched

    def _extract_single(
        self, article: Dict, variant_id: str, clinvar_info: Optional[Dict]
    ) -> Dict:
        """Single-article extraction (kept for backward compatibility; not used by batch_extract)."""
        raw_text = (
            article.get("text_for_llm")
            or article.get("full_text")
            or article.get("abstract")
            or ""
        )
        aliases = get_variant_aliases(variant_id, clinvar_info)

        context = self._prepare_context(raw_text, aliases, article)

        if len(context.strip()) < 50:
            return self._get_empty_result("Text too short")

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

        if llm_result is None or llm_result.get("disease") in ("Not specified", None, ""):
            regex_result = self._regex_extract(raw_text, aliases, article)
            if regex_result.get("disease") not in ("Not specified", None):
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
        - [TITLE] + [ABSTRACT] (always included)
        - [RELEVANT SECTIONS] — large windows around each alias occurrence
        - [FULL TEXT START] — large prefix when no alias is matched
        Caps total length at self.max_context_length (in chars).

        Tuned for DGX Spark / MedGemma-27B with n=5 self-consistency:
        - WINDOW = 5000 chars/hit so an entire variant table row + headers
          (typically 1–3 KB) is always inside the window
        - When no alias is matched, take a much larger prefix because some
          papers only mention the variant deep in supplementary text
        """
        ABSTRACT_CAP = 3_500
        WINDOW = 5_000  # chars around each hit (was 1_500; too small for tables)
        NO_HIT_FALLBACK = 25_000  # when zero alias hits, scan more of the paper

        parts: List[str] = []
        budget = self.max_context_length
        pmid = article.get("pmid", "?")

        # Title
        title = (article.get("title") or "").strip()
        if title:
            parts.append(f"[TITLE]\n{title}")
            budget -= len(title) + 10

        # Abstract (always included — most reliable summary)
        abstract = (article.get("abstract") or "").strip()
        if abstract:
            snip = abstract[:ABSTRACT_CAP]
            parts.append(f"[ABSTRACT]\n{snip}")
            budget -= len(snip) + 12

        if budget <= 0 or not text:
            ctx = "\n\n".join(parts)
            logger.info(
                f"  [ctx] PMID {pmid}: title+abstract only "
                f"({len(ctx)} chars, no full text)"
            )
            return ctx

        # Find all alias hit positions in the full text
        text_lower = text.lower()
        hit_positions: List[int] = []
        per_alias_hits: Dict[str, int] = {}
        for alias in aliases:
            alias_lower = alias.lower().strip()
            if not alias_lower or len(alias_lower) < 3:
                # Skip ultra-short aliases that match too aggressively (e.g. "AD")
                continue
            idx = 0
            count = 0
            while True:
                pos = text_lower.find(alias_lower, idx)
                if pos == -1:
                    break
                hit_positions.append(pos)
                idx = pos + 1
                count += 1
            if count:
                per_alias_hits[alias] = count

        hit_positions = sorted(set(hit_positions))

        if not hit_positions:
            # No alias found — take a generous prefix of the full text
            snippet = text[: min(budget, NO_HIT_FALLBACK)]
            parts.append(f"[FULL TEXT START]\n{snippet}")
            ctx = "\n\n".join(parts)
            logger.info(
                f"  [ctx] PMID {pmid}: NO alias hits in {len(text)} chars of text "
                f"-> sending {len(snippet)} char prefix (total ctx={len(ctx)})"
            )
            return ctx

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
        truncated = False
        if len(combined) > budget:
            combined = combined[:budget]
            truncated = True

        parts.append(f"[RELEVANT SECTIONS (around variant mentions)]\n{combined}")
        ctx = "\n\n".join(parts)

        top_aliases = sorted(per_alias_hits.items(), key=lambda kv: -kv[1])[:5]
        logger.info(
            f"  [ctx] PMID {pmid}: {len(hit_positions)} alias hits "
            f"({len(windows)} merged windows), context={len(ctx)} chars"
            f"{' (TRUNCATED)' if truncated else ''}, top aliases="
            f"{[f'{a}x{n}' for a, n in top_aliases]}"
        )
        return ctx

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
6. **CRITICAL: Always use FULL disease names, never abbreviations alone**. For example:
   - Write "Dilated Cardiomyopathy (DCM)" NOT just "DCM"
   - Write "Hypertrophic Cardiomyopathy (HCM)" NOT just "HCM"
   - Write "Tibial Muscular Dystrophy (TMD)" NOT just "TMD"
   - Write "Limb-Girdle Muscular Dystrophy (LGMD)" NOT just "LGMD"

### COMMON TTN-RELATED DISEASES (use these full names)
- **Cardiac diseases:**
  - Dilated Cardiomyopathy (DCM)
  - Hypertrophic Cardiomyopathy (HCM)
  - Restrictive Cardiomyopathy (RCM)
  - Left Ventricular Noncompaction (LVNC)
- **Skeletal muscle diseases:**
  - Tibial Muscular Dystrophy (TMD) — NOT "TTN Muscular Dystrophy"
  - Limb-Girdle Muscular Dystrophy (LGMD)
  - Hereditary Myopathy with Early Respiratory Failure (HMERF)
  - Centronuclear Myopathy
- **Both**
  - Early-onset myopathy with fatal cardiomyopathy

## OUTPUT FORMAT
Respond with ONLY a JSON object (no extra text):
```json
{{
  "reasoning": "<1-2 sentences on what you found>",
  "disease": "<FULL disease name with abbreviation in parentheses, e.g. 'Dilated Cardiomyopathy (DCM)', or 'Not specified'>",
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

    # ── B1: Self-consistency majority voting ─────────────────────────────────

    @staticmethod
    def _vote_field(values: List, default: str = "Not specified"):
        """
        Pick the most-common informative value across n samples.
        "Not specified" / empty are only returned if every sample agrees on them.
        """
        cleaned = [
            v for v in values
            if v is not None and str(v).strip() not in ("", "Not specified")
        ]
        if not cleaned:
            return default
        counter = Counter(str(v).strip() for v in cleaned)
        return counter.most_common(1)[0][0]

    def _majority_vote(
        self,
        samples: List[Dict],
        raw_text: Optional[str] = None,
        article: Optional[Dict] = None,
    ) -> Dict:
        """
        Majority-vote across n parsed samples for each field independently.

        For evidence_sentence we pick (in order of preference):
          1. A sample whose evidence is grounded in the source text AND
             whose disease matches the voted disease — longest first.
          2. Any sample whose disease matches the voted disease — longest first.
          3. Any sample — longest evidence first.

        This prevents the case where a longer-but-ungrounded quote ("...|...|...")
        beats a shorter-but-grounded quote, then gets nuked by B4 evidence
        verification downstream.
        """
        if not samples:
            return self._get_empty_result("No samples to vote on")
        if len(samples) == 1:
            return samples[0]

        voted_disease = self._vote_field([s.get("disease") for s in samples])
        voted_tissue = self._vote_field([s.get("tissue_affected") for s in samples])
        voted_age = self._vote_field([s.get("age_onset") for s in samples])
        voted_inh = self._vote_field([s.get("inheritance") for s in samples])

        # patient_count: median of integer-valued samples
        counts: List[int] = []
        for s in samples:
            try:
                counts.append(int(s.get("patient_count", 0) or 0))
            except (ValueError, TypeError):
                continue
        if counts:
            counts.sort()
            voted_count = counts[len(counts) // 2]
        else:
            voted_count = 0

        # Candidate samples whose disease matches the vote (preferred)
        disease_winners = [
            s for s in samples
            if str(s.get("disease", "")).strip() == voted_disease
        ]
        candidates = disease_winners or list(samples)

        # If we can verify grounding, prefer grounded evidence
        chosen = None
        if (
            LLM_VERIFY_EVIDENCE
            and raw_text is not None
            and article is not None
        ):
            grounded = []
            for s in candidates:
                ev = str(s.get("evidence_sentence", "") or "").strip()
                if not ev:
                    continue
                ok, score = self._verify_evidence_grounding(ev, raw_text, article)
                if ok:
                    grounded.append((score, len(ev), s))
            if grounded:
                # higher score first, then longer evidence
                grounded.sort(key=lambda t: (t[0], t[1]), reverse=True)
                chosen = grounded[0][2]

        if chosen is None:
            # Fall back to longest evidence among disease-matching samples
            ranked = sorted(
                candidates,
                key=lambda s: len(str(s.get("evidence_sentence", "") or "")),
                reverse=True,
            )
            chosen = ranked[0]

        return {
            "reasoning": str(chosen.get("reasoning", "")).strip(),
            "disease": voted_disease,
            "tissue_affected": voted_tissue,
            "age_onset": voted_age,
            "inheritance": voted_inh,
            "patient_count": voted_count,
            "evidence_sentence": str(chosen.get("evidence_sentence", "")).strip(),
        }

    # ── B4: Evidence grounding verification ──────────────────────────────────

    @staticmethod
    def _normalize_for_match(s: str) -> str:
        """
        Loose normalisation that survives docling/OCR table noise.

        - Lower-case
        - Treat all non-alphanumeric chars (incl. `|`, `>`, `.`, `:`) as a single
          space, so `c.77989C>T`, `c.77989C . T`, and `c.77989C | T` all collapse
          to the same `c 77989c t` token sequence.
        - Collapse runs of whitespace.
        """
        s = s.lower()
        s = re.sub(r"[^a-z0-9]+", " ", s)
        s = re.sub(r"\s+", " ", s)
        return s.strip()

    def _verify_evidence_grounding(
        self, evidence: str, raw_text: str, article: Dict
    ) -> Tuple[bool, float]:
        """
        Check that `evidence` is supported by the source text.
        Returns (is_grounded, overlap_score).

        Why this is tricky: LLMs often quote evidence by stitching together
        non-adjacent table cells, e.g. "c.77989C . T p.Q23354X" — the original
        text has "c 77989c t 326 c 70060c t 275 p q23354x" with extra tokens
        between the quoted parts. A naive 5-gram check returns 0 even though
        every quoted token appears in the source.

        Strategy:
          1. Normalize whitespace + casing for both.
          2. Direct substring match → grounded with score 1.0.
          3. Multi-scale n-gram overlap (n=2..5). Take the max overlap ratio
             across scales. Variant identifiers like `77989c` are typically
             rare 1-2-token anchors, so requiring rare tokens to appear helps.
          4. "Anchor coverage": fraction of *informative* tokens (length>=4 or
             contains a digit) that appear in haystack. If >= threshold, grounded.
          5. Drop empties / "Not specified" / very short (<3 tokens) evidence.
        """
        if not evidence:
            return False, 0.0
        ev = self._normalize_for_match(evidence)
        if not ev or ev in ("not specified", "n/a", "none"):
            return False, 0.0

        haystack = " ".join(
            filter(None, [
                article.get("title", ""),
                article.get("abstract", ""),
                raw_text or "",
            ])
        )
        haystack_norm = self._normalize_for_match(haystack)
        if not haystack_norm:
            return False, 0.0

        # Step 1: direct substring match
        if ev in haystack_norm:
            return True, 1.0

        ev_tokens = ev.split()
        if len(ev_tokens) < 3:
            return False, 0.0

        # Step 2: multi-scale n-gram overlap, take best ratio across n=2..5
        best_ratio = 0.0
        for n in range(2, 6):
            if len(ev_tokens) < n:
                break
            ngrams = {
                " ".join(ev_tokens[i : i + n])
                for i in range(len(ev_tokens) - n + 1)
            }
            if not ngrams:
                continue
            hits = sum(1 for ng in ngrams if ng in haystack_norm)
            ratio = hits / len(ngrams)
            best_ratio = max(best_ratio, ratio)

        # Step 3: informative-token coverage (catches stitched-cell quotes)
        # An informative token: length >= 4 OR contains a digit (e.g. "77989c",
        # "q23354x", "1044", "p", "c"). Common stop-tokens like "the" / "and"
        # add no signal so we ignore them.
        informative = [
            t for t in ev_tokens
            if (len(t) >= 4 or any(ch.isdigit() for ch in t))
        ]
        anchor_ratio = 0.0
        if informative:
            anchor_hits = sum(1 for t in informative if f" {t} " in f" {haystack_norm} ")
            anchor_ratio = anchor_hits / len(informative)

        # Final score = best of the two signals
        score = max(best_ratio, anchor_ratio)
        return score >= LLM_EVIDENCE_NGRAM_OVERLAP_THRESHOLD, score

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