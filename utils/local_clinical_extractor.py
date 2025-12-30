#!/usr/bin/env python3
"""
Improved Local Clinical Information Extractor
Fixes:
- Moves 'Reasoning' to top (Chain-of-Thought) for better accuracy
- Removes brittle JSON forcing
- Handles specific variant isolation strictly (IMPROVED)
"""

import logging
import json
import re
import os
from datetime import datetime
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

class LocalClinicalExtractor:
    def __init__(self, backend: str = "vllm", model_name: str = "meta-llama/Llama-3.2-3B", tensor_parallel_size: int = 2):
        self.backend = backend
        self.model_name = model_name
        self.tensor_parallel_size = tensor_parallel_size
        self.is_small_model = "3B" in model_name.upper() or "3b" in model_name.lower()
        self.max_context_length = 12000 if not self.is_small_model else 6000
        
        logger.info(f"Initializing {backend} backend with model: {model_name}")
        
        if backend == "vllm":
            self._init_vllm()
        elif backend == "ollama":
            self._init_ollama()
        else:
            raise ValueError(f"Unsupported backend: {backend}")
    
    def _init_vllm(self):
        try:
            from vllm import LLM, SamplingParams
            self.llm = LLM(
                model=self.model_name,
                tensor_parallel_size=self.tensor_parallel_size,
                trust_remote_code=True,
                gpu_memory_utilization=0.9,
                max_model_len=8192,
                enforce_eager=True
            )
            # Slightly higher temp to allow reasoning flow, but still focused
            # 增加 max_tokens 避免回覆被截斷，並添加 stop 條件避免重複
            self.sampling_params = SamplingParams(
                temperature=0.7, 
                top_p=0.95, 
                max_tokens=2048,
                stop=["\n\n\n", "```\n\n"]  # 添加 stop 條件避免過度生成
            )
        except Exception as e:
            logger.error(f"vLLM init failed: {e}")
            raise

    def _init_ollama(self):
        import ollama
        self.ollama_client = ollama

    def _log_prompt_response(self, pmid: str, prompt: str, response: str, variant_id: str = None):
        """Log prompt and LLaMA response to a file for debugging/analysis."""
        log_base_dir = os.path.join(os.path.dirname(__file__), "..", "llm_logs")
        
        # 如果有 variant_id，創建以變異命名的子資料夾
        if variant_id:
            log_dir = os.path.join(log_base_dir, variant_id)
        else:
            log_dir = log_base_dir
        
        os.makedirs(log_dir, exist_ok=True)
        
        # 文件名只用 PMID
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
            f.write("\n")
            f.write("=" * 80 + "\n")
        
        logger.info(f"Logged prompt/response to: {log_file}")

    def batch_extract(self, articles: List[Dict], variant_info: Dict, clinvar_info: Optional[Dict] = None) -> List[Dict]:
        variant_id = f"{variant_info['chrom']}-{variant_info['pos']}-{variant_info['ref']}-{variant_info['alt']}"
        logger.info(f"Starting extraction for {len(articles)} articles. Target: {variant_id}")
        
        enriched = []
        for article in articles:
            # Skip empty articles
            if not any([article.get('title'), article.get('abstract'), article.get('full_text')]):
                continue

            info = self._extract_single(article, variant_id, clinvar_info)
            article['clinical_info'] = info
            
            # Log result for debugging
            if info.get('disease') != "Not specified":
                logger.info(f"  [+] PMID {article.get('pmid')}: Found {info.get('disease')} ({info.get('tissue_affected')})")
            else:
                logger.info(f"  [-] PMID {article.get('pmid')}: No relevant info")
                
            enriched.append(article)
        return enriched

    def _extract_single(self, article: Dict, variant_id: str, clinvar_info: Optional[Dict]) -> Dict:
        text = article.get('text_for_llm') or article.get('full_text') or article.get('abstract') or ""
        aliases = self._get_aliases(variant_id, clinvar_info)
        relevant_text = self._prepare_context(text, aliases)
        
        if len(relevant_text) < 50:
            return self._get_empty_result("Text too short")

        prompt = self._build_cot_prompt(relevant_text, aliases)
        pmid = article.get('pmid', 'unknown')
        
        try:
            if self.backend == "vllm":
                raw = self._generate_vllm(prompt)
            else:
                raw = self._generate_ollama(prompt)
            
            # Log prompt and response to file
            self._log_prompt_response(pmid, prompt, raw, variant_id)
            
            data = self._robust_json_parse(raw)
            data['raw_llm_output'] = raw
            return data
        except Exception as e:
            logger.error(f"Extraction error: {e}")
            return self._get_empty_result(f"Error: {str(e)}")

    def _get_aliases(self, variant_id: str, clinvar_info: Optional[Dict]) -> List[str]:
        """
        Generates a comprehensive list of known variant notations/aliases.
        This greatly aids the LLM in identifying the target variant in the text.
        """
        aliases = [variant_id] # e.g., 2-178612477-T-A
        
        if clinvar_info:
            # Use any HGVS found during ClinVar search/scrape
            # Combined list of all HGVS
            if 'hgvs' in clinvar_info:
                h = clinvar_info['hgvs']
                if isinstance(h, list): aliases.extend(h)
                else: aliases.append(h)
            
            # Nucleotide HGVS (e.g., NM_001267550.2:c.3208G>A)
            if 'hgvs_nucleotide' in clinvar_info:
                h = clinvar_info['hgvs_nucleotide']
                if isinstance(h, list): aliases.extend(h)
                else: aliases.append(h)
            
            # Protein HGVS (e.g., NP_001254479.2:p.Glu1070Lys)
            if 'hgvs_protein' in clinvar_info:
                h = clinvar_info['hgvs_protein']
                if isinstance(h, list): aliases.extend(h)
                else: aliases.append(h)
            
            # rsID (e.g., rs1057518195)
            if 'rsid' in clinvar_info:
                aliases.append(clinvar_info['rsid'])
            
            # Include VCV ID if available (VCVXXXXXXXXX)
            # if 'vcv' in clinvar_info:
            #     aliases.append(clinvar_info['vcv'])
        
        parts = variant_id.split('-')
        if len(parts) == 4:
            c, p, r, a = parts
            p_str = str(p)
            
            # Alias 1: Minimal VCF-like notation
            aliases.append(f"{p_str}{r}>{a}") # e.g., 178612477T>A
            
            # Alias 2: Chromosome:Position notation (with and without chr prefix)
            aliases.append(f"chr{c}:{p_str}")
            aliases.append(f"{c}:{p_str}")
            
            # Alias 3: Canonical HGVS genomic notation (GRCh38, positive strand)
            # Format: NC_000002.12:g.{pos}{ref}>{alt}
            # hgvs_g = f"NC_000002.12:g.{p_str}{r}>{a}"
            # aliases.append(hgvs_g)

            # Alias 4: Common text format (e.g., T178612477A)
            # aliases.append(f"{r}{p_str}{a}")

        # Filter out Nones, empty strings, and duplicates
        return list(set([str(a).strip() for a in aliases if a]))

    def _prepare_context(self, text: str, aliases: List[str]) -> str:
        if len(text) <= self.max_context_length: return text
        
        text_lower = text.lower()
        best_pos = -1
        for alias in aliases:
            pos = text_lower.find(alias.lower())
            if pos != -1:
                best_pos = pos
                break
        
        if best_pos != -1:
            start = max(0, best_pos - (self.max_context_length // 2))
            end = min(len(text), start + self.max_context_length)
            return "... " + text[start:end] + " ..."
        return text[:self.max_context_length]

    def _build_cot_prompt(self, text: str, aliases: List[str]) -> str:
        # Create a clearly formatted list of all aliases
        target_variant_list = "\n- " + "\n- ".join(aliases)
        # print(f"=====================target_variant_list: {target_variant_list}=====================")
        # Enhanced Prompt with STRICT filtering instructions
        return f"""<|begin_of_text|><|start_header_id|>system<|end_header_id|>
You are a precision medical extraction AI.
Your task is to extract clinical information specifically and only for the target variant mentioned below.

TARGET VARIANT ALIASES (Must be one of these in the text):
{target_variant_list}

STRICT INSTRUCTIONS:
1. **IDENTIFICATION:** Scan the entire 'Article Text' ONLY for clear and explicit mentions of the TARGET VARIANT ALIASES.
2. **FILTERING:** IGNORE ALL information related to *other variants, other genes, or general disease descriptions* not explicitly linked to the TARGET VARIANT.
3. **EVIDENCE:** The 'evidence_sentence' MUST be a **direct quote** from the text that explicitly mentions one of the TARGET VARIANT ALIASES and a corresponding clinical finding (e.g., phenotype, tissue affected). If no such sentence is found, stop.
4. **FALLBACK:** If no patient data is found that explicitly mentions the target variant, set all fields (disease, tissue_affected, etc.) to "Not specified" and 'patient_count' to 0.
5. You must output a valid JSON object.

JSON STRUCTURE (Reasoning First):
{{
  "reasoning": "Step-by-step thinking: Did the text explicitly mention any of the TARGET VARIANT ALIASES? If yes, what specific patient details were associated with that variant?",
  "evidence_sentence": "Direct quote from text verifying the finding",
  "disease": "Specific diagnosis of input variant (e.g. Dilated Cardiomyopathy, muscular dystrophy)",
  "tissue_affected": ONLY "Cardiac", "Skeletal", "Both", or "Not specified",
  "age_onset": "Congenital or Adult",
  "inheritance": "Pattern (e.g. Autosomal Dominant, Sporadic)",
  "patient_count": Integer (0 if none found)
}}
Notice the tissue_affected can only be "Cardiac", "Skeletal", "Both", or "Not specified".
<|eot_id|><|start_header_id|>user<|end_header_id|>
Article Text:
{text}
<|eot_id|><|start_header_id|>assistant<|end_header_id|>
```json"""

    def _generate_vllm(self, prompt: str) -> str:
        outputs = self.llm.generate([prompt], self.sampling_params)
        return outputs[0].outputs[0].text

    def _generate_ollama(self, prompt: str) -> str:
        return self.ollama_client.generate(model=self.model_name, prompt=prompt)['response']

    def _robust_json_parse(self, text: str) -> Dict:
        # If the model didn't output ```json, it might just be the object
        # or it might be ```json ... ```
        clean_text = text.strip()
        
        # Remove markdown code blocks if present
        if clean_text.startswith("```json"):
            clean_text = clean_text[7:]
        if clean_text.startswith("```"):
            clean_text = clean_text[3:]
        if clean_text.endswith("```"):
            clean_text = clean_text[:-3]
        
        clean_text = clean_text.strip()
        
        # Ensure it starts with {
        start = clean_text.find('{')
        end = clean_text.rfind('}')
        
        if start != -1 and end != -1:
            json_str = clean_text[start:end+1]
            try:
                return json.loads(json_str)
            except json.JSONDecodeError as e:
                # Try cleaning trailing commas which LLMs often output
                try:
                    json_str_clean = re.sub(r',\s*}', '}', json_str)
                    json_str_clean = re.sub(r',\s*]', ']', json_str_clean)
                    return json.loads(json_str_clean)
                except json.JSONDecodeError as e2:
                    logger.warning(f"JSON parse error: {e2}")
                    logger.debug(f"Raw LLM output (first 500 chars): {text[:500]}")

        # Log when we can't find JSON structure at all
        if start == -1 or end == -1:
            logger.warning(f"No JSON found in LLM output. Raw output (first 300 chars): {text[:300]}")
        
        return self._get_empty_result("JSON Parsing Failed")

    def _get_empty_result(self, reason: str) -> Dict:
        return {
            "reasoning": reason,
            "disease": "Not specified",
            "tissue_affected": "Not specified",
            "age_onset": "Not specified",
            "inheritance": "Not specified",
            "evidence_sentence": "N/A",
            "patient_count": 0
        }

    def aggregate_stats(self, articles: List[Dict]) -> Dict:
        stats = {
            'total_articles': len(articles),
            'articles_with_extraction': 0,
            'total_patients': 0,
            'total_families': 0,
            'disease_distribution': {},
            'inheritance_distribution': {},
            'age_onset_distribution': {},
            'tissue_affected': {'cardiac': 0, 'skeletal_muscle': 0, 'both': 0},
            'cardiac_phenotypes': {},
            'skeletal_phenotypes': {}
        }
        for a in articles:
            info = a.get('clinical_info', {})
            # Only count if extraction was successful and meaningful
            if info.get('disease') not in ["Not specified", None, "JSON Parsing Failed"]:
                stats['articles_with_extraction'] += 1
                try: stats['total_patients'] += int(info.get('patient_count', 0))
                except: pass
                
                d = info.get('disease', 'Unknown')
                stats['disease_distribution'][d] = stats['disease_distribution'].get(d, 0) + 1
                
                i = info.get('inheritance', 'Unknown')
                stats['inheritance_distribution'][i] = stats['inheritance_distribution'].get(i, 0) + 1
                
                t = info.get('tissue_affected', 'Not specified').lower()
                if 'cardiac' in t: stats['tissue_affected']['cardiac'] += 1
                if 'skeletal' in t: stats['tissue_affected']['skeletal_muscle'] += 1
                if 'both' in t: stats['tissue_affected']['both'] += 1
                
        return stats