
import logging
import sys
from utils.local_clinical_extractor import LocalClinicalExtractor
from unittest.mock import MagicMock

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_context_extraction():
    # Mock init to avoid loading model
    original_init = LocalClinicalExtractor.__init__
    LocalClinicalExtractor.__init__ = lambda self, **kwargs: None
    
    extractor = LocalClinicalExtractor()
    # Restore init just in case
    LocalClinicalExtractor.__init__ = original_init
    
    # Mock logger since we bypassed init
    extractor.backend = "vllm"
    extractor.model_name = "test"
    
    # Mock data
    variant_id = "2-178528273-G-C"
    aliases = [variant_id, "c.123G>C", "p.Arg41Gln"]
    
    # Create a long text with two variants
    text = """
    Introduction
    Dilated cardiomyopathy is a serious condition. We studied several variants.
    
    Results
    Patient A carried the variant c.456T>A (Variant B). This patient showed no symptoms and had normal heart function.
    
    On the other hand, Patient B carried the variant c.123G>C (Target Variant). This patient presented with severe dilated cardiomyopathy at age 45.
    The phenotype was characterized by left ventricular enlargement.
    
    Discussion
    The variant c.456T>A appears to be benign. However, c.123G>C is likely pathogenic.
    """
    
    # Pad text to exceed 3000 chars to trigger smart context
    padding = "bla " * 600
    full_text = text + "\n" + padding
    
    print(f"Full text length: {len(full_text)}")
    
    # Test extraction (Paragraph based)
    relevant_context = extractor._get_relevant_context(full_text, aliases)
    
    print("\n--- Relevant Context ---")
    print(relevant_context)
    print("------------------------")
    
    # Verification
    if "c.123G>C" in relevant_context:
        print("SUCCESS: Target variant found in context.")
    else:
        print("FAILURE: Target variant NOT found in context.")
        
    if "c.456T>A" in relevant_context:
        print("NOTE: Distractor variant also found (might be in same paragraph).")
        
    # Test prompt generation (Few-Shot)
    prompt = extractor._create_few_shot_prompt(relevant_context, variant_id, aliases)
    print("\n--- Generated Prompt ---")
    print(prompt[:500] + "...")
    print("------------------------")
    
    if "Example 1 (Positive Case)" in prompt:
        print("SUCCESS: Prompt contains few-shot examples.")
    else:
        print("FAILURE: Prompt missing few-shot examples.")
        
    if "IGNORE total cohort sizes" in prompt:
        print("SUCCESS: Prompt contains cohort ignore instruction.")
    else:
        print("FAILURE: Prompt missing cohort ignore instruction.")

    # Test Case 2: Tricky Cohort vs Variant
    print("\n--- Test Case 2: Cohort vs Variant ---")
    tricky_text = """
    Methods: We recruited a total cohort of 500 patients (mean age 60 +/- 10 years).
    
    Results:
    Genetic testing revealed that only 2 patients carried the variant c.123G>C (Target Variant).
    These two patients were much younger (ages 25 and 28).
    """
    
    # Verify prompt generation for this case
    prompt2 = extractor._create_few_shot_prompt(tricky_text, variant_id, aliases)
    if "Example 2 (Tricky Case" in prompt2:
        print("SUCCESS: Prompt contains tricky example.")
    else:
        print("FAILURE: Prompt missing tricky example.")

if __name__ == "__main__":
    try:
        test_context_extraction()
    except Exception as e:
        print(f"Error: {e}")
