
import sys
import os
sys.path.append('/home/ryan910702/ttn_agent')
from utils.local_clinical_extractor import LocalClinicalExtractor

def test_prompt_creation():
    # Mock backend to avoid loading model
    class MockExtractor(LocalClinicalExtractor):
        def __init__(self):
            self.backend = "mock"
            self.model_name = "mock"
            self.tensor_parallel_size = 1
            
    extractor = MockExtractor()
    
    text = "This is a sample text about variant c.123G>A."
    variant_id = "c.123G>A"
    aliases = ["c.123G>A", "p.Trp41*"]
    
    prompt = extractor._create_few_shot_prompt(text, variant_id, aliases)
    print("--- Generated Prompt ---")
    print(prompt)
    print("------------------------")

if __name__ == "__main__":
    test_prompt_creation()
