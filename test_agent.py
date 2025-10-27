#!/usr/bin/env python3
"""
Test script for TTN Variant AI Agent
Tests individual components and the complete pipeline
"""

import sys
from pathlib import Path

# Test variants
TEST_VARIANTS = [
    "2-178527121-A-G",  # Pathogenic
    "2-178527628-C-T",  # Benign
    "2-178612477-T-A",  # Example
]


def test_imports():
    """Test if all required modules can be imported"""
    print("Testing imports...")
    try:
        from utils.variant_parser import parse_variant
        from utils.evo2_predictor import Evo2Predictor
        from utils.pubmed_search import PubMedSearcher
        from utils.image_generator import ImageGenerator
        from utils.html_report import HTMLReportGenerator
        print("All imports successful")
        return True
    except ImportError as e:
        print(f"Import error: {e}")
        return False


def test_variant_parser():
    """Test variant parser"""
    print("\nTesting variant parser...")
    from utils.variant_parser import parse_variant, validate_ttn_variant
    
    try:
        # Test valid variants
        for variant in TEST_VARIANTS:
            parsed = parse_variant(variant)
            is_ttn = validate_ttn_variant(parsed)
            print(f"{variant}: {parsed['variant_id']} (TTN: {is_ttn})")
        
        # Test invalid variant
        try:
            parse_variant("invalid-variant")
            print("Should have raised error for invalid variant")
            return False
        except Exception:
            pass  # Expected
        
        print("Variant parser working correctly")
        return True
    except Exception as e:
        print(f"Variant parser error: {e}")
        return False


def test_evo2_predictor():
    """Test Evo2 predictor"""
    print("\nTesting Evo2 predictor...")
    from utils.variant_parser import parse_variant
    from utils.evo2_predictor import Evo2Predictor
    
    try:
        predictor = Evo2Predictor()
        variant_info = parse_variant(TEST_VARIANTS[0])
        
        print(f"  Testing with {variant_info['variant_id']}...")
        result = predictor.predict(variant_info)
        
        if result.get('success'):
            print(f"  Prediction: {result['prediction']}")
            print(f"  Delta score: {result['delta_score']:.6f}")
            print("Evo2 predictor working")
            return True
        else:
            print(f"Evo2 prediction failed: {result.get('error')}")
            return False
    except Exception as e:
        print(f"Evo2 predictor error: {e}")
        return False


def test_pubmed_search():
    """Test PubMed search"""
    print("\nTesting PubMed search...")
    from utils.variant_parser import parse_variant
    from utils.pubmed_search import PubMedSearcher
    
    try:
        searcher = PubMedSearcher()
        variant_info = parse_variant(TEST_VARIANTS[0])
        
        print(f"Searching for {variant_info['variant_id']}...")
        results = searcher.search(variant_info)
        
        print(f"Found {len(results)} articles")
        if results:
            print(f"  Example: {results[0]['title'][:60]}...")
        
        print("PubMed search working")
        return True
    except Exception as e:
        print(f"PubMed search error: {e}")
        return False


def test_image_generator():
    """Test image generator"""
    print("\nTesting image generator...")
    from utils.variant_parser import parse_variant
    from utils.image_generator import ImageGenerator
    
    try:
        generator = ImageGenerator()
        variant_info = parse_variant(TEST_VARIANTS[0])
        
        print(f"  Generating schematic for {variant_info['variant_id']}...")
        image_path = generator.generate_titin_schematic(variant_info)
        
        if image_path.exists():
            print(f"Image saved to: {image_path}")
            print("Image generator working")
            return True
        else:
            print("Image not generated")
            return False
    except Exception as e:
        print(f"Image generator error: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_complete_pipeline():
    """Test complete pipeline"""
    print("\nTesting complete pipeline...")
    
    try:
        import subprocess
        result = subprocess.run(
            [sys.executable, "main.py", TEST_VARIANTS[0], "--skip-evo2"],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0:
            print("Complete pipeline working")
            return True
        else:
            print(f"Pipeline failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"Pipeline error: {e}")
        return False


def main():
    """Run all tests"""
    print("=" * 70)
    print("TTN Variant AI Agent - Component Tests")
    print("=" * 70)
    
    tests = [
        ("Imports", test_imports),
        ("Variant Parser", test_variant_parser),
        ("Evo2 Predictor", test_evo2_predictor),
        ("PubMed Search", test_pubmed_search),
        ("Image Generator", test_image_generator),
        ("Complete Pipeline", test_complete_pipeline),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"{name} crashed: {e}")
            results.append((name, False))
    
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    for name, success in results:
        status = "PASS" if success else "FAIL"