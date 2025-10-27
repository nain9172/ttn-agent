"""
Variant Parser Module
Parses variant notation and validates input
"""

import re
from typing import Dict


class VariantParseError(Exception):
    """Custom exception for variant parsing errors"""
    pass


def parse_variant(variant_string: str) -> Dict[str, str]:
    """
    Parse variant string in format: chromosome-position-ref-alt
    
    Args:
        variant_string: Variant notation (e.g., "2-178612477-T-A")
    
    Returns:
        Dictionary with parsed variant information
    
    Raises:
        VariantParseError: If variant format is invalid
    """
    variant_string = variant_string.strip()
    
    # Try different separators
    for sep in ['-', '_', ':']:
        if variant_string.count(sep) == 3:
            parts = variant_string.split(sep)
            break
    else:
        raise VariantParseError(
            f"Invalid variant format: {variant_string}. "
            "Expected format: chromosome-position-ref-alt (e.g., 2-178612477-T-A)"
        )
    
    if len(parts) != 4:
        raise VariantParseError(
            f"Invalid variant format: {variant_string}. "
            "Expected 4 parts: chromosome, position, ref, alt"
        )
    
    chrom, pos, ref, alt = parts
    
    # Validate chromosome
    chrom = chrom.replace('chr', '').upper()
    valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M', 'MT']
    if chrom not in valid_chroms:
        raise VariantParseError(f"Invalid chromosome: {chrom}")
    
    # Validate position
    try:
        pos = int(pos)
        if pos <= 0:
            raise ValueError
    except ValueError:
        raise VariantParseError(f"Invalid position: {parts[1]}")
    
    # Validate ref and alt bases
    valid_bases = set('ACGT')
    ref = ref.upper()
    alt = alt.upper()
    
    if not all(b in valid_bases for b in ref):
        raise VariantParseError(f"Invalid reference base: {ref}")
    if not all(b in valid_bases for b in alt):
        raise VariantParseError(f"Invalid alternate base: {alt}")
    
    # For SNVs, ref and alt should be single bases
    if len(ref) != 1 or len(alt) != 1:
        raise VariantParseError(
            "Currently only single nucleotide variants (SNVs) are supported"
        )
    
    return {
        'chrom': chrom,
        'pos': pos,
        'ref': ref,
        'alt': alt,
        'variant_id': f"{chrom}-{pos}-{ref}-{alt}"
    }


def format_variant_hgvs(variant_info: Dict[str, str]) -> str:
    """
    Format variant in HGVS notation
    
    Args:
        variant_info: Parsed variant dictionary
    
    Returns:
        HGVS formatted string
    """
    return (
        f"NC_000002.12:g.{variant_info['pos']}"
        f"{variant_info['ref']}>{variant_info['alt']}"
    )


def validate_ttn_variant(variant_info: Dict[str, str]) -> bool:
    """
    Check if variant is within TTN gene region
    
    Args:
        variant_info: Parsed variant dictionary
    
    Returns:
        True if variant is in TTN region, False otherwise
    """
    from config import TTN_GENE_INFO
    
    if variant_info['chrom'] != TTN_GENE_INFO['chromosome']:
        return False
    
    pos = variant_info['pos']
    return TTN_GENE_INFO['start'] <= pos <= TTN_GENE_INFO['end']