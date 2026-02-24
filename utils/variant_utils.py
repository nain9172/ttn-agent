#!/usr/bin/env python3
"""
Variant Utilities
共享的變異處理工具函數
"""

import re
from typing import Dict, List, Optional

# 氨基酸三字母到單字母的轉換表
AA_THREE_TO_ONE = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*', 'Stop': '*'
}


def _convert_protein_hgvs_to_single_letter(hgvs: str) -> Optional[str]:
    """
    將三字母氨基酸 HGVS 轉換為單字母版本
    例如: p.Ile35947Asn -> I35947N
    """
    # 匹配蛋白質變異格式: p.Xxx123Yyy 或 NP_xxx:p.Xxx123Yyy
    pattern = r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Stop|Ter)'
    match = re.search(pattern, hgvs)
    
    if match:
        ref_aa_three = match.group(1)
        position = match.group(2)
        alt_aa_three = match.group(3)
        
        # 轉換為單字母
        ref_aa_one = AA_THREE_TO_ONE.get(ref_aa_three)
        alt_aa_one = AA_THREE_TO_ONE.get(alt_aa_three, alt_aa_three)
        
        if ref_aa_one and alt_aa_one:
            return f"{ref_aa_one}{position}{alt_aa_one}"
    
    return None


def get_variant_aliases(variant_id: str, clinvar_info: Optional[Dict]) -> List[str]:
    """
    從 ClinVar 信息中提取變異的所有別名
    這些信息已經在 ClinVar 頁面抓取時獲得，不需要重新生成
    
    特別注意：會自動生成單字母和三字母氨基酸格式的變異
    
    Args:
        variant_id: 變異 ID (格式: chrom-pos-ref-alt)
        clinvar_info: ClinVar 信息字典（包含 HGVS、rsID 等）
    
    Returns:
        變異別名列表
    """
    aliases = [variant_id]  # e.g., 2-178612477-T-A
    
    if clinvar_info:
        # Use any HGVS found during ClinVar search/scrape
        # Combined list of all HGVS
        if 'hgvs' in clinvar_info:
            h = clinvar_info['hgvs']
            if isinstance(h, list):
                aliases.extend(h)
                # 生成單字母版本
                for hgvs in h:
                    single_letter = _convert_protein_hgvs_to_single_letter(str(hgvs))
                    if single_letter:
                        aliases.append(single_letter)
            else:
                aliases.append(h)
                single_letter = _convert_protein_hgvs_to_single_letter(str(h))
                if single_letter:
                    aliases.append(single_letter)
        
        # Nucleotide HGVS (e.g., NM_001267550.2:c.3208G>A)
        if 'hgvs_nucleotide' in clinvar_info:
            h = clinvar_info['hgvs_nucleotide']
            if isinstance(h, list):
                aliases.extend(h)
            else:
                aliases.append(h)
        
        # Protein HGVS (e.g., NP_001254479.2:p.Glu1070Lys)
        if 'hgvs_protein' in clinvar_info:
            h = clinvar_info['hgvs_protein']
            if isinstance(h, list):
                aliases.extend(h)
                # 生成單字母版本
                for hgvs in h:
                    single_letter = _convert_protein_hgvs_to_single_letter(str(hgvs))
                    if single_letter:
                        aliases.append(single_letter)
            else:
                aliases.append(h)
                single_letter = _convert_protein_hgvs_to_single_letter(str(h))
                if single_letter:
                    aliases.append(single_letter)
        
        # rsID (e.g., rs1057518195)
        if 'rsid' in clinvar_info:
            aliases.append(clinvar_info['rsid'])
    
    # 從 variant_id 生成簡化格式
    parts = variant_id.split('-')
    if len(parts) == 4:
        c, p, r, a = parts
        p_str = str(p)
        
        # Alias 1: Minimal VCF-like notation
        aliases.append(f"{p_str}{r}>{a}")  # e.g., 178612477T>A
        
        # Alias 2: Chromosome:Position notation (with and without chr prefix)
        aliases.append(f"chr{c}:{p_str}")
        aliases.append(f"{c}:{p_str}")
    
    # Filter out Nones, empty strings, and duplicates
    return list(set([str(a).strip() for a in aliases if a]))
