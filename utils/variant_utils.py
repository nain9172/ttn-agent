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


def _strip_hgvs_prefix(hgvs: str) -> Optional[str]:
    """
    剝掉 HGVS accession 前綴，只保留 p./c./n./g./m. 開頭的那段。

    例如:
        NP_597681.4:p.Trp19400Ter   -> p.Trp19400Ter
        NM_133378.4:c.77115G>A      -> c.77115G>A
        ENST00000589042.6:c.123A>G  -> c.123A>G
        LRG_391t1:c.456C>T          -> c.456C>T
        p.Trp19400Ter               -> None  (沒有前綴可剝)
    """
    if not hgvs:
        return None
    # 抓出第一個 p./c./n./g./m./r. 開頭的子字串
    m = re.search(r'\b([pcngmr]\.[^\s,;]+)', hgvs)
    if not m:
        return None
    stripped = m.group(1)
    return stripped if stripped != hgvs.strip() else None


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
    aliases: List[str] = [variant_id]  # e.g., 2-178612477-T-A

    def _add_hgvs(hgvs_value):
        """把單一 HGVS 字串展開：原字串 + 剝前綴版本 + 單字母蛋白質版本。"""
        if not hgvs_value:
            return
        h = str(hgvs_value).strip()
        if not h:
            return
        aliases.append(h)
        # 剝掉 NP_xxx:/NM_xxx:/ENST xxx: 前綴，只留 p./c./n./g./m. 部分
        stripped = _strip_hgvs_prefix(h)
        if stripped:
            aliases.append(stripped)
        # 三字母 → 單字母蛋白質格式
        single_letter = _convert_protein_hgvs_to_single_letter(h)
        if single_letter:
            aliases.append(single_letter)

    if clinvar_info:
        for key in ('hgvs', 'hgvs_nucleotide', 'hgvs_protein'):
            if key not in clinvar_info:
                continue
            value = clinvar_info[key]
            if isinstance(value, list):
                for item in value:
                    _add_hgvs(item)
            else:
                _add_hgvs(value)

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
