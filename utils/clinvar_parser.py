"""
ClinVar Parser Module
使用 NCBI E-utilities API 查詢 ClinVar 變異資訊（可選功能）
支援網頁抓取作為備用方法
"""

import logging
import time
import re
from typing import Dict, List, Optional, Set

try:
    from Bio import Entrez
except ImportError:
    print("Warning: Biopython not installed. ClinVar API access will not work.")
    Entrez = None

try:
    import requests
    from bs4 import BeautifulSoup
    SCRAPING_AVAILABLE = True
except ImportError:
    print("Warning: requests or beautifulsoup4 not installed. Web scraping will not work.")
    SCRAPING_AVAILABLE = False

from config import PUBMED_EMAIL, PUBMED_API_KEY

logger = logging.getLogger(__name__)

# 設置 Entrez email
if Entrez:
    Entrez.email = PUBMED_EMAIL
    if PUBMED_API_KEY:
        Entrez.api_key = PUBMED_API_KEY

# 需要排除的通用引用 PMID（非變異特定的文獻）
EXCLUDED_PMIDS = {
    '28492532',  # 大型 TTN 研究，通用引用
    '25741868',  # 通用 TTN 文獻
    '35802134',  # 用戶指定排除
    '34012068',  # 用戶指定排除
    '20301486',  # 用戶指定排除
    # 可以根據需要添加更多
}


class ClinVarParser:
    """
    ClinVar API 解析器（可選功能）
    
    注意：此類嘗試從 ClinVar API 獲取變異資訊，但如果失敗會優雅地退回到標準 PubMed 搜索。
    不需要任何本地文件或手動下載。
    """
    
    # 鹼基互補配對
    COMPLEMENT = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G'
    }
    
    def __init__(self):
        if not Entrez:
            logger.info("Biopython 未安裝 - ClinVar API 功能被禁用（將使用標準 PubMed 搜索）")
        self.use_api = bool(Entrez)
    
    def _get_complement_base(self, base: str) -> str:
        """獲取互補鹼基"""
        return self.COMPLEMENT.get(base.upper(), base)
    
    def parse_variant(self, variant_info: Dict[str, str]) -> Optional[Dict]:
        """
        從 ClinVar API 查詢變異資訊（可選）
        
        Args:
            variant_info: 變異資訊字典（包含 chrom, pos, ref, alt）
        
        Returns:
            如果成功，返回包含 PubMed IDs、疾病類型等的字典
            如果失敗或未啟用，返回 None（系統會使用標準 PubMed 搜索）
        """
        if not self.use_api:
            logger.info("ClinVar API 未啟用，跳過 ClinVar 查詢")
            return None
        
        logger.info(f"嘗試從 ClinVar API 查詢變異: {variant_info['variant_id']}")
        
        try:
            # 使用 E-utilities 搜索 ClinVar
            result = self._search_clinvar_api(variant_info)
            if result:
                return result
        
        except Exception as e:
            logger.warning(f"ClinVar API 查詢失敗: {e}")
        
        # API 失敗或未找到結果時，嘗試使用網頁抓取方法
        if SCRAPING_AVAILABLE:
            logger.info("嘗試使用網頁抓取方法查詢 ClinVar...")
            try:
                result = self._search_clinvar_by_web(variant_info)
                if result:
                    logger.info("網頁抓取成功獲取 ClinVar 資訊")
                    return result
            except Exception as web_e:
                logger.warning(f"網頁抓取也失敗: {web_e}")
        
        logger.info("將使用標準 PubMed 搜索")
        return None
    
    def _search_clinvar_api(self, variant_info: Dict[str, str]) -> Optional[Dict]:
        """使用 NCBI E-utilities API 精確搜索 ClinVar"""
        try:
            # TTN 基因位於負鏈，需要：
            # 1. 位置 -1
            # 2. 轉換為互補鹼基
            
            # 原始值
            original_pos = variant_info['pos']
            original_ref = variant_info['ref'].upper()
            original_alt = variant_info['alt'].upper()
            
            # 負鏈轉換
            negative_strand_pos = original_pos - 1
            negative_strand_ref = self._get_complement_base(original_ref)
            negative_strand_alt = self._get_complement_base(original_alt)
            
            # 構建 ClinVar 搜索查詢（使用負鏈坐標）
            # 格式：NC_000002.12:178527120:A:G
            clinvar_query = f"NC_000002.12:{negative_strand_pos}:{negative_strand_ref}:{negative_strand_alt}"
            
            logger.info(f"ClinVar 搜索（負鏈轉換）:")
            logger.info(f"  原始: chr{variant_info['chrom']}:{original_pos}:{original_ref}>{original_alt}")
            logger.info(f"  轉換: {clinvar_query}")
            
            search_terms = [
                # 精確搜索（加雙引號）
                f'"{clinvar_query}"',
                # 也嘗試不帶 NC 前綴的精確搜索
                # f'"{variant_info["chrom"]}:{negative_strand_pos}:{negative_strand_ref}:{negative_strand_alt}"',
                # # 原始格式（作為後備）
                # f"{variant_info['chrom']}: {original_pos} (GRCh38) AND TTN[gene]",
            ]
            
            for search_term in search_terms:
                logger.info(f"嘗試搜索: {search_term}")
                
                time.sleep(0.34)  # API 限制
                handle = Entrez.esearch(
                    db="clinvar",
                    term=search_term,
                    retmax=5
                )
                search_results = Entrez.read(handle, validate=False)
                handle.close()
                
                id_list = search_results.get("IdList", [])
                
                if not id_list:
                    logger.debug(f"  未找到結果，嘗試下一個格式")
                    continue
                
                logger.info(f"找到 {len(id_list)} 個 ClinVar 記錄")
                
                # 驗證每個記錄，找到真正匹配的變異
                for clinvar_id in id_list:
                    logger.info(f"  檢查 ClinVar ID: {clinvar_id}")
                    
                    try:
                        # Step 1: 先取得 esummary，從中取出 VCV accession
                        # 注意：Entrez.efetch 必須用 VCV accession（如 VCV000012654）才會回傳完整 XML，
                        #       直接用內部 ID（如 12654）會回傳空的 <set/>
                        time.sleep(0.34)
                        handle = Entrez.esummary(db="clinvar", id=clinvar_id)
                        summary = Entrez.read(handle, validate=False)
                        handle.close()
                        
                        vcv_accession = None
                        try:
                            doc_first = summary['DocumentSummarySet']['DocumentSummary'][0]
                            vcv_accession = doc_first.get('accession')  # e.g., VCV000012654
                        except Exception:
                            vcv_accession = None
                        
                        if not vcv_accession:
                            # 後備：依規則手動建構 VCV accession
                            try:
                                vcv_accession = f"VCV{int(clinvar_id):09d}"
                            except Exception:
                                vcv_accession = clinvar_id
                        
                        # Step 2: 用 VCV accession efetch 取得完整 XML（含所有 isoform HGVS）
                        time.sleep(0.34)
                        handle = Entrez.efetch(
                            db="clinvar",
                            id=vcv_accession,
                            rettype="vcv",
                            retmode="xml"
                        )
                        xml_record = handle.read()
                        handle.close()
                        
                        # 檢查是否為空 XML（避免無效驗證）
                        if isinstance(xml_record, bytes):
                            xml_check = xml_record.decode('utf-8', errors='replace')
                        else:
                            xml_check = xml_record
                        xml_is_empty = (len(xml_check) < 500 or '<VariationArchive' not in xml_check)
                        
                        # Step 3: 驗證變異是否匹配（若 XML 為空，跳過 XML 驗證，靠網頁驗證）
                        if not xml_is_empty:
                            if not self._verify_variant_match(xml_record, variant_info):
                                logger.info(f"  ClinVar ID {clinvar_id} 不匹配輸入變異，跳過")
                                continue
                            logger.info(f"  ✓ ClinVar ID {clinvar_id} (acc={vcv_accession}) 匹配輸入變異")
                        else:
                            logger.warning(
                                f"  efetch 回傳空 XML（accession={vcv_accession}），跳過 XML 驗證"
                            )
                        
                        # Step 4: 解析 summary 並獲取 PubMed IDs
                        result = self._parse_clinvar_summary(summary, variant_info, [clinvar_id])
                        if result:
                            result['clinvar_url'] = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/"
                            result['clinvar_variation_id'] = clinvar_id
                            result['clinvar_accession'] = vcv_accession
                            
                            # Step 5: 從 esummary + XML 抽取所有 HGVS（含所有 isoform）與 Protein change
                            summary_hgvs = self._extract_hgvs_from_esummary(summary)
                            xml_hgvs = self._extract_hgvs_from_xml(xml_record) if not xml_is_empty else {
                                'nucleotide': [], 'protein': [], 'protein_change': []
                            }
                            
                            # 結構化來源（XML <HGVSlist>/<ProteinChange> + esummary）的結果。
                            # 這兩個來源都只描述「本變異本身」，不會包含 co-occurring variants。
                            structured_nuc = list(dict.fromkeys(
                                xml_hgvs.get('nucleotide', []) + summary_hgvs.get('nucleotide', [])
                            ))
                            structured_prot = list(dict.fromkeys(
                                xml_hgvs.get('protein', []) + summary_hgvs.get('protein', [])
                            ))
                            structured_pc = list(dict.fromkeys(
                                xml_hgvs.get('protein_change', []) + summary_hgvs.get('protein_change', [])
                            ))
                            has_structured = bool(structured_nuc or structured_prot or structured_pc)
                            
                            existing_nuc = result.get('hgvs_nucleotide') or []
                            existing_prot = result.get('hgvs_protein') or []
                            existing_pc = result.get('protein_change') or []
                            
                            if has_structured:
                                # 結構化來源有資料 → 直接覆寫 web-scrape 版本，
                                # 避免 ClinVar 頁面 SCV submission Comment 區段帶進來的
                                # co-occurring variants 污染 alias 清單。
                                merged_nuc = structured_nuc
                                merged_prot = structured_prot
                                merged_pc = structured_pc
                            else:
                                # XML/esummary 都失敗（罕見：efetch 回傳空 XML 且 esummary
                                # 也沒帶 hgvs 欄位）→ 退而保留 web-scrape fallback。
                                merged_nuc = list(dict.fromkeys(existing_nuc))
                                merged_prot = list(dict.fromkeys(existing_prot))
                                merged_pc = list(dict.fromkeys(existing_pc))
                            
                            if merged_nuc:
                                result['hgvs_nucleotide'] = merged_nuc
                            if merged_prot:
                                result['hgvs_protein'] = merged_prot
                            if merged_pc:
                                result['protein_change'] = merged_pc
                            if merged_nuc or merged_prot:
                                result['hgvs'] = merged_nuc + merged_prot
                            
                            # 也記下 rsID（從 esummary 抽到的）
                            sum_rsid = summary_hgvs.get('rsid')
                            if sum_rsid and not result.get('rsid'):
                                result['rsid'] = sum_rsid
                            
                            logger.info(
                                f"  HGVS 整合: {len(merged_nuc)} 核苷酸, "
                                f"{len(merged_prot)} 蛋白質, {len(merged_pc)} Protein change "
                                f"(source={'structured' if has_structured else 'web-scrape-fallback'}; "
                                f"XML: {len(xml_hgvs.get('nucleotide', []))}/{len(xml_hgvs.get('protein', []))}/"
                                f"{len(xml_hgvs.get('protein_change', []))}, "
                                f"summary: {len(summary_hgvs.get('nucleotide', []))}/{len(summary_hgvs.get('protein', []))}/"
                                f"{len(summary_hgvs.get('protein_change', []))}, "
                                f"web: {len(existing_nuc)}/{len(existing_prot)}/{len(existing_pc)})"
                            )
                            
                            if result.get('pmid_list'):
                                logger.info(f"成功從 ClinVar 提取資訊（含 {len(result['pmid_list'])} 個 PubMed IDs）")
                            else:
                                logger.info(f"成功從 ClinVar 提取資訊（但該記錄無關聯的 PubMed 文章）")
                            return result
                        else:
                            logger.debug(f"  無法解析 ClinVar 記錄 {clinvar_id}")
                            continue
                    
                    except Exception as e:
                        logger.warning(f"  處理記錄 {clinvar_id} 時出錯: {e}")
                        import traceback
                        logger.debug(traceback.format_exc())
                        continue
            
            logger.warning(f"ClinVar 中未找到精確匹配的變異 {variant_info['variant_id']}")
            return None
            
        except Exception as e:
            logger.warning(f"ClinVar API 查詢失敗: {e}")
            return None
    
    def _verify_variant_match_from_page(self, soup: BeautifulSoup, variant_info: Dict[str, str]) -> bool:
        """
        從 ClinVar 網頁驗證變異是否匹配輸入的變異
        
        Args:
            soup: BeautifulSoup 對象
            variant_info: 輸入的變異資訊
        
        Returns:
            True 如果匹配，False 否則
        """
        try:
            # 取得需要匹配的資訊
            original_pos = variant_info['pos']
            original_ref = variant_info['ref'].upper()
            original_alt = variant_info['alt'].upper()
            
            # 負鏈轉換
            negative_strand_pos = original_pos - 1
            negative_strand_ref = self._get_complement_base(original_ref)
            negative_strand_alt = self._get_complement_base(original_alt)
            
            page_text = soup.get_text().upper()
            
            # 檢查位置（兩種可能性）
            pos_match = (str(original_pos) in page_text or 
                        str(negative_strand_pos) in page_text)
            
            # 檢查參考和替代鹼基
            # 在 ClinVar 頁面中可能以不同格式出現
            # 例如：NC_000002.12:g.178621118C>G 或 VCF 格式
            ref_alt_patterns = [
                f"{negative_strand_ref}>{negative_strand_alt}",  # 負鏈
                f"{negative_strand_ref}>{negative_strand_alt}".lower(),
                f"{original_ref}>{original_alt}",  # 原始
                f"{original_ref}>{original_alt}".lower(),
                f"REFERENCEALLELEVCF>{negative_strand_ref}<",
                f"ALTERNATEALLELEVCF>{negative_strand_alt}<",
            ]
            
            variant_match = any(pattern.upper() in page_text for pattern in ref_alt_patterns)
            
            if pos_match and variant_match:
                logger.info(f"    ✓ 頁面變異匹配: pos={original_pos}/{negative_strand_pos}, {original_ref}/{negative_strand_ref}>{original_alt}/{negative_strand_alt}")
                return True
            else:
                logger.debug(f"    變異不匹配: pos_match={pos_match}, variant_match={variant_match}")
                return False
            
        except Exception as e:
            logger.warning(f"驗證頁面變異時出錯: {e}")
            # 如果無法驗證，保守地假設不匹配
            return False
    
    def _verify_variant_match(self, xml_record: str, variant_info: Dict[str, str]) -> bool:
        """
        驗證 ClinVar 記錄是否與輸入的變異匹配
        
        ClinVar VCV XML 真正的變異座標寫在 <SequenceLocation> 標籤的「屬性」上，
        而不是子元素的 text。例子：
            <SequenceLocation Assembly="GRCh38" forDisplay="true"
                Chr="2" positionVCF="178785999"
                referenceAlleleVCF="C" alternateAlleleVCF="A"/>
        
        重要：positionVCF 是 1-based VCF 座標，與輸入的 variant_info['pos'] 直接對應。
              （不需要 -1，那是 SPDI 才需要的轉換）
              referenceAlleleVCF/alternateAlleleVCF 是正鏈表示。
              對於負鏈基因（如 TTN），這會是輸入 ref/alt 的互補鹼基。
        
        驗證策略：
            候選 (pos, ref, alt) 三元組
            1. 正鏈：(original_pos, original_ref, original_alt)
            2. 負鏈：(original_pos, complement(ref), complement(alt))  ← TTN 常見
            只要任一組與 XML 中任一 SequenceLocation 屬性完全一致 → 匹配
        """
        try:
            import re as _re
            
            original_pos = variant_info['pos']
            original_ref = variant_info['ref'].upper()
            original_alt = variant_info['alt'].upper()
            
            neg_ref = self._get_complement_base(original_ref)
            neg_alt = self._get_complement_base(original_alt)
            
            # 候選的 (pos, ref, alt) 組合
            candidates = [
                (str(original_pos), original_ref, original_alt),       # 正鏈
                (str(original_pos), neg_ref, neg_alt),                 # 負鏈
            ]
            
            if isinstance(xml_record, bytes):
                xml_str = xml_record.decode('utf-8', errors='replace')
            else:
                xml_str = xml_record or ''
            
            if not xml_str or '<VariationArchive' not in xml_str:
                logger.debug("  XML 內容為空或非 VCV 格式，無法驗證")
                return False
            
            # 抓出所有 <SequenceLocation ... /> 元素的屬性字串
            seqloc_attrs_list: List[Dict[str, str]] = []
            for m in _re.finditer(r'<SequenceLocation\s+([^/>]+)/?>', xml_str):
                attrs_str = m.group(1)
                attrs: Dict[str, str] = {}
                for kv in _re.finditer(r'(\w+)\s*=\s*"([^"]*)"', attrs_str):
                    attrs[kv.group(1)] = kv.group(2)
                seqloc_attrs_list.append(attrs)
            
            if not seqloc_attrs_list:
                logger.debug("  XML 中沒有 SequenceLocation 標籤")
                return False
            
            # 逐一比對：只看 GRCh38 的 forDisplay/current 紀錄優先，但其他也可以接受
            # 先收集 GRCh38 的紀錄；若沒有匹配再 fallback 到任何 assembly
            grch38_locs = [a for a in seqloc_attrs_list if a.get('Assembly') == 'GRCh38']
            check_order = grch38_locs + [a for a in seqloc_attrs_list if a not in grch38_locs]
            
            for attrs in check_order:
                xml_pos = attrs.get('positionVCF') or attrs.get('start')
                xml_ref = (attrs.get('referenceAlleleVCF') or '').upper()
                xml_alt = (attrs.get('alternateAlleleVCF') or '').upper()
                
                if not xml_pos or not xml_ref or not xml_alt:
                    continue
                
                for cand_pos, cand_ref, cand_alt in candidates:
                    if (xml_pos == cand_pos and 
                        xml_ref == cand_ref.upper() and 
                        xml_alt == cand_alt.upper()):
                        logger.info(
                            f"  ✓ 變異匹配: assembly={attrs.get('Assembly', '?')}, "
                            f"pos={xml_pos}, {xml_ref}>{xml_alt}"
                        )
                        return True
            
            # 沒匹配上 - 印出 XML 中實際的內容以利除錯
            sample_locs = [
                (a.get('Assembly', '?'), a.get('positionVCF', a.get('start', '?')),
                 a.get('referenceAlleleVCF', '?'), a.get('alternateAlleleVCF', '?'))
                for a in seqloc_attrs_list[:3]
            ]
            logger.debug(
                f"  XML 中無匹配紀錄。候選: {candidates}, "
                f"XML 中前幾筆 SequenceLocation: {sample_locs}"
            )
            return False
            
        except Exception as e:
            logger.debug(f"驗證變異匹配時出錯: {e}")
            return False
    
    def _parse_clinvar_summary(self, summary: Dict, variant_info: Dict[str, str], id_list: List) -> Optional[Dict]:
        """解析 ClinVar summary 結果並獲取 PubMed IDs"""
        try:
            if not summary or 'DocumentSummarySet' not in summary:
                return None
            
            doc_set = summary['DocumentSummarySet']
            if 'DocumentSummary' not in doc_set or not doc_set['DocumentSummary']:
                return None
            
            doc = doc_set['DocumentSummary'][0]
            
            result = {
                'pmid_list': [],
                'conditions': [],
                'variant_impact': 'Not specified',
                'clinical_significance': 'Not specified',
                'review_status': 'Not specified'
            }
            
            # 提取臨床意義
            if 'clinical_significance' in doc:
                sig = doc['clinical_significance']
                if 'description' in sig:
                    result['clinical_significance'] = sig['description']
            
            # 提取疾病/表型
            if 'trait_set' in doc:
                traits = doc['trait_set']
                if isinstance(traits, list):
                    for trait in traits:
                        if 'trait_name' in trait:
                            result['conditions'].append(trait['trait_name'])
            
            # 判斷影響類型
            result['variant_impact'] = self._determine_impact_type(result['conditions'])
            
            # 嘗試從 ClinVar 記錄獲取相關的 PubMed IDs
            # 優先使用網頁抓取方法
            if SCRAPING_AVAILABLE and id_list:
                try:
                    logger.info("優先使用網頁抓取方法獲取 PubMed IDs 和 dbSNP ID...")
                    scraped_pmids = self._scrape_pmids_from_clinvar_page(id_list[0], result)
                    if scraped_pmids:
                        result['pmid_list'] = scraped_pmids
                        logger.info(f"從網頁抓取到 {len(scraped_pmids)} 個 PubMed IDs")
                    else:
                        logger.info(f"網頁抓取未找到 PubMed IDs，嘗試 elink API 作為備用...")
                        # 如果網頁抓取沒有結果，嘗試 elink API 作為備用
                        if Entrez:
                            try:
                                time.sleep(0.34)
                                handle = Entrez.elink(
                                    dbfrom="clinvar",
                                    db="pubmed",
                                    id=id_list[0]
                                )
                                link_results = Entrez.read(handle, validate=False)
                                handle.close()
                                
                                # 提取 PubMed IDs
                                if link_results and len(link_results) > 0:
                                    linksetdb = link_results[0].get('LinkSetDb', [])
                                    for linkdb in linksetdb:
                                        if linkdb.get('LinkName') == 'clinvar_pubmed':
                                            links = linkdb.get('Link', [])
                                            pmids = [link['Id'] for link in links]
                                            
                                            # 過濾排除的 PMID
                                            original_count = len(pmids)
                                            pmids = [pmid for pmid in pmids if pmid not in EXCLUDED_PMIDS]
                                            excluded_count = original_count - len(pmids)
                                            
                                            if excluded_count > 0:
                                                logger.info(f"排除 {excluded_count} 個通用引用 PMID")
                                            
                                            result['pmid_list'] = pmids[:20]  # 限制最多20篇
                                            logger.info(f"從 ClinVar elink 提取到 {len(result['pmid_list'])} 個 PubMed IDs")
                                            break
                            except Exception as api_error:
                                logger.warning(f"elink API 備用方法失敗: {api_error}")
                except Exception as scrape_error:
                    logger.warning(f"網頁抓取失敗: {scrape_error}")
                    logger.debug(f"詳細錯誤: {scrape_error}", exc_info=True)
                    # 如果網頁抓取失敗，嘗試 elink API 作為備用
                    if Entrez and id_list:
                        try:
                            logger.info("網頁抓取失敗，嘗試 elink API 作為備用...")
                            time.sleep(0.34)
                            handle = Entrez.elink(
                                dbfrom="clinvar",
                                db="pubmed",
                                id=id_list[0]
                            )
                            link_results = Entrez.read(handle, validate=False)
                            handle.close()
                            
                            # 提取 PubMed IDs
                            if link_results and len(link_results) > 0:
                                linksetdb = link_results[0].get('LinkSetDb', [])
                                for linkdb in linksetdb:
                                    if linkdb.get('LinkName') == 'clinvar_pubmed':
                                        links = linkdb.get('Link', [])
                                        pmids = [link['Id'] for link in links]
                                        
                                        # 過濾排除的 PMID
                                        original_count = len(pmids)
                                        pmids = [pmid for pmid in pmids if pmid not in EXCLUDED_PMIDS]
                                        excluded_count = original_count - len(pmids)
                                        
                                        if excluded_count > 0:
                                            logger.info(f"排除 {excluded_count} 個通用引用 PMID")
                                        
                                        result['pmid_list'] = pmids[:20]  # 限制最多20篇
                                        logger.info(f"從 ClinVar elink 提取到 {len(result['pmid_list'])} 個 PubMed IDs")
                                        break
                        except Exception as api_error:
                            logger.warning(f"elink API 備用方法也失敗: {api_error}")
            elif Entrez and id_list:
                # 如果網頁抓取不可用，使用 elink API
                try:
                    logger.info("網頁抓取不可用，使用 elink API...")
                    time.sleep(0.34)
                    handle = Entrez.elink(
                        dbfrom="clinvar",
                        db="pubmed",
                        id=id_list[0]
                    )
                    link_results = Entrez.read(handle, validate=False)
                    handle.close()
                    
                    # 提取 PubMed IDs
                    if link_results and len(link_results) > 0:
                        linksetdb = link_results[0].get('LinkSetDb', [])
                        for linkdb in linksetdb:
                            if linkdb.get('LinkName') == 'clinvar_pubmed':
                                links = linkdb.get('Link', [])
                                pmids = [link['Id'] for link in links]
                                
                                # 過濾排除的 PMID
                                original_count = len(pmids)
                                pmids = [pmid for pmid in pmids if pmid not in EXCLUDED_PMIDS]
                                excluded_count = original_count - len(pmids)
                                
                                if excluded_count > 0:
                                    logger.info(f"排除 {excluded_count} 個通用引用 PMID")
                                
                                result['pmid_list'] = pmids[:20]  # 限制最多20篇
                                logger.info(f"從 ClinVar elink 提取到 {len(result['pmid_list'])} 個 PubMed IDs")
                                break
                except Exception as e:
                    logger.warning(f"從 ClinVar 獲取 PubMed IDs 時出錯: {e}")
                    logger.debug(f"詳細錯誤: {e}", exc_info=True)
            
            logger.info(f"ClinVar 解析成功: {result['clinical_significance']}, {result['variant_impact']}")
            return result
        
        except Exception as e:
            logger.warning(f"解析 ClinVar summary 失敗: {e}")
            return None
    
    def _search_clinvar_by_web(self, variant_info: Dict[str, str]) -> Optional[Dict]:
        """
        通過網頁搜索 ClinVar 變異（當 API 失敗時使用）
        
        Args:
            variant_info: 變異資訊字典
        
        Returns:
            包含 PubMed IDs 和疾病資訊的字典，或 None
        """
        if not SCRAPING_AVAILABLE:
            return None
        
        try:
            # TTN 基因位於負鏈，需要轉換
            original_pos = variant_info['pos']
            original_ref = variant_info['ref'].upper()
            original_alt = variant_info['alt'].upper()
            
            # 負鏈轉換
            negative_strand_pos = original_pos - 1
            negative_strand_ref = self._get_complement_base(original_ref)
            negative_strand_alt = self._get_complement_base(original_alt)
            
            # 構建 ClinVar 搜索 URL
            # 使用精確的座標格式搜索
            search_query = f"NC_000002.12:{negative_strand_pos}:{negative_strand_ref}:{negative_strand_alt}"
            search_url = f'https://www.ncbi.nlm.nih.gov/clinvar/?term="{search_query}"'
            print(search_url)
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            }
            
            logger.info(f"網頁搜索 URL: {search_url}")
            
            # 訪問搜索頁面
            response = requests.get(search_url, headers=headers, timeout=15)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # 檢查是否有搜索結果
            result_count_elem = soup.find('h3', class_='result_count')
            if result_count_elem:
                logger.debug(f"搜索結果: {result_count_elem.get_text(strip=True)}")
            
            # 嘗試多種方式從搜索結果中找到變異連結
            variation_link = None
            candidate_links = []
            
            # 方法1: 尋找 docsum 區塊（ClinVar 搜索結果的標準結構）
            for docsum in soup.find_all(['div', 'dl'], class_=re.compile(r'docsum|rprt')):
                link = docsum.find('a', href=re.compile(r'/clinvar/variation/(\d+)'))
                if link:
                    candidate_links.append(link['href'])
                    logger.debug(f"  在 docsum 中找到連結: {link['href']}")
            
            # 方法2: 尋找包含 "NM_" 的區塊（表示基因轉錄本）
            if not candidate_links:
                for elem in soup.find_all(text=re.compile(r'NM_')):
                    parent = elem.find_parent(['div', 'tr', 'td'])
                    if parent:
                        link = parent.find('a', href=re.compile(r'/clinvar/variation/(\d+)'))
                        if link:
                            candidate_links.append(link['href'])
                            logger.debug(f"  在 NM_ 區塊中找到連結: {link['href']}")
            
            # 方法3: 尋找所有 variation 連結
            if not candidate_links:
                all_variation_links = soup.find_all('a', href=re.compile(r'/clinvar/variation/(\d+)'))
                for link in all_variation_links:
                    href = link['href']
                    # 排除導航連結（通常不包含完整路徑）
                    if 'variation/' in href:
                        candidate_links.append(href)
                        logger.debug(f"  找到 variation 連結: {href}")
            
            # 方法4: 尋找 VCV ID 文字
            if not candidate_links:
                vcv_pattern = re.compile(r'VCV(\d+)')
                for elem in soup.find_all(text=vcv_pattern):
                    vcv_match = vcv_pattern.search(elem)
                    if vcv_match:
                        variation_id = vcv_match.group(1)
                        candidate_links.append(f"/clinvar/variation/{variation_id}/")
                        logger.debug(f"  從文字中找到 VCV{variation_id}")
            
            if candidate_links:
                # 遍歷所有候選連結，驗證變異是否匹配
                for variation_link in candidate_links:
                    try:
                        # 提取 variation ID
                        variation_id_match = re.search(r'/clinvar/variation/(\d+)', variation_link)
                        if not variation_id_match:
                            logger.debug(f"  無法從連結中提取 variation ID: {variation_link}")
                            continue
                        
                        variation_id = variation_id_match.group(1)
                        logger.info(f"  檢查 ClinVar variation ID: {variation_id}")
                        
                        # 詳細頁面 URL
                        detail_url = f"https://www.ncbi.nlm.nih.gov{variation_link}" if not variation_link.startswith('http') else variation_link
                        response = requests.get(detail_url, headers=headers, timeout=15)
                        response.raise_for_status()
                        
                        soup = BeautifulSoup(response.content, 'html.parser')
                        
                        # 驗證此頁面的變異是否與輸入的變異匹配
                        if not self._verify_variant_match_from_page(soup, variant_info):
                            logger.info(f"  變異不匹配，繼續檢查下一個候選...")
                            continue
                        
                        logger.info(f"  ✓ 找到匹配的 ClinVar variation ID: {variation_id}")
                        
                        # 提取資訊
                        result = {
                            'pmid_list': [],
                            'conditions': [],
                            'variant_impact': 'Not specified',
                            'clinical_significance': 'Not specified',
                            'review_status': 'Not specified'
                        }
                        
                        # 提取 dbSNP ID (rsid)
                        rsid = self._extract_rsid_from_soup(soup)
                        if rsid:
                            result['rsid'] = rsid
                            logger.info(f"找到 dbSNP ID: {rsid}")
                        
                        # 提取臨床意義
                        clinical_sig_elem = soup.find('span', class_=re.compile(r'clinical-significance'))
                        if not clinical_sig_elem:
                            # 嘗試其他方式
                            for elem in soup.find_all(['span', 'div', 'strong']):
                                text = elem.get_text(strip=True)
                                if text in ['Pathogenic', 'Likely pathogenic', 'Uncertain significance', 
                                           'Benign', 'Likely benign', 'Conflicting']:
                                    result['clinical_significance'] = text
                                    logger.debug(f"找到臨床意義: {text}")
                                    break
                        else:
                            result['clinical_significance'] = clinical_sig_elem.get_text(strip=True)
                        
                        # 提取疾病/表型（從標題或條件區塊）
                        condition_elems = soup.find_all(['span', 'div'], class_=re.compile(r'condition|trait'))
                        for elem in condition_elems:
                            text = elem.get_text(strip=True)
                            if text and len(text) > 5:  # 過濾太短的文字
                                result['conditions'].append(text)
                        
                        # 如果沒找到，嘗試從頁面標題提取
                        if not result['conditions']:
                            title = soup.find('h1')
                            if title:
                                title_text = title.get_text(strip=True)
                                # 嘗試從標題中提取疾病名稱
                                if 'cardiomyopathy' in title_text.lower() or 'muscular dystrophy' in title_text.lower():
                                    result['conditions'].append(title_text)
                        
                        # 判斷影響類型
                        result['variant_impact'] = self._determine_impact_type(result['conditions'])
                        
                        # 提取 PMID（只從 germline classification 區塊）
                        pmid_list = self._extract_germline_pmids_from_soup(soup)
                        result['pmid_list'] = pmid_list

                        # 儲存 ClinVar 頁面連結
                        result['clinvar_url'] = detail_url
                        result['clinvar_variation_id'] = variation_id

                        logger.info(f"網頁抓取結果: {len(pmid_list)} 個 PMIDs, 臨床意義: {result['clinical_significance']}")

                        return result if result['pmid_list'] or result['clinical_significance'] != 'Not specified' else None
                        
                    except Exception as e:
                        logger.warning(f"  處理候選連結時出錯: {e}")
                        continue
                
                logger.warning("所有候選變異都不匹配輸入的變異")
                return None
            else:
                logger.warning("無法在搜索結果中找到 ClinVar variation 連結")
                logger.debug(f"頁面內容摘要（前500字元）: {soup.get_text()[:500]}")
                return None
            
        except Exception as e:
            logger.warning(f"網頁搜索失敗: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return None
    
    def _extract_germline_pmids_from_soup(self, soup: BeautifulSoup) -> List[str]:
        """
        從已解析的 BeautifulSoup 對象中提取 germline classification 區塊的 PMIDs
        
        Args:
            soup: BeautifulSoup 對象
        
        Returns:
            PMID 列表
        """
        # 尋找 "Citations for germline classification of this variant" 標題
        germline_citation_section = None
        
        # 嘗試多種方式找到這個標題
        for h3 in soup.find_all('h3'):
            if h3.get_text(strip=True).startswith('Citations for germline classification'):
                germline_citation_section = h3
                logger.debug(f"找到 germline citation 標題 (h3)")
                break
        
        if not germline_citation_section:
            for heading in soup.find_all(['h2', 'h3', 'h4']):
                if 'Citations for germline classification' in heading.get_text():
                    germline_citation_section = heading
                    logger.debug(f"找到 germline citation 標題 ({heading.name})")
                    break
        
        if not germline_citation_section:
            germline_citation_section = soup.find(id='new-germline-citation')
            if germline_citation_section:
                logger.debug(f"找到 germline citation 區塊 (by id)")
        
        if not germline_citation_section:
            logger.warning(f"無法在 ClinVar 頁面中找到 'Citations for germline classification' 標題")
            # 回退：提取所有 PMID
            pmid_links = soup.find_all('a', href=re.compile(r'pubmed[./](\d+)'))
        else:
            # 找到包含這個標題的父容器
            citation_container = germline_citation_section.find_parent('div', class_='new-citation-table')
            
            if not citation_container:
                citation_container = germline_citation_section
            
            # 在這個容器內查找所有 PubMed 連結
            pmid_links = citation_container.find_all('a', href=re.compile(r'pubmed[./](\d+)'))
            logger.debug(f"在 germline classification 區塊中找到 {len(pmid_links)} 個 PubMed 連結")
        
        seen_pmids = set()
        pmid_list = []
        
        for link in pmid_links:
            pmid_match = re.search(r'pubmed[./](\d+)', link['href'])
            if pmid_match:
                pmid = pmid_match.group(1)
                
                # 避免重複
                if pmid in seen_pmids:
                    continue
                seen_pmids.add(pmid)
                
                # 排除通用引用
                if pmid in EXCLUDED_PMIDS:
                    logger.debug(f"排除通用引用 PMID: {pmid}")
                    continue
                
                pmid_list.append(pmid)
        
        return pmid_list[:20]  # 限制最多20篇
    
    def _scrape_pmids_from_clinvar_page(self, clinvar_id: str, result: Optional[Dict] = None) -> List[str]:
        """
        從 ClinVar 頁面抓取 PubMed IDs（備用方法）
        當 API 失敗或沒有返回結果時使用
        
        只提取 "Citations for germline classification of this variant" 部分之後的 PMID
        
        Args:
            clinvar_id: ClinVar variation ID
            result: 可選的結果字典，如提供會同時提取 rsid 並存入
        
        Returns:
            PMID 列表
        """
        if not SCRAPING_AVAILABLE:
            return []
        
        try:
            url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/"
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            }
            
            response = requests.get(url, headers=headers, timeout=15)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # 同時提取 rsid（如果提供了 result 字典）
            if result is not None:
                rsid = self._extract_rsid_from_soup(soup)
                if rsid:
                    result['rsid'] = rsid
                    logger.info(f"從 ClinVar 頁面提取到 dbSNP ID: {rsid}")
                
                # 同時提取 HGVS 表示法（含 Protein change 單字母多 isoform）
                # ⚠️ 這份 HGVS 來自整頁 regex 掃描，會把 ClinVar 變異頁 SCV submission
                # Comment 區段裡列的「co-occurring variants」也抓進來（e.g. rs138060032
                # 頁面會出現 "TTN c.87491C>T (p.Pro29164Leu); c.69086_69095del10
                # (p.Arg23029ThrfsX9)" — 那些不是本變異）。下游 _parse_clinvar_summary
                # 之後的 Step 5（XML/esummary 整合，line ~227）會根據是否有結構化來源
                # 決定要不要覆寫這份髒結果，所以此處仍寫入 result 作為 fallback。
                hgvs_data = self._extract_hgvs_from_soup(soup)
                if hgvs_data['nucleotide'] or hgvs_data['protein'] or hgvs_data.get('protein_change'):
                    result['hgvs'] = hgvs_data['nucleotide'] + hgvs_data['protein']
                    result['hgvs_nucleotide'] = hgvs_data['nucleotide']
                    result['hgvs_protein'] = hgvs_data['protein']
                    if hgvs_data.get('protein_change'):
                        result['protein_change'] = hgvs_data['protein_change']
            
            # 使用共享的提取方法
            pmid_list = self._extract_germline_pmids_from_soup(soup)
            
            logger.info(f"從 germline classification 區塊提取到 {len(pmid_list)} 個 PMID")
            return pmid_list
            
        except Exception as e:
            logger.debug(f"網頁抓取失敗: {e}")
            return []
    
    def _extract_rsid_from_soup(self, soup: BeautifulSoup) -> Optional[str]:
        """
        從 ClinVar 頁面提取 dbSNP ID (rs number)
        
        Args:
            soup: BeautifulSoup 對象
        
        Returns:
            dbSNP ID (如 rs375159973) 或 None
        """
        try:
            # 方法1: 從 Variant Details 區塊中尋找 (常見格式)
            # 尋找包含 "Identifiers" 或 dbSNP 連結的區塊
            dbsnp_pattern = re.compile(r'rs\d+')
            
            # 嘗試從 dbSNP 連結中提取
            dbsnp_links = soup.find_all('a', href=re.compile(r'dbsnp|ncbi\.nlm\.nih\.gov/snp/rs'))
            for link in dbsnp_links:
                text = link.get_text(strip=True)
                match = dbsnp_pattern.search(text)
                if match:
                    rsid = match.group(0)
                    logger.debug(f"從 dbSNP 連結提取到 rsid: {rsid}")
                    return rsid
                # 也檢查 href 中的 rsid
                href = link.get('href', '')
                match = dbsnp_pattern.search(href)
                if match:
                    rsid = match.group(0)
                    logger.debug(f"從 dbSNP href 提取到 rsid: {rsid}")
                    return rsid
            
            # 方法2: 從頁面文字中尋找 rs 開頭的 ID
            # 尋找 "Identifiers" 區塊
            for heading in soup.find_all(['h3', 'h4', 'dt', 'th']):
                text = heading.get_text(strip=True).lower()
                if 'identifier' in text or 'dbsnp' in text:
                    parent = heading.find_parent(['div', 'dl', 'table', 'section'])
                    if parent:
                        parent_text = parent.get_text()
                        match = dbsnp_pattern.search(parent_text)
                        if match:
                            rsid = match.group(0)
                            logger.debug(f"從 Identifiers 區塊提取到 rsid: {rsid}")
                            return rsid
            
            # 方法3: 從頁面標題或主要資訊區塊中尋找
            # ClinVar 頁面通常在 Variant Details 區塊列出 rsid
            variant_details = soup.find(['div', 'section'], class_=re.compile(r'variant-details|summary'))
            if variant_details:
                details_text = variant_details.get_text()
                match = dbsnp_pattern.search(details_text)
                if match:
                    rsid = match.group(0)
                    logger.debug(f"從 variant details 區塊提取到 rsid: {rsid}")
                    return rsid
            
            # 方法4: 直接在整個頁面中搜索 (作為備用)
            # 找到所有看起來像 rs number 的文字
            full_text = soup.get_text()
            matches = dbsnp_pattern.findall(full_text)
            if matches:
                # 取第一個匹配的 rsid
                rsid = matches[0]
                logger.debug(f"從頁面全文提取到 rsid: {rsid}")
                return rsid
            
            logger.debug("無法在 ClinVar 頁面中找到 dbSNP ID")
            return None
            
        except Exception as e:
            logger.warning(f"提取 dbSNP ID 時出錯: {e}")
            return None
    
    def _extract_hgvs_from_xml(self, xml_record) -> Dict[str, List[str]]:
        """
        從 ClinVar VCV XML 記錄抽取所有 isoform 的 HGVS 與 Protein change。
        
        ClinVar XML 結構（簡化）:
            <HGVSlist>
              <HGVS>
                <NucleotideExpression sequenceAccessionVersion="NM_..." change="c....">
                  <Expression>NM_...:c.....</Expression>
                </NucleotideExpression>
                <ProteinExpression sequenceAccessionVersion="NP_..." change="p....">
                  <Expression>NP_...:p.....</Expression>
                </ProteinExpression>
              </HGVS>
              ...
            </HGVSlist>
            <ProteinChange>W34072R</ProteinChange>
            <ProteinChange>W25132R</ProteinChange>
            ...
        
        Returns:
            包含 'nucleotide'、'protein'、'protein_change' 列表的字典
        """
        result = {
            'nucleotide': [],
            'protein': [],
            'protein_change': [],
        }
        
        try:
            if isinstance(xml_record, bytes):
                xml_str = xml_record.decode('utf-8', errors='replace')
            else:
                xml_str = xml_record or ''
            
            if not xml_str or '<VariationArchive' not in xml_str:
                return result
            
            def _decode(s):
                return (s.replace('&gt;', '>')
                          .replace('&lt;', '<')
                          .replace('&amp;', '&')
                          .replace('&quot;', '"'))
            
            # 1) 解析所有 <Expression>...</Expression> 標籤
            for m in re.finditer(r'<Expression[^>]*>([^<]+)</Expression>', xml_str):
                expr = _decode(m.group(1).strip())
                if not expr:
                    continue
                if re.search(r':[cgnr]\.', expr):
                    if expr not in result['nucleotide']:
                        result['nucleotide'].append(expr)
                elif ':p.' in expr or expr.startswith('p.'):
                    if expr not in result['protein']:
                        result['protein'].append(expr)
            
            # 2) 解析所有 <ProteinChange>...</ProteinChange>（單字母多 isoform）
            for m in re.finditer(r'<ProteinChange[^>]*>([^<]+)</ProteinChange>', xml_str):
                pc = m.group(1).strip()
                if pc and pc not in result['protein_change']:
                    result['protein_change'].append(pc)
            
            # 3) 備用：從 NucleotideExpression/ProteinExpression 的 change 屬性也抽
            for m in re.finditer(r'<NucleotideExpression[^>]*\bchange="([^"]+)"', xml_str):
                change = _decode(m.group(1).strip())
                if change and change not in result['nucleotide']:
                    result['nucleotide'].append(change)
            for m in re.finditer(r'<ProteinExpression[^>]*\bchange="([^"]+)"', xml_str):
                change = _decode(m.group(1).strip())
                if change and change not in result['protein']:
                    result['protein'].append(change)
            
            if result['nucleotide'] or result['protein'] or result['protein_change']:
                logger.debug(
                    f"  從 XML 抽取到 HGVS: {len(result['nucleotide'])} 核苷酸, "
                    f"{len(result['protein'])} 蛋白質, "
                    f"{len(result['protein_change'])} Protein change"
                )
            
            return result
            
        except Exception as e:
            logger.warning(f"從 XML 抽取 HGVS 時出錯: {e}")
            return result
    
    def _extract_hgvs_from_esummary(self, summary) -> Dict[str, List[str]]:
        """
        從 ClinVar esummary 結果抽取 HGVS / Protein change / rsID。
        
        esummary 的主要 HGVS 欄位:
            - title: "NM_001267550.2(TTN):c.107840T>A (p.Ile35947Asn)"
            - protein_change: "I33379N, I34306N, I35947N, ..."  (單字母多 isoform)
            - variation_set[0].variation_name: 同 title
            - variation_set[0].cdna_change: "c.107840T>A"
            - variation_set[0].canonical_spdi: "NC_000002.12:178527147:A:T"
            - variation_set[0].aliases: 其他別名 list
            - variation_set[0].variation_xrefs: dbSNP 等 xref
        
        Returns:
            dict with keys 'nucleotide', 'protein', 'protein_change', 'rsid'
        """
        result = {
            'nucleotide': [],
            'protein': [],
            'protein_change': [],
            'rsid': None,
        }
        
        try:
            if not summary or 'DocumentSummarySet' not in summary:
                return result
            docs = summary['DocumentSummarySet'].get('DocumentSummary') or []
            if not docs:
                return result
            doc = docs[0]
            
            nucleotide_pattern = re.compile(
                r'N[MCGR]_\d+\.\d+:[cgnr]\.\d+[A-Za-z]*[>_]?[A-Za-z]*(?:del|ins|dup)?[A-Za-z0-9]*'
            )
            protein_pattern = re.compile(
                r'NP_\d+\.\d+:p\.[A-Za-z]{3}\d+(?:[A-Za-z]{3,4}|\*|Ter)(?:fs)?(?:\*\d+)?'
            )
            bare_nucleotide_pattern = re.compile(
                r'(?<![A-Za-z0-9_:])[cgnr]\.\d+[A-Za-z]*[>_]?[A-Za-z]*(?:del|ins|dup)?[A-Za-z0-9]*'
            )
            bare_protein_three_pattern = re.compile(
                r'(?<![A-Za-z0-9_:])p\.[A-Za-z]{3}\d+(?:[A-Za-z]{3,4}|\*|Ter)(?:fs)?(?:\*\d+)?'
            )
            
            def _decode(s):
                if not s:
                    return ''
                return str(s).replace('&gt;', '>').replace('&lt;', '<').replace('&amp;', '&')
            
            def _add(s, target):
                s = _decode(s).strip()
                s = re.sub(r'[,;)\]]+$', '', s)
                if s and s not in target:
                    target.append(s)
            
            # 1) title
            title = _decode(doc.get('title', ''))
            if title:
                for m in nucleotide_pattern.findall(title):
                    _add(m, result['nucleotide'])
                for m in protein_pattern.findall(title):
                    _add(m, result['protein'])
                for m in bare_nucleotide_pattern.findall(title):
                    _add(m, result['nucleotide'])
                for m in bare_protein_three_pattern.findall(title):
                    _add(m, result['protein'])
            
            # 2) protein_change（逗號分隔字串）
            pc_field = doc.get('protein_change', '')
            if pc_field:
                for item in re.split(r'[,;\s]+', str(pc_field)):
                    item = item.strip()
                    if re.match(r'^[ARNDCQEGHILKMFPSTWYV]\d+[ARNDCQEGHILKMFPSTWYV\*]$', item):
                        if item not in result['protein_change']:
                            result['protein_change'].append(item)
            
            # 3) variation_set
            for vs in (doc.get('variation_set') or []):
                name = _decode(vs.get('variation_name', ''))
                if name:
                    for m in nucleotide_pattern.findall(name):
                        _add(m, result['nucleotide'])
                    for m in protein_pattern.findall(name):
                        _add(m, result['protein'])
                    for m in bare_nucleotide_pattern.findall(name):
                        _add(m, result['nucleotide'])
                    for m in bare_protein_three_pattern.findall(name):
                        _add(m, result['protein'])
                
                cdna = _decode(vs.get('cdna_change', ''))
                if cdna:
                    _add(cdna, result['nucleotide'])
                
                spdi = _decode(vs.get('canonical_spdi', ''))
                if spdi:
                    _add(spdi, result['nucleotide'])
                
                for alias in (vs.get('aliases') or []):
                    a = _decode(alias)
                    if nucleotide_pattern.search(a) or bare_nucleotide_pattern.search(a):
                        _add(a, result['nucleotide'])
                    elif protein_pattern.search(a) or bare_protein_three_pattern.search(a):
                        _add(a, result['protein'])
                
                # rsID
                for xref in (vs.get('variation_xrefs') or []):
                    src = xref.get('db_source', '').lower()
                    if 'dbsnp' in src:
                        rs_id_val = xref.get('db_id', '')
                        if rs_id_val:
                            rs_full = f"rs{rs_id_val}" if not str(rs_id_val).startswith('rs') else str(rs_id_val)
                            result['rsid'] = rs_full
                            break
            
            return result
            
        except Exception as e:
            logger.warning(f"從 esummary 抽取 HGVS 時出錯: {e}")
            return result
    
    def _extract_hgvs_from_soup(self, soup: BeautifulSoup) -> Dict[str, List[str]]:
        """
        從 ClinVar 頁面提取 HGVS 表示法（核苷酸、蛋白質、單字母 Protein change）
        
        Args:
            soup: BeautifulSoup 對象
        
        Returns:
            包含 'nucleotide'、'protein'、'protein_change' 列表的字典
        """
        hgvs_data = {
            'nucleotide': [],      # e.g., NM_001267550.2:c.3208G>A
            'protein': [],         # e.g., NP_001254479.2:p.Glu1070Lys
            'protein_change': [],  # e.g., W34072R, W25132R （單字母多 isoform）
        }
        
        try:
            # HGVS patterns
            # 含 RefSeq 前綴
            nucleotide_pattern = re.compile(
                r'N[MCGR]_\d+\.\d+:[cgnr]\.\d+[A-Za-z]*[>_]?[A-Za-z]*(?:del|ins|dup)?[A-Za-z0-9]*'
            )
            protein_pattern = re.compile(
                r'NP_\d+\.\d+:p\.[A-Za-z]{3}\d+(?:[A-Za-z]{3,4}|\*|Ter)(?:fs)?(?:\*\d+)?'
            )
            # 裸 HGVS（不含 RefSeq 前綴）
            bare_nucleotide_pattern = re.compile(
                r'(?<![A-Za-z0-9_:])[cgnr]\.\d+[A-Za-z]*[>_]?[A-Za-z]*(?:del|ins|dup)?[A-Za-z0-9]*'
            )
            bare_protein_three_pattern = re.compile(
                r'(?<![A-Za-z0-9_:])p\.[A-Za-z]{3}\d+(?:[A-Za-z]{3,4}|\*|Ter)(?:fs)?(?:\*\d+)?'
            )
            # 單字母 protein change (e.g., W34072R)
            single_letter_protein_pattern = re.compile(
                r'(?<![A-Za-z0-9])([ARNDCQEGHILKMFPSTWYV])(\d{2,6})([ARNDCQEGHILKMFPSTWYV\*])(?![A-Za-z0-9])'
            )
            
            full_text = soup.get_text()
            
            # 提取核苷酸 HGVS（含 RefSeq）
            for match in nucleotide_pattern.findall(full_text):
                clean_match = re.sub(r'[,;)\]]+$', '', match)
                if clean_match and clean_match not in hgvs_data['nucleotide']:
                    hgvs_data['nucleotide'].append(clean_match)
            
            # 提取蛋白質 HGVS（含 RefSeq）
            for match in protein_pattern.findall(full_text):
                clean_match = re.sub(r'[,;)\]]+$', '', match)
                if clean_match and clean_match not in hgvs_data['protein']:
                    hgvs_data['protein'].append(clean_match)
            
            # 提取裸 HGVS（不含 RefSeq）
            for match in bare_nucleotide_pattern.findall(full_text):
                clean_match = re.sub(r'[,;)\]]+$', '', match)
                if clean_match and clean_match not in hgvs_data['nucleotide']:
                    hgvs_data['nucleotide'].append(clean_match)
            for match in bare_protein_three_pattern.findall(full_text):
                clean_match = re.sub(r'[,;)\]]+$', '', match)
                if clean_match and clean_match not in hgvs_data['protein']:
                    hgvs_data['protein'].append(clean_match)
            
            # 提取 "Protein change" 區塊（單字母多 isoform）
            protein_change_match = re.search(
                r'Protein\s*change\s*[:\s]*([A-Z\d,\s\*]+?)(?=\n|Other\s*names|Canonical|Global|Allele|Links|HGVS|$)',
                full_text,
                re.IGNORECASE
            )
            if protein_change_match:
                pc_text = protein_change_match.group(1)
                for m in single_letter_protein_pattern.finditer(pc_text):
                    single_letter = f"{m.group(1)}{m.group(2)}{m.group(3)}"
                    if single_letter not in hgvs_data['protein_change']:
                        hgvs_data['protein_change'].append(single_letter)
            
            # 方法2: 從 HGVS / Protein change 區塊額外抽取
            for heading in soup.find_all(['dt', 'th', 'strong', 'b']):
                text = heading.get_text(strip=True).lower()
                if 'hgvs' in text or 'name' in text or 'protein change' in text:
                    parent = heading.find_parent(['div', 'dl', 'table', 'tr', 'dd'])
                    if parent:
                        parent_text = parent.get_text()
                        for m in nucleotide_pattern.findall(parent_text):
                            clean_m = re.sub(r'[,;)\]]+$', '', m)
                            if clean_m and clean_m not in hgvs_data['nucleotide']:
                                hgvs_data['nucleotide'].append(clean_m)
                        for m in protein_pattern.findall(parent_text):
                            clean_m = re.sub(r'[,;)\]]+$', '', m)
                            if clean_m and clean_m not in hgvs_data['protein']:
                                hgvs_data['protein'].append(clean_m)
                        if 'protein change' in text:
                            for m in single_letter_protein_pattern.finditer(parent_text):
                                single = f"{m.group(1)}{m.group(2)}{m.group(3)}"
                                if single not in hgvs_data['protein_change']:
                                    hgvs_data['protein_change'].append(single)
            
            if hgvs_data['nucleotide'] or hgvs_data['protein'] or hgvs_data['protein_change']:
                logger.info(
                    f"從 ClinVar 頁面提取到 HGVS: "
                    f"{len(hgvs_data['nucleotide'])} 核苷酸, "
                    f"{len(hgvs_data['protein'])} 蛋白質, "
                    f"{len(hgvs_data['protein_change'])} Protein change"
                )
            
            return hgvs_data
            
        except Exception as e:
            logger.warning(f"提取 HGVS 時出錯: {e}")
            return hgvs_data
    
    def get_rsid_from_clinvar(self, variant_info: Dict[str, str]) -> Optional[str]:
        """
        從 ClinVar 獲取變異的 dbSNP ID (rs number)
        
        此方法可以獨立調用，用於查詢特定變異的 dbSNP ID 以便後續用於 LitVar 搜索
        
        Args:
            variant_info: 變異資訊字典（包含 chrom, pos, ref, alt）
        
        Returns:
            dbSNP ID (如 rs375159973) 或 None
        """
        if not SCRAPING_AVAILABLE:
            logger.warning("網頁抓取功能不可用，無法獲取 dbSNP ID")
            return None
        
        try:
            # 構建搜索查詢
            original_pos = variant_info['pos']
            original_ref = variant_info['ref'].upper()
            original_alt = variant_info['alt'].upper()
            
            # 負鏈轉換 (TTN 基因位於負鏈)
            negative_strand_pos = original_pos - 1
            negative_strand_ref = self._get_complement_base(original_ref)
            negative_strand_alt = self._get_complement_base(original_alt)
            
            search_query = f"NC_000002.12:{negative_strand_pos}:{negative_strand_ref}:{negative_strand_alt}"
            search_url = f'https://www.ncbi.nlm.nih.gov/clinvar/?term="{search_query}"'
            
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            }
            
            logger.info(f"搜索 ClinVar 以獲取 dbSNP ID: {search_url}")
            
            # 訪問搜索頁面
            response = requests.get(search_url, headers=headers, timeout=15)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # 尋找變異詳情頁面連結
            variation_link = None
            for link in soup.find_all('a', href=re.compile(r'/clinvar/variation/(\d+)')):
                variation_link = link['href']
                break
            
            if not variation_link:
                logger.warning("無法找到 ClinVar 變異詳情頁面")
                return None
            
            # 訪問詳細頁面
            detail_url = f"https://www.ncbi.nlm.nih.gov{variation_link}" if not variation_link.startswith('http') else variation_link
            response = requests.get(detail_url, headers=headers, timeout=15)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # 提取 rsid
            rsid = self._extract_rsid_from_soup(soup)
            
            if rsid:
                logger.info(f"成功從 ClinVar 獲取 dbSNP ID: {rsid}")
            else:
                logger.warning("無法從 ClinVar 頁面提取 dbSNP ID")
            
            return rsid
            
        except Exception as e:
            logger.warning(f"從 ClinVar 獲取 dbSNP ID 失敗: {e}")
            return None
    
    def _determine_impact_type(self, conditions: List[str]) -> str:
        """判斷變異影響類型"""
        text = ' '.join(conditions).lower()
        
        cardiac_keywords = [
            'cardiomyopathy', 'dilated cardiomyopathy', 'cardiac', 
            'heart', 'dcm', 'hypertrophic cardiomyopathy'
        ]
        
        skeletal_keywords = [
            'muscular dystrophy', 'limb-girdle', 'skeletal', 'myopathy',
            'lgmd', 'muscle weakness', 'tibial muscular dystrophy',
            'nemaline myopathy'
        ]
        
        has_cardiac = any(keyword in text for keyword in cardiac_keywords)
        has_skeletal = any(keyword in text for keyword in skeletal_keywords)
        
        if has_cardiac and has_skeletal:
            return 'Both cardiac and skeletal'
        elif has_cardiac:
            return 'Cardiac'
        elif has_skeletal:
            return 'Skeletal muscle'
        else:
            return 'Not specified'
