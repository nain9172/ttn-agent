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
            return self._search_clinvar_api(variant_info)
        
        except Exception as e:
            logger.warning(f"ClinVar API 查詢失敗: {e}")
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
                
                logger.info(f"✅ 找到 {len(id_list)} 個 ClinVar 記錄")
                
                # 使用第一個記錄（通常是最相關的）
                clinvar_id = id_list[0]
                logger.info(f"  使用 ClinVar ID: {clinvar_id}")
                
                try:
                    # 獲取 summary 信息
                    time.sleep(0.34)
                    handle = Entrez.esummary(db="clinvar", id=clinvar_id)
                    summary = Entrez.read(handle, validate=False)
                    handle.close()
                    
                    # 解析並獲取 PubMed IDs
                    result = self._parse_clinvar_summary(summary, variant_info, [clinvar_id])
                    if result:
                        # 只要能解析出結果就返回，不管是否有 PubMed IDs
                        if result.get('pmid_list'):
                            logger.info(f"✅ 成功從 ClinVar 提取資訊（含 {len(result['pmid_list'])} 個 PubMed IDs）")
                        else:
                            logger.info(f"✅ 成功從 ClinVar 提取資訊（但該記錄無關聯的 PubMed 文章）")
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
    
    def _verify_variant_match(self, xml_record: str, variant_info: Dict[str, str]) -> bool:
        """驗證 ClinVar 記錄是否與輸入的變異匹配（使用負鏈轉換後的坐標）"""
        try:
            import xml.etree.ElementTree as XMLTree
            
            # 使用負鏈轉換後的坐標進行驗證（因為 ClinVar 中存儲的是負鏈坐標）
            original_pos = variant_info['pos']
            original_ref = variant_info['ref'].upper()
            original_alt = variant_info['alt'].upper()
            
            # 負鏈轉換
            negative_strand_pos = original_pos - 1
            negative_strand_ref = self._get_complement_base(original_ref)
            negative_strand_alt = self._get_complement_base(original_alt)
            
            # 同時檢查原始和轉換後的坐標
            ref = negative_strand_ref
            alt = negative_strand_alt
            pos = str(negative_strand_pos)
            
            # 將字符串轉換為 XML 對象
            if isinstance(xml_record, bytes):
                xml_record = xml_record.decode('utf-8')
            
            root = XMLTree.fromstring(xml_record)
            
            # ClinVar XML 中的變異信息通常在 SimpleAllele 或 VariationArchive 中
            # 查找所有可能包含 allele 信息的標籤
            
            # 檢查位置是否匹配
            position_found = False
            for start_elem in root.iter('Start'):
                if start_elem.text == pos:
                    position_found = True
                    break
            
            if not position_found:
                for pos_elem in root.iter('Position'):
                    if pos in str(pos_elem.text):
                        position_found = True
                        break
            
            if not position_found:
                logger.debug(f"  位置不匹配: 期望 {pos}")
                return False
            
            # 檢查 ref 和 alt alleles
            # ClinVar 使用 ReferenceAllele 和 AlternateAllele 標籤
            ref_found = False
            alt_found = False
            
            for ref_elem in root.iter('ReferenceAlleleVCF'):
                if ref_elem.text and ref_elem.text.upper() == ref:
                    ref_found = True
                    break
            
            for alt_elem in root.iter('AlternateAlleleVCF'):
                if alt_elem.text and alt_elem.text.upper() == alt:
                    alt_found = True
                    break
            
            # 也檢查其他可能的標籤
            if not ref_found:
                text_content = xml_record.upper()
                ref_patterns = [
                    f'REFERENCEALLELE>{ref}<',
                    f'REFERENCEALLELEVCF>{ref}<',
                    f'"REFERENCE":"{ref}"'
                ]
                ref_found = any(pattern in text_content for pattern in ref_patterns)
            
            if not alt_found:
                text_content = xml_record.upper()
                alt_patterns = [
                    f'ALTERNATEALLELE>{alt}<',
                    f'ALTERNATEALLELEVCF>{alt}<',
                    f'"ALTERNATE":"{alt}"'
                ]
                alt_found = any(pattern in text_content for pattern in alt_patterns)
            
            if position_found and ref_found and alt_found:
                logger.info(f"  ✅ 變異匹配: pos={pos}, ref={ref}, alt={alt}")
                return True
            else:
                # 如果負鏈坐標不匹配，也嘗試原始坐標（有些記錄可能使用正鏈）
                logger.debug(f"  負鏈坐標不匹配，嘗試原始坐標...")
                logger.debug(f"  不匹配: pos={position_found}, ref={ref_found}, alt={alt_found}")
                
                # 嘗試原始坐標
                pos_original = str(original_pos)
                ref_original = original_ref
                alt_original = original_alt
                
                position_found_orig = False
                for start_elem in root.iter('Start'):
                    if start_elem.text == pos_original:
                        position_found_orig = True
                        break
                
                if not position_found_orig:
                    for pos_elem in root.iter('Position'):
                        if pos_original in str(pos_elem.text):
                            position_found_orig = True
                            break
                
                ref_found_orig = False
                for ref_elem in root.iter('ReferenceAlleleVCF'):
                    if ref_elem.text and ref_elem.text.upper() == ref_original:
                        ref_found_orig = True
                        break
                
                alt_found_orig = False
                for alt_elem in root.iter('AlternateAlleleVCF'):
                    if alt_elem.text and alt_elem.text.upper() == alt_original:
                        alt_found_orig = True
                        break
                
                if position_found_orig and ref_found_orig and alt_found_orig:
                    logger.info(f"  ✅ 變異匹配（原始坐標）: pos={pos_original}, ref={ref_original}, alt={alt_original}")
                    return True
                else:
                    logger.debug(f"  原始坐標也不匹配")
                    return False
            
        except Exception as e:
            logger.debug(f"驗證變異匹配時出錯: {e}")
            # 如果無法解析，使用簡單的字符串匹配作為後備
            xml_str = xml_record if isinstance(xml_record, str) else xml_record.decode('utf-8')
            return (pos in xml_str and 
                    ref.upper() in xml_str.upper() and 
                    alt.upper() in xml_str.upper())
    
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
            # 方法 1: 使用 elink API（主要方法）
            try:
                if id_list:
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
                        
                        if not result['pmid_list']:
                            logger.info(f"ClinVar 記錄 {id_list[0]} 沒有關聯的 PubMed 文章（elink 方法）")
                            # 嘗試備用方法：網頁抓取
                            if SCRAPING_AVAILABLE:
                                logger.info("嘗試使用網頁抓取方法獲取 PubMed IDs...")
                                scraped_pmids = self._scrape_pmids_from_clinvar_page(id_list[0])
                                if scraped_pmids:
                                    result['pmid_list'] = scraped_pmids
                                    logger.info(f"從網頁抓取到 {len(scraped_pmids)} 個 PubMed IDs")
            except Exception as e:
                logger.warning(f"從 ClinVar 獲取 PubMed IDs 時出錯: {e}")
                logger.debug(f"詳細錯誤: {e}", exc_info=True)
                
                # 如果 API 失敗，嘗試網頁抓取
                if SCRAPING_AVAILABLE and id_list:
                    logger.info("API 方法失敗，嘗試網頁抓取...")
                    try:
                        scraped_pmids = self._scrape_pmids_from_clinvar_page(id_list[0])
                        if scraped_pmids:
                            result['pmid_list'] = scraped_pmids
                            logger.info(f"從網頁抓取到 {len(scraped_pmids)} 個 PubMed IDs")
                    except Exception as scrape_error:
                        logger.debug(f"網頁抓取也失敗: {scrape_error}")
            
            logger.info(f"ClinVar 解析成功: {result['clinical_significance']}, {result['variant_impact']}")
            return result
        
        except Exception as e:
            logger.warning(f"解析 ClinVar summary 失敗: {e}")
            return None
    
    def _scrape_pmids_from_clinvar_page(self, clinvar_id: str) -> List[str]:
        """
        從 ClinVar 頁面抓取 PubMed IDs（備用方法）
        當 API 失敗或沒有返回結果時使用
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
            
            # 尋找所有 PubMed 連結
            pmid_links = soup.find_all('a', href=re.compile(r'pubmed[./](\d+)'))
            
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
            
        except Exception as e:
            logger.debug(f"網頁抓取失敗: {e}")
            return []
    
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
