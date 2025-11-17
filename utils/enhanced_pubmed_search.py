"""
Enhanced PubMed Search Module with Full Text Support
增強版 PubMed 搜尋模組，支持全文獲取
"""

import logging
import time
import re
import requests
from typing import Dict, List, Optional
from bs4 import BeautifulSoup

try:
    from Bio import Entrez
except ImportError:
    print("Warning: Biopython not installed. PubMed search will not work.")
    print("Install with: pip install biopython")
    Entrez = None

from config import (
    PUBMED_EMAIL,
    PUBMED_API_KEY,
    PUBMED_MAX_RESULTS,
    PHENOTYPE_CATEGORIES
)

logger = logging.getLogger(__name__)

# Set Entrez email (required by NCBI)
if Entrez:
    Entrez.email = PUBMED_EMAIL
    if PUBMED_API_KEY:
        Entrez.api_key = PUBMED_API_KEY


class EnhancedPubMedSearcher:
    """Enhanced PubMed literature searcher with full text support"""
    
    # 鹼基互補配對（TTN 基因位於負鏈）
    COMPLEMENT = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G'
    }
    
    def __init__(self, try_full_text: bool = True, max_text_length: int = 8000):
        """
        Initialize enhanced PubMed searcher
        
        Args:
            try_full_text: 是否嘗試獲取全文（默認 True）
            max_text_length: 最大文本長度（默認 8000 字元，適合 LLM）
        """
        self.max_results = PUBMED_MAX_RESULTS
        self.try_full_text = try_full_text
        self.max_text_length = max_text_length
        
        if not Entrez:
            logger.warning("Biopython not available - PubMed search disabled")
    
    def _get_complement_base(self, base: str) -> str:
        """獲取互補鹼基（用於負鏈基因）"""
        return self.COMPLEMENT.get(base.upper(), base)
    
    def _get_doi_from_pubmed(self, pmid: str) -> Optional[str]:
        """從 PubMed 獲取 DOI"""
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid,
                rettype="medline",
                retmode="xml"
            )
            article = Entrez.read(handle)
            handle.close()
            
            # 提取 DOI
            if article and 'PubmedArticle' in article:
                article_ids = article['PubmedArticle'][0]['PubmedData']['ArticleIdList']
                for article_id in article_ids:
                    if article_id.attributes['IdType'] == 'doi':
                        return str(article_id)
            
            return None
            
        except Exception as e:
            logger.debug(f"Error getting DOI for PMID {pmid}: {e}")
            return None
    
    def _try_fetch_full_text(self, pmid: str, doi: Optional[str] = None) -> Optional[str]:
        """
        嘗試獲取全文
        
        策略：
        1. 使用 PubMed Central (PMC) 開放獲取
        2. 使用 Unpaywall API (如果有 DOI)
        3. 使用 Europe PMC
        
        Returns:
            完整文本或 None（如果無法獲取）
        """
        # 策略 1: 檢查 PubMed Central
        pmc_text = self._try_fetch_from_pmc(pmid)
        if pmc_text:
            logger.info(f"PMID {pmid}: 從 PMC 獲取全文 ({len(pmc_text)} 字元)")
            return pmc_text
        
        # 策略 2: 使用 Europe PMC
        eupmc_text = self._try_fetch_from_eupmc(pmid)
        if eupmc_text:
            logger.info(f"PMID {pmid}: 從 Europe PMC 獲取全文 ({len(eupmc_text)} 字元)")
            return eupmc_text
        
        # 策略 3: 如果有 DOI，嘗試 Unpaywall
        if doi:
            unpaywall_text = self._try_fetch_from_unpaywall(doi)
            if unpaywall_text:
                logger.info(f"PMID {pmid}: 從 Unpaywall 獲取全文 ({len(unpaywall_text)} 字元)")
                return unpaywall_text
        
        logger.debug(f"PMID {pmid}: 無法獲取全文，將使用摘要")
        return None
    
    def _try_fetch_from_pmc(self, pmid: str) -> Optional[str]:
        """從 PubMed Central 獲取全文"""
        try:
            # 先查詢是否在 PMC 中
            handle = Entrez.elink(
                dbfrom="pubmed",
                db="pmc",
                id=pmid,
                linkname="pubmed_pmc"
            )
            result = Entrez.read(handle)
            handle.close()
            
            # 檢查是否有 PMC ID
            if not result[0]['LinkSetDb']:
                return None
            
            pmc_id = result[0]['LinkSetDb'][0]['Link'][0]['Id']
            
            # 獲取全文
            handle = Entrez.efetch(
                db="pmc",
                id=pmc_id,
                rettype="xml"
            )
            
            # 解析 XML
            xml_content = handle.read()
            handle.close()
            
            # 使用 BeautifulSoup 提取文本
            soup = BeautifulSoup(xml_content, 'xml')
            
            # 提取標題
            title = soup.find('article-title')
            title_text = title.get_text() if title else ""
            
            # 提取摘要
            abstract = soup.find('abstract')
            abstract_text = abstract.get_text() if abstract else ""
            
            # 提取正文
            body = soup.find('body')
            body_text = body.get_text() if body else ""
            
            # 組合文本
            full_text = f"{title_text}\n\n{abstract_text}\n\n{body_text}"
            
            # 清理文本
            full_text = self._clean_text(full_text)
            
            # 限制長度
            if len(full_text) > self.max_text_length:
                # 保留標題 + 摘要 + 部分正文
                full_text = full_text[:self.max_text_length] + "\n\n[Text truncated due to length...]"
            
            return full_text if full_text.strip() else None
            
        except Exception as e:
            logger.debug(f"Error fetching from PMC: {e}")
            return None
    
    def _try_fetch_from_eupmc(self, pmid: str) -> Optional[str]:
        """從 Europe PMC 獲取全文"""
        try:
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/fullTextXML"
            response = requests.get(url, timeout=10)
            
            if response.status_code != 200:
                return None
            
            # 解析 XML
            soup = BeautifulSoup(response.content, 'xml')
            
            # 提取文本
            title = soup.find('article-title')
            title_text = title.get_text() if title else ""
            
            abstract = soup.find('abstract')
            abstract_text = abstract.get_text() if abstract else ""
            
            body = soup.find('body')
            body_text = body.get_text() if body else ""
            
            full_text = f"{title_text}\n\n{abstract_text}\n\n{body_text}"
            full_text = self._clean_text(full_text)
            
            # 限制長度
            if len(full_text) > self.max_text_length:
                full_text = full_text[:self.max_text_length] + "\n\n[Text truncated...]"
            
            return full_text if full_text.strip() else None
            
        except Exception as e:
            logger.debug(f"Error fetching from Europe PMC: {e}")
            return None
    
    def _try_fetch_from_unpaywall(self, doi: str) -> Optional[str]:
        """從 Unpaywall API 獲取開放獲取版本"""
        try:
            # Unpaywall API
            email = PUBMED_EMAIL or "your.email@example.com"
            url = f"https://api.unpaywall.org/v2/{doi}?email={email}"
            
            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                return None
            
            data = response.json()
            
            # 檢查是否開放獲取
            if not data.get('is_oa'):
                return None
            
            # 獲取 PDF 或 HTML 連結
            best_oa_location = data.get('best_oa_location')
            if not best_oa_location:
                return None
            
            pdf_url = best_oa_location.get('url_for_pdf')
            landing_url = best_oa_location.get('url_for_landing_page')
            
            # 嘗試獲取文本（這裡簡化處理，實際可能需要 PDF 解析）
            # 暫時返回 None，因為需要額外的 PDF 解析庫
            logger.debug(f"Found OA version at {pdf_url or landing_url}, but PDF parsing not implemented")
            return None
            
        except Exception as e:
            logger.debug(f"Error fetching from Unpaywall: {e}")
            return None
    
    def _clean_text(self, text: str) -> str:
        """清理文本"""
        # 移除多餘空白
        text = re.sub(r'\s+', ' ', text)
        # 移除特殊字符
        text = re.sub(r'[\x00-\x08\x0B-\x0C\x0E-\x1F\x7F]', '', text)
        return text.strip()
    
    def search(self, variant_info: Dict[str, str], pmid_list: List[str] = None) -> List[Dict]:
        """
        Search PubMed and fetch articles with full text support
        
        Args:
            variant_info: Variant information dictionary
            pmid_list: Optional list of PubMed IDs to fetch directly
        
        Returns:
            List of article dictionaries with full text when available
        """
        if not Entrez:
            logger.error("Biopython not installed - cannot search PubMed")
            return []
        
        logger.info("搜索 PubMed（增強版 - 支持全文）...")
        results = []
        
        try:
            # 使用提供的 PMID 列表
            if not pmid_list:
                logger.warning("No PMID list provided")
                return []
            
            logger.info(f"準備獲取 {len(pmid_list)} 篇文章")
            if self.try_full_text:
                logger.info("將嘗試獲取全文（如果可用）")
            
            # 批次獲取基本資訊
            time.sleep(0.5)
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid_list,
                rettype="medline",
                retmode="xml"
            )
            articles = Entrez.read(handle)
            handle.close()
            
            # 處理每篇文章
            for idx, article in enumerate(articles['PubmedArticle'], 1):
                try:
                    logger.info(f"處理文章 {idx}/{len(pmid_list)}...")
                    
                    # 解析基本資訊
                    parsed = self._parse_article_enhanced(article, variant_info)
                    if parsed:
                        results.append(parsed)
                    
                    # 避免請求過快
                    if idx < len(pmid_list):
                        time.sleep(0.5)
                        
                except Exception as e:
                    logger.warning(f"解析文章時出錯: {e}")
                    continue
            
            logger.info(f"成功獲取 {len(results)} 篇文章")
            
            # 統計全文獲取情況
            full_text_count = sum(1 for r in results if r.get('has_full_text'))
            if full_text_count > 0:
                logger.info(f"其中 {full_text_count} 篇獲取了全文")
            
        except Exception as e:
            logger.error(f"PubMed 搜索錯誤: {e}", exc_info=True)
        
        return results
    
    def _parse_article_enhanced(self, article: Dict, variant_info: Dict[str, str]) -> Dict:
        """解析文章並嘗試獲取全文"""
        try:
            medline = article['MedlineCitation']
            
            # 提取基本資訊
            pmid = str(medline['PMID'])
            article_data = medline['Article']
            
            title = article_data.get('ArticleTitle', '')
            
            # 提取摘要
            abstract = article_data.get('Abstract', {}).get('AbstractText', [''])
            if isinstance(abstract, list):
                abstract = ' '.join(str(a) for a in abstract)
            else:
                abstract = str(abstract)
            
            # 提取 DOI
            doi = None
            if 'PubmedData' in article:
                article_ids = article['PubmedData'].get('ArticleIdList', [])
                for article_id in article_ids:
                    if article_id.attributes.get('IdType') == 'doi':
                        doi = str(article_id)
                        break
            
            # 嘗試獲取全文
            full_text = None
            has_full_text = False
            
            if self.try_full_text:
                full_text = self._try_fetch_full_text(pmid, doi)
                if full_text:
                    has_full_text = True
            
            # 使用全文或摘要進行分析
            text_for_analysis = full_text if full_text else abstract
            
            # 提取作者
            authors = []
            author_list = article_data.get('AuthorList', [])
            for author in author_list[:3]:
                last = author.get('LastName', '')
                init = author.get('Initials', '')
                if last:
                    authors.append(f"{last} {init}".strip())
            
            # 提取期刊和年份
            journal = article_data.get('Journal', {})
            journal_title = journal.get('Title', '')
            pub_date = journal.get('JournalIssue', {}).get('PubDate', {})
            year = pub_date.get('Year', '')
            
            # 構建 PubMed 和 DOI 連結
            pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            doi_link = f"https://doi.org/{doi}" if doi else None
            
            return {
                'pmid': pmid,
                'title': title,
                'authors': ', '.join(authors) + (' et al.' if len(author_list) > 3 else ''),
                'journal': journal_title,
                'year': year,
                'abstract': abstract,
                'full_text': full_text,  # 完整文本（如果有）
                'text_for_llm': text_for_analysis,  # 用於 LLM 分析的文本
                'has_full_text': has_full_text,
                'doi': doi,
                'pubmed_link': pubmed_link,
                'doi_link': doi_link,
                'phenotype': 'Not specified',  # 將由 LLM 提取
                'inheritance': 'Not specified',  # 將由 LLM 提取
                'age_onset': 'Not specified'  # 將由 LLM 提取
            }
            
        except Exception as e:
            logger.warning(f"解析文章時出錯: {e}")
            return None