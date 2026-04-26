"""
Enhanced PubMed Search Module with Full Text Support
增強版 PubMed 搜尋模組，支持全文獲取
"""

import logging
import time
import re
import threading
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
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
    PHENOTYPE_CATEGORIES,
    OUTPUT_DIR
)

# Try to import docling for PDF processing
try:
    from config import ENABLE_DOCLING_PDF, DOCLING_MAX_PRIORITY_LENGTH, DOWNLOAD_SUPPLEMENTARY_FILES
except ImportError:
    ENABLE_DOCLING_PDF = False
    DOCLING_MAX_PRIORITY_LENGTH = 8000
    DOWNLOAD_SUPPLEMENTARY_FILES = True

DOCLING_AVAILABLE = False
try:
    from utils.docling_pdf_processor import DoclingPDFProcessor, download_pdf_from_pmc
    DOCLING_AVAILABLE = True
except ImportError:
    pass

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
    
    # Lock for Entrez API calls (biopython is not thread-safe)
    _entrez_lock = threading.Lock()

    def __init__(self, try_full_text: bool = True, max_text_length: int = 8000, use_docling: bool = True):
        """
        Initialize enhanced PubMed searcher

        Args:
            try_full_text: 是否嘗試獲取全文（默認 True）
            max_text_length: 最大文本長度（默認 8000 字元，適合 LLM）
            use_docling: 是否使用 docling 處理 PDF（默認 True）
        """
        self.max_results = PUBMED_MAX_RESULTS
        self.try_full_text = try_full_text
        self.max_text_length = max_text_length
        self.use_docling = use_docling and DOCLING_AVAILABLE and ENABLE_DOCLING_PDF
        
        # 初始化 docling 處理器
        self.docling_processor = None
        if self.use_docling:
            try:
                self.docling_processor = DoclingPDFProcessor(
                    max_priority_content_length=DOCLING_MAX_PRIORITY_LENGTH,
                    output_dir=OUTPUT_DIR / "markdown"
                )
                logger.info("Docling PDF 處理器已啟用")
            except Exception as e:
                logger.warning(f"無法初始化 Docling: {e}")
                self.use_docling = False
        
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
        """從 PubMed Central 獲取全文（使用 lock 保護 Entrez，以支援多執行緒）"""
        try:
            # Entrez 非 thread-safe，使用 class-level lock 序列化所有呼叫
            with self.__class__._entrez_lock:
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
            with self.__class__._entrez_lock:
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
    
    def _try_fetch_with_docling(self, pmid: str, variant_aliases: Optional[List[str]] = None) -> Optional[dict]:
        """
        使用 Docling 下載並處理 PDF
        
        優先提取：Tables、Results、Supplementary Data
        如果提供 variant_aliases，表格會被過濾以只包含目標變異
        
        Args:
            pmid: PubMed ID
            variant_aliases: 目標變異的別名列表（用於過濾表格）
            
        Returns:
            包含優先內容的字典，或 None
        """
        if not self.docling_processor:
            return None
        
        try:
            # 嘗試從 PMC 下載 PDF
            pdf_path = download_pdf_from_pmc(pmid, output_dir=OUTPUT_DIR / "pdfs")
            
            if not pdf_path:
                logger.debug(f"PMID {pmid}: 無法從 PMC 下載 PDF")
                return None
            
            # 使用 docling 處理 PDF（傳遞 variant_aliases 用於表格過濾）
            logger.info(f"PMID {pmid}: 使用 Docling 處理 PDF...")
            if variant_aliases:
                logger.info(f"  使用 {len(variant_aliases)} 個變異別名進行表格過濾")
            result = self.docling_processor.process_pdf_for_llm(
                pdf_path, 
                include_full_text=False, 
                pmid=pmid,
                download_supplementary=DOWNLOAD_SUPPLEMENTARY_FILES,
                variant_aliases=variant_aliases
            )
            
            if result and result.get('priority_content'):
                logger.info(f"PMID {pmid}: Docling 提取成功")
                logger.info(f"  - 表格數量: {result.get('tables_count', 0)}")
                logger.info(f"  - Supplementary 連結: {len(result.get('supplementary_links', []))}")
                logger.info(f"  - 有 Results: {result.get('has_results', False)}")
                return result
            
            return None
            
        except Exception as e:
            logger.debug(f"Docling 處理失敗 (PMID {pmid}): {e}")
            return None
    
    def _clean_text(self, text: str) -> str:
        """清理文本"""
        # 移除多餘空白
        text = re.sub(r'\s+', ' ', text)
        # 移除特殊字符
        text = re.sub(r'[\x00-\x08\x0B-\x0C\x0E-\x1F\x7F]', '', text)
        return text.strip()
    
    def search(self, variant_info: Dict[str, str], pmid_list: List[str] = None, variant_aliases: Optional[List[str]] = None) -> List[Dict]:
        """
        Search PubMed and fetch articles with full text support
        
        Args:
            variant_info: Variant information dictionary
            pmid_list: Optional list of PubMed IDs to fetch directly
            variant_aliases: 變異別名列表（從 ClinVar 獲取，用於表格過濾）
        
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
            
            # 批次獲取所有文章 metadata（一次 API call）
            time.sleep(0.5)
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid_list,
                rettype="medline",
                retmode="xml"
            )
            articles = Entrez.read(handle)
            handle.close()

            # Step 1: 解析所有文章的基本 metadata（本地操作，不需要 I/O）
            metas = []
            for article in articles['PubmedArticle']:
                try:
                    meta = self._parse_article_metadata(article)
                    if meta:
                        metas.append(meta)
                except Exception as e:
                    logger.warning(f"解析 metadata 時出錯: {e}")

            logger.info(f"解析了 {len(metas)} 篇文章 metadata")

            # Step 2: 並行獲取全文 / Docling（I/O bound，可安全並行）
            if self.try_full_text and metas:
                # 最多 4 個 workers 以符合 NCBI 速率限制（API key 最多 10 req/s）
                max_workers = min(len(metas), 4)
                logger.info(f"以 {max_workers} 個 worker 並行獲取全文...")

                def _fetch_worker(meta):
                    try:
                        return self._enrich_with_full_text(meta, variant_aliases)
                    except Exception as e:
                        logger.warning(f"全文獲取失敗 (PMID {meta.get('pmid')}): {e}")
                        return meta  # 回傳僅有 metadata 的版本

                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    future_to_meta = {executor.submit(_fetch_worker, m): m for m in metas}
                    for future in as_completed(future_to_meta):
                        enriched = future.result()
                        if enriched:
                            results.append(enriched)
            else:
                results = metas

            logger.info(f"成功獲取 {len(results)} 篇文章")

            # 統計全文獲取情況
            full_text_count = sum(1 for r in results if r.get('has_full_text'))
            if full_text_count > 0:
                logger.info(f"其中 {full_text_count} 篇獲取了全文")
            
        except Exception as e:
            logger.error(f"PubMed 搜索錯誤: {e}", exc_info=True)
        
        return results
    
    def _parse_article_metadata(self, article: Dict) -> Optional[Dict]:
        """
        Phase 1: 從已下載的 Entrez 結果中解析基本 metadata（純本地操作，無 I/O）。
        傳回的 dict 包含 pmid, title, abstract, doi, authors, journal, year 等欄位，
        full_text / docling_content 欄位留空，等待 _enrich_with_full_text 填入。
        """
        try:
            medline = article['MedlineCitation']
            pmid = str(medline['PMID'])
            article_data = medline['Article']

            title = article_data.get('ArticleTitle', '')

            abstract = article_data.get('Abstract', {}).get('AbstractText', [''])
            if isinstance(abstract, list):
                abstract = ' '.join(str(a) for a in abstract)
            else:
                abstract = str(abstract)

            doi = None
            if 'PubmedData' in article:
                for article_id in article['PubmedData'].get('ArticleIdList', []):
                    if article_id.attributes.get('IdType') == 'doi':
                        doi = str(article_id)
                        break

            authors = []
            author_list = article_data.get('AuthorList', [])
            for author in author_list[:3]:
                last = author.get('LastName', '')
                init = author.get('Initials', '')
                if last:
                    authors.append(f"{last} {init}".strip())

            journal = article_data.get('Journal', {})
            journal_title = journal.get('Title', '')
            pub_date = journal.get('JournalIssue', {}).get('PubDate', {})
            year = pub_date.get('Year', '')

            pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            doi_link = f"https://doi.org/{doi}" if doi else None

            return {
                'pmid': pmid,
                'title': title,
                'authors': ', '.join(authors) + (' et al.' if len(author_list) > 3 else ''),
                'journal': journal_title,
                'year': year,
                'abstract': abstract,
                'full_text': None,
                'text_for_llm': abstract,  # 先以摘要為預設，_enrich_with_full_text 會覆蓋
                'has_full_text': False,
                'has_docling': False,
                'docling_content': None,
                'doi': doi,
                'pubmed_link': pubmed_link,
                'doi_link': doi_link,
                'phenotype': 'Not specified',
                'inheritance': 'Not specified',
                'age_onset': 'Not specified',
            }

        except Exception as e:
            logger.warning(f"解析 metadata 時出錯: {e}")
            return None

    def _enrich_with_full_text(
        self, meta: Dict, variant_aliases: Optional[List[str]] = None
    ) -> Dict:
        """
        Phase 2: 為已有 metadata 的文章獲取全文 / Docling 內容（I/O bound，可並行）。
        直接修改並回傳傳入的 meta dict。
        """
        pmid = meta['pmid']
        doi = meta.get('doi')

        full_text = None
        has_full_text = False
        docling_content = None
        has_docling = False

        # 優先 Docling（使用 requests 下載 PDF，thread-safe）
        if self.use_docling and self.docling_processor:
            docling_result = self._try_fetch_with_docling(pmid, variant_aliases)
            if docling_result:
                docling_content = docling_result.get('priority_content', '')
                has_docling = True
                has_full_text = True
                logger.info(f"PMID {pmid}: 使用 Docling 提取優先內容 ({len(docling_content)} 字元)")

        # Docling 失敗時改用傳統全文獲取（PMC XML / Europe PMC）
        if not has_docling:
            full_text = self._try_fetch_full_text(pmid, doi)
            if full_text:
                has_full_text = True

        # 優先級：docling > full_text > abstract
        if docling_content:
            text_for_analysis = docling_content
        elif full_text:
            text_for_analysis = full_text
        else:
            text_for_analysis = meta['abstract']

        meta['full_text'] = full_text
        meta['text_for_llm'] = text_for_analysis
        meta['has_full_text'] = has_full_text
        meta['has_docling'] = has_docling
        meta['docling_content'] = docling_content
        return meta

    def _parse_article_enhanced(self, article: Dict, variant_info: Dict[str, str], variant_aliases: Optional[List[str]] = None) -> Dict:
        """解析文章並嘗試獲取全文（保留以向後相容；新流程請使用 _parse_article_metadata + _enrich_with_full_text）"""
        meta = self._parse_article_metadata(article)
        if meta is None:
            return None
        if self.try_full_text:
            return self._enrich_with_full_text(meta, variant_aliases)
        return meta