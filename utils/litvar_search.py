#!/usr/bin/env python3
"""
LitVar2 Integration Module
Uses the LitVar2 API to find publications mentioning specific variants by rsID.
Now with full text fetching support!
"""

import logging
import requests
import time
import re
from typing import Dict, List, Optional
from bs4 import BeautifulSoup
from Bio import Entrez

from config import ENABLE_FULL_TEXT_FETCH, MAX_TEXT_LENGTH

try:
    from config import LITVAR_MAX_RESULTS as _DEFAULT_LITVAR_MAX
except ImportError:
    _DEFAULT_LITVAR_MAX = 50

# Try to import docling for PDF processing
try:
    from config import ENABLE_DOCLING_PDF, DOCLING_MAX_PRIORITY_LENGTH, OUTPUT_DIR, DOWNLOAD_SUPPLEMENTARY_FILES
except ImportError:
    ENABLE_DOCLING_PDF = False
    DOCLING_MAX_PRIORITY_LENGTH = 8000
    OUTPUT_DIR = Path("./outputs")
    DOWNLOAD_SUPPLEMENTARY_FILES = True

DOCLING_AVAILABLE = False
try:
    from utils.docling_pdf_processor import DoclingPDFProcessor, download_pdf_from_pmc
    DOCLING_AVAILABLE = True
except ImportError:
    pass

logger = logging.getLogger(__name__)

# Set Entrez email (required by NCBI)
Entrez.email = "variant_agent@example.com"


class LitVarSearcher:
    def __init__(self, try_full_text: bool = None, use_docling: bool = True):
        self.base_url = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get"
        self.litvar_web_url = "https://www.ncbi.nlm.nih.gov/research/litvar2"
        # Use config setting if not specified
        self.try_full_text = try_full_text if try_full_text is not None else ENABLE_FULL_TEXT_FETCH
        self.max_text_length = MAX_TEXT_LENGTH
        self.use_docling = use_docling and DOCLING_AVAILABLE and ENABLE_DOCLING_PDF
        
        # 初始化 docling 處理器
        self.docling_processor = None
        if self.use_docling:
            try:
                self.docling_processor = DoclingPDFProcessor(
                    max_priority_content_length=DOCLING_MAX_PRIORITY_LENGTH,
                    output_dir=OUTPUT_DIR / "markdown"
                )
                logger.info("LitVar2: Docling PDF 處理器已啟用")
            except Exception as e:
                logger.warning(f"LitVar2: 無法初始化 Docling: {e}")
                self.use_docling = False
    
    def search_by_rsid(self, rsid: str) -> List[str]:
        """
        Search LitVar2 for publications by rsID.
        Returns a list of PMIDs.
        """
        if not rsid or not rsid.startswith('rs'):
            logger.debug(f"Invalid rsID format: {rsid}")
            return []
            
        try:
            variant_id = f"litvar@{rsid}%23%23"
            url = f"{self.base_url}/{variant_id}/publications"
            
            logger.info(f"LitVar2: Searching publications for {rsid}")
            
            response = requests.get(url, timeout=15)
            
            if response.status_code == 400:
                logger.debug(f"LitVar2: Variant {rsid} not found in database")
                return []
            elif response.status_code != 200:
                logger.warning(f"LitVar2 API returned {response.status_code}")
                return []
                
            data = response.json()
            pmids = data.get('pmids', [])
            
            logger.info(f"LitVar2: Found {len(pmids)} publications for {rsid}")
            return [str(pmid) for pmid in pmids]
            
        except Exception as e:
            logger.error(f"LitVar2 search error: {e}")
            return []
    
    def _try_fetch_with_docling(self, pmid: str, variant_aliases: Optional[List[str]] = None) -> Optional[dict]:
        """使用 Docling 處理 PDF（支持變異別名過濾）"""
        if not self.docling_processor:
            return None
            
        try:
            pdf_path = download_pdf_from_pmc(pmid, output_dir=OUTPUT_DIR / "pdfs")
            if not pdf_path:
                return None
                
            result = self.docling_processor.process_pdf_for_llm(
                pdf_path, 
                include_full_text=False, 
                pmid=pmid,
                download_supplementary=DOWNLOAD_SUPPLEMENTARY_FILES,
                variant_aliases=variant_aliases
            )
            if result and result.get('priority_content'):
                return result
            return None
        except Exception as e:
            logger.debug(f"LitVar2: Docling 處理失敗 (PMID {pmid}): {e}")
            return None

    def _try_fetch_full_text(self, pmid: str) -> Optional[str]:
        """
        Try to fetch full text from PMC or Europe PMC.
        """
        # Try PMC first
        pmc_text = self._try_fetch_from_pmc(pmid)
        if pmc_text:
            return pmc_text
        
        # Try Europe PMC
        eupmc_text = self._try_fetch_from_eupmc(pmid)
        if eupmc_text:
            return eupmc_text
        
        return None
    
    def _try_fetch_from_pmc(self, pmid: str) -> Optional[str]:
        """Fetch full text from PubMed Central"""
        try:
            # Get PMC ID from PMID
            handle = Entrez.elink(
                dbfrom="pubmed",
                db="pmc",
                id=pmid,
                linkname="pubmed_pmc"
            )
            result = Entrez.read(handle)
            handle.close()
            
            # Extract PMC ID
            pmc_id = None
            for link_set in result:
                if 'LinkSetDb' in link_set:
                    for link_db in link_set['LinkSetDb']:
                        if link_db.get('LinkName') == 'pubmed_pmc':
                            links = link_db.get('Link', [])
                            if links:
                                pmc_id = links[0]['Id']
                                break
            
            if not pmc_id:
                return None
            
            # Fetch full text from PMC
            handle = Entrez.efetch(
                db="pmc",
                id=pmc_id,
                rettype="xml",
                retmode="xml"
            )
            content = handle.read()
            handle.close()
            
            # Parse XML
            soup = BeautifulSoup(content, 'xml')
            
            # Extract text
            title = soup.find('article-title')
            title_text = title.get_text() if title else ""
            
            abstract = soup.find('abstract')
            abstract_text = abstract.get_text() if abstract else ""
            
            body = soup.find('body')
            body_text = body.get_text() if body else ""
            
            full_text = f"{title_text}\n\n{abstract_text}\n\n{body_text}"
            full_text = self._clean_text(full_text)
            
            if len(full_text) > self.max_text_length:
                full_text = full_text[:self.max_text_length] + "\n\n[Text truncated...]"
            
            return full_text if full_text.strip() else None
            
        except Exception as e:
            logger.debug(f"PMC fetch failed for PMID {pmid}: {e}")
            return None
    
    def _try_fetch_from_eupmc(self, pmid: str) -> Optional[str]:
        """Fetch full text from Europe PMC"""
        try:
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmid}/fullTextXML"
            response = requests.get(url, timeout=15)
            
            if response.status_code != 200:
                return None
            
            soup = BeautifulSoup(response.content, 'xml')
            
            title = soup.find('article-title')
            title_text = title.get_text() if title else ""
            
            abstract = soup.find('abstract')
            abstract_text = abstract.get_text() if abstract else ""
            
            body = soup.find('body')
            body_text = body.get_text() if body else ""
            
            full_text = f"{title_text}\n\n{abstract_text}\n\n{body_text}"
            full_text = self._clean_text(full_text)
            
            if len(full_text) > self.max_text_length:
                full_text = full_text[:self.max_text_length] + "\n\n[Text truncated...]"
            
            return full_text if full_text.strip() else None
            
        except Exception as e:
            logger.debug(f"Europe PMC fetch failed for PMID {pmid}: {e}")
            return None
    
    def _clean_text(self, text: str) -> str:
        """Clean and normalize text"""
        text = re.sub(r'\s+', ' ', text)
        text = re.sub(r'\n\s*\n', '\n\n', text)
        return text.strip()
    
    def fetch_article_details(self, pmids: List[str], variant_aliases: Optional[List[str]] = None) -> List[Dict]:
        """
        Fetch article details from PubMed for a list of PMIDs.
        Now includes full text fetching when enabled!
        
        Args:
            pmids: List of PubMed IDs
            variant_aliases: Optional list of variant aliases for table filtering
        """
        if not pmids:
            return []
            
        articles = []
        
        try:
            # Fetch details using Entrez
            pmid_str = ','.join(pmids)
            handle = Entrez.efetch(db='pubmed', id=pmid_str, rettype='xml', retmode='xml')
            records = Entrez.read(handle)
            handle.close()
            
            for article in records.get('PubmedArticle', []):
                try:
                    medline = article.get('MedlineCitation', {})
                    article_data = medline.get('Article', {})
                    
                    # Extract PMID
                    pmid = str(medline.get('PMID', ''))
                    
                    # Extract title
                    title = article_data.get('ArticleTitle', 'No title')
                    
                    # Extract abstract
                    abstract_data = article_data.get('Abstract', {})
                    abstract_texts = abstract_data.get('AbstractText', [])
                    if isinstance(abstract_texts, list):
                        abstract = ' '.join([str(t) for t in abstract_texts])
                    else:
                        abstract = str(abstract_texts)
                    
                    # Extract journal
                    journal_info = article_data.get('Journal', {})
                    journal = journal_info.get('Title', 'Unknown')
                    
                    # Extract year
                    pub_date = article_data.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
                    year = pub_date.get('Year', '')
                    
                    # Extract authors
                    author_list = article_data.get('AuthorList', [])
                    authors = []
                    for author in author_list[:3]:
                        last = author.get('LastName', '')
                        first = author.get('ForeName', '')
                        if last:
                            authors.append(f"{last} {first}".strip())
                    author_str = ', '.join(authors)
                    if len(author_list) > 3:
                        author_str += ' et al.'
                    
                    # Try to fetch full text if enabled
                    full_text = None
                    has_full_text = False
                    docling_content = None
                    has_docling = False
                    
                    docling_full_markdown = None  # 給 deterministic alias gating 用的完整 markdown
                    if self.try_full_text:
                        # 優先嘗試 Docling（傳遞 variant_aliases 進行表格過濾）
                        if self.use_docling:
                            docling_result = self._try_fetch_with_docling(pmid, variant_aliases)
                            if docling_result:
                                docling_content = docling_result.get('priority_content', '')
                                docling_full_markdown = docling_result.get('full_markdown', '') or None
                                has_docling = True
                                has_full_text = True
                                logger.info(f"LitVar2 PMID {pmid}: 使用 Docling 提取優先內容 ({len(docling_content)} 字元)")
                        
                        # 如果 Docling 沒開或失敗，嘗試傳統方式
                        if not has_docling:
                            full_text = self._try_fetch_full_text(pmid)
                            if full_text:
                                has_full_text = True
                                logger.info(f"LitVar2 PMID {pmid}: 獲取全文 ({len(full_text)} 字元)")
                            else:
                                logger.debug(f"LitVar2 PMID {pmid}: 無法獲取全文，使用摘要")
                    
                    # Use docling > full text or abstract for LLM analysis
                    if docling_content:
                        text_for_llm = docling_content
                    elif full_text:
                        text_for_llm = full_text
                    else:
                        text_for_llm = abstract
                    
                    # gate_text：給 deterministic alias gating 用，盡量完整
                    # （docling 的 full_markdown > full_text > priority content > abstract）。
                    # 注意：這份不會塞進 LLM prompt，所以「完整 markdown 很大」也不會吃 token。
                    gate_text = (
                        docling_full_markdown
                        or full_text
                        or docling_content
                        or abstract
                        or ""
                    )

                    articles.append({
                        'pmid': pmid,
                        'title': title,
                        'abstract': abstract,
                        'full_text': full_text,
                        'docling_content': docling_content,
                        'has_docling': has_docling,
                        'text_for_llm': text_for_llm,  # This is what the LLM uses!
                        'gate_text': gate_text,         # Used only by deterministic alias gate
                        'has_full_text': has_full_text,
                        'year': year,
                        'journal': journal,
                        'authors': author_str,
                        'pubmed_link': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        'source': 'LitVar2',
                    })
                    
                except Exception as e:
                    logger.debug(f"Error parsing article: {e}")
                    continue
            
            # Log summary
            full_text_count = sum(1 for a in articles if a.get('has_full_text'))
            if full_text_count > 0:
                logger.info(f"LitVar2: 其中 {full_text_count}/{len(articles)} 篇獲取了全文")
                    
        except Exception as e:
            logger.error(f"Error fetching article details: {e}")
            
        return articles
    
    def search_multiple_formats(self, variant_info: Dict, clinvar_info: Optional[Dict] = None, variant_aliases: Optional[List[str]] = None, max_results: int = _DEFAULT_LITVAR_MAX) -> List[Dict]:
        """
        Search LitVar2 using rsID from ClinVar info.
        
        Args:
            variant_info: 變異信息
            clinvar_info: ClinVar 信息（包含 rsID）
            variant_aliases: 變異別名列表（從 ClinVar 獲取，用於表格過濾）
            max_results: 最大結果數量
        """
        all_pmids = []
        
        if clinvar_info and clinvar_info.get('rsid'):
            rsid = clinvar_info['rsid']
            logger.info(f"LitVar2: Using rsID {rsid} from ClinVar")
            pmids = self.search_by_rsid(rsid)
            all_pmids.extend(pmids)
        else:
            logger.info("LitVar2: No rsID available from ClinVar, skipping LitVar search")
            return []
        
        if not all_pmids:
            logger.info("LitVar2: No publications found")
            return []
        
        # 用 dict.fromkeys 做保序去重，保留 LitVar API 給的原始排序（通常較相關的在前），
        # 不要用 set() 否則順序會被打亂、超過 max_results 時挑到的不是真正的「前 N 篇」。
        unique_pmids = list(dict.fromkeys(all_pmids))
        total_unique = len(unique_pmids)
        if total_unique > max_results:
            logger.warning(
                f"LitVar2: rsID returned {total_unique} unique PMIDs, "
                f"truncating to LITVAR_MAX_RESULTS={max_results}. "
                f"Increase LITVAR_MAX_RESULTS in config.py to fetch all."
            )
            unique_pmids = unique_pmids[:max_results]
        logger.info(f"LitVar2: Fetching details for {len(unique_pmids)} publications")
        
        articles = self.fetch_article_details(unique_pmids, variant_aliases)
        
        logger.info(f"LitVar2: Successfully retrieved {len(articles)} articles")
        return articles
    
    # Keep old methods for backward compatibility
    def search_variant(self, variant_query: str, max_results: int = 50) -> List[Dict]:
        """Legacy method"""
        if variant_query.startswith('rs'):
            pmids = self.search_by_rsid(variant_query)
            if pmids:
                return self.fetch_article_details(pmids[:max_results])
        return []
    
    def _generate_variant_queries(self, variant_info: Dict, clinvar_info: Optional[Dict] = None) -> List[str]:
        """Legacy method"""
        queries = []
        if clinvar_info and clinvar_info.get('rsid'):
            queries.append(clinvar_info['rsid'])
        return queries
