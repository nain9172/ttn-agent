"""
PubMed Search Module (Updated)
"""

import logging
import time
from typing import Dict, List
try:
    from Bio import Entrez
except ImportError:
    Entrez = None

from config import PUBMED_EMAIL, PUBMED_API_KEY, PUBMED_MAX_RESULTS

logger = logging.getLogger(__name__)

if Entrez:
    Entrez.email = PUBMED_EMAIL
    if PUBMED_API_KEY: Entrez.api_key = PUBMED_API_KEY

class PubMedSearcher:
    def __init__(self, try_full_text=False, max_text_length=8000):
        self.max_results = PUBMED_MAX_RESULTS
    
    def search(self, variant_info: Dict[str, str], pmid_list: List[str] = None) -> List[Dict]:
        if not Entrez: return []
        
        results = []
        target_pmids = pmid_list if pmid_list else []
        
        # If no PMIDs provided, search by variant terms
        if not target_pmids:
            query = self._build_query(variant_info)
            try:
                handle = Entrez.esearch(db="pubmed", term=query, retmax=20)
                record = Entrez.read(handle)
                target_pmids = record["IdList"]
            except Exception as e:
                logger.error(f"PubMed search failed: {e}")
                return []

        if not target_pmids:
            return []

        logger.info(f"Fetching details for {len(target_pmids)} PMIDs...")
        try:
            handle = Entrez.efetch(db="pubmed", id=target_pmids, rettype="medline", retmode="xml")
            articles = Entrez.read(handle)
            
            for art in articles['PubmedArticle']:
                parsed = self._parse_article(art)
                if parsed: results.append(parsed)
                
        except Exception as e:
            logger.error(f"PubMed fetch failed: {e}")
            
        return results

    def _build_query(self, v: Dict) -> str:
        # Construct a strict query first
        c, p, r, a = v['chrom'], v['pos'], v['ref'], v['alt']
        return f"(TTN[Gene]) AND ({p}[Text Word] OR {p-1}[Text Word])"

    def _parse_article(self, article: Dict) -> Dict:
        try:
            medline = article['MedlineCitation']
            pmid = str(medline['PMID'])
            data = medline['Article']
            
            abstract = data.get('Abstract', {}).get('AbstractText', [''])
            if isinstance(abstract, list): abstract = ' '.join(str(x) for x in abstract)
            
            return {
                'pmid': pmid,
                'title': data.get('ArticleTitle', ''),
                'abstract': abstract,
                'authors': 'Et al.', # Simplified for brevity
                'journal': data.get('Journal', {}).get('Title', ''),
                'year': data.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', 'N/A'),
                'pubmed_link': f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/",
                'full_text': abstract # Fallback
            }
        except:
            return None