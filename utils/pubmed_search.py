"""
PubMed Search Module
Searches PubMed for variant-related literature
"""

import logging
import time
from typing import Dict, List

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


class PubMedSearcher:
    """PubMed literature searcher"""
    
    def __init__(self):
        self.max_results = PUBMED_MAX_RESULTS
        if not Entrez:
            logger.warning("Biopython not available - PubMed search disabled")
    
    def search(self, variant_info: Dict[str, str]) -> List[Dict]:
        """
        Search PubMed for variant-related articles
        
        Args:
            variant_info: Variant information dictionary
        
        Returns:
            List of article dictionaries with extracted information
        """
        if not Entrez:
            logger.error("Biopython not installed - cannot search PubMed")
            return []
        
        logger.info("Searching PubMed...")
        results = []
        
        # Build search query
        query = self._build_query(variant_info)
        logger.info(f"Search query: {query}")
        
        try:
            # Search PubMed
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=self.max_results,
                sort="relevance"
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            pmid_list = search_results["IdList"]
            logger.info(f"Found {len(pmid_list)} articles")
            
            if not pmid_list:
                return results
            
            # Fetch article details
            time.sleep(0.5)  # Be nice to NCBI servers
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid_list,
                rettype="medline",
                retmode="xml"
            )
            articles = Entrez.read(handle)
            handle.close()
            
            # Parse articles
            for article in articles['PubmedArticle']:
                try:
                    parsed = self._parse_article(article, variant_info)
                    if parsed:
                        results.append(parsed)
                except Exception as e:
                    logger.warning(f"Error parsing article: {e}")
                    continue
            
            logger.info(f"Successfully parsed {len(results)} articles")
            
        except Exception as e:
            logger.error(f"PubMed search error: {e}", exc_info=True)
        
        return results
    
    def _build_query(self, variant_info: Dict[str, str]) -> str:
        """Build PubMed search query"""
        # Base query with TTN
        query_parts = [
            "TTN[Gene]",
            "titin",
            "cardiomyopathy",
            "myopathy"
        ]
        
        # Add position info
        query_parts.append(f"chr{variant_info['chrom']}:{variant_info['pos']}")
        
        # Combine with OR
        query = " OR ".join(query_parts)
        
        # Limit to last 10 years for relevance
        query += " AND (\"2014\"[PDAT] : \"3000\"[PDAT])"
        
        return query
    
    def _parse_article(self, article: Dict, variant_info: Dict[str, str]) -> Dict:
        """Parse PubMed article and extract relevant information"""
        try:
            medline = article['MedlineCitation']
            
            # Extract basic info
            pmid = str(medline['PMID'])
            article_data = medline['Article']
            
            title = article_data.get('ArticleTitle', '')
            abstract = article_data.get('Abstract', {}).get('AbstractText', [''])
            if isinstance(abstract, list):
                abstract = ' '.join(str(a) for a in abstract)
            else:
                abstract = str(abstract)
            
            # Extract authors
            authors = []
            author_list = article_data.get('AuthorList', [])
            for author in author_list[:3]:  # First 3 authors
                last = author.get('LastName', '')
                init = author.get('Initials', '')
                if last:
                    authors.append(f"{last} {init}".strip())
            
            # Extract journal and year
            journal = article_data.get('Journal', {})
            journal_title = journal.get('Title', '')
            pub_date = journal.get('JournalIssue', {}).get('PubDate', {})
            year = pub_date.get('Year', '')
            
            # Extract phenotype
            phenotype = self._extract_phenotype(title + ' ' + abstract)
            
            # Extract inheritance pattern
            inheritance = self._extract_inheritance(title + ' ' + abstract)
            
            # Extract age of onset
            age_onset = self._extract_age_onset(title + ' ' + abstract)
            
            # Build PubMed link
            pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            
            return {
                'pmid': pmid,
                'title': title,
                'authors': ', '.join(authors) + (' et al.' if len(author_list) > 3 else ''),
                'journal': journal_title,
                'year': year,
                'abstract': abstract[:500] + '...' if len(abstract) > 500 else abstract,
                'phenotype': phenotype,
                'inheritance': inheritance,
                'age_onset': age_onset,
                'pubmed_link': pubmed_link
            }
            
        except Exception as e:
            logger.warning(f"Error parsing article: {e}")
            return None
    
    def _extract_phenotype(self, text: str) -> str:
        """Extract phenotype from text"""
        text = text.lower()
        
        for category, keywords in PHENOTYPE_CATEGORIES.items():
            for keyword in keywords:
                if keyword.lower() in text:
                    if category == "heart":
                        return "Cardiac"
                    elif category == "skeletal_muscle":
                        return "Skeletal muscle"
                    elif category == "both":
                        return "Both cardiac and skeletal"
        
        return "Not specified"
    
    def _extract_inheritance(self, text: str) -> str:
        """Extract inheritance pattern from text"""
        text = text.lower()
        
        patterns = {
            "Autosomal dominant": ["autosomal dominant", "ad inheritance"],
            "Autosomal recessive": ["autosomal recessive", "ar inheritance"],
            "X-linked": ["x-linked", "x linked"],
            "De novo": ["de novo", "sporadic"]
        }
        
        for pattern, keywords in patterns.items():
            for keyword in keywords:
                if keyword in text:
                    return pattern
        
        return "Not specified"
    
    def _extract_age_onset(self, text: str) -> str:
        """Extract age of onset from text"""
        text = text.lower()
        
        age_patterns = {
            "Congenital": ["congenital", "birth", "neonatal"],
            "Infantile": ["infant", "infancy"],
            "Childhood": ["childhood", "pediatric", "juvenile"],
            "Adult-onset": ["adult onset", "adulthood"],
            "Late-onset": ["late onset", "elderly"]
        }
        
        for age, keywords in age_patterns.items():
            for keyword in keywords:
                if keyword in text:
                    return age
        
        return "Not specified"