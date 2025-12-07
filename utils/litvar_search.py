#!/usr/bin/env python3
"""
LitVar Integration Module
Fixed: Handles 404 errors gracefully
"""

import logging
import requests
import time
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

class LitVarSearcher:
    def __init__(self):
        self.base_url = "[https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/api/v1](https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/api/v1)"
        self.litvar_web_url = "[https://www.ncbi.nlm.nih.gov/research/litvar2](https://www.ncbi.nlm.nih.gov/research/litvar2)"
        
    def search_variant(self, variant_query: str, max_results: int = 50) -> List[Dict]:
        try:
            url = f"{self.base_url}/search"
            params = {'query': variant_query, 'limit': max_results}
            
            with requests.Session() as session:
                response = session.get(url, params=params, timeout=10)
                
                # 404 means variant not found in LitVar, which is expected for rare variants
                if response.status_code == 404:
                    return []
                    
                response.raise_for_status()
                data = response.json()
                return self._parse_litvar_response(data, variant_query)
            
        except Exception as e:
            logger.debug(f"LitVar search '{variant_query}' returned no results or error: {e}")
            return []
    
    def search_multiple_formats(self, variant_info: Dict, clinvar_info: Optional[Dict] = None, max_results: int = 50) -> List[Dict]:
        queries = self._generate_variant_queries(variant_info, clinvar_info)
        logger.info(f"LitVar: Searching {len(queries)} formats...")
        
        all_articles = []
        seen_pmids = set()
        
        for query in queries:
            articles = self.search_variant(query, max_results=max_results)
            for article in articles:
                pmid = article.get('pmid')
                if pmid and pmid not in seen_pmids:
                    seen_pmids.add(pmid)
                    all_articles.append(article)
            time.sleep(0.2)
        
        logger.info(f"LitVar: Found {len(all_articles)} unique articles")
        return all_articles
    
    def _generate_variant_queries(self, variant_info: Dict, clinvar_info: Optional[Dict] = None) -> List[str]:
        queries = []
        c, p, r, a = variant_info['chrom'], variant_info['pos'], variant_info['ref'], variant_info['alt']
        
        queries.append(f"chr{c}:{p}{r}>{a}")
        
        if clinvar_info:
            if 'rsid' in clinvar_info: queries.append(clinvar_info['rsid'])
            if 'hgvs' in clinvar_info:
                h = clinvar_info['hgvs']
                if isinstance(h, list): queries.extend(h)
                else: queries.append(h)
                
        return list(set([q for q in queries if q]))

    def _parse_litvar_response(self, data: Dict, variant_query: str) -> List[Dict]:
        articles = []
        for result in data.get('results', []):
            try:
                articles.append({
                    'pmid': result.get('pmid'),
                    'title': result.get('title', 'No title'),
                    'abstract': result.get('abstract', ''),
                    'year': result.get('pub_date', '')[:4] if result.get('pub_date') else '',
                    'journal': result.get('journal', 'Unknown'),
                    'authors': result.get('authors', 'Unknown'),
                    'pubmed_link': f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){result.get('pmid')}/",
                    'source': 'LitVar',
                    'snippet': result.get('snippet', ''),
                })
            except: pass
        return articles