#!/usr/bin/env python3
"""
Supplementary Files Downloader
從 PubMed Central 或出版商網站下載 supplementary materials
"""

import logging
import requests
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from bs4 import BeautifulSoup
from Bio import Entrez

logger = logging.getLogger(__name__)

# 設定 Entrez email
try:
    from config import PUBMED_EMAIL
    Entrez.email = PUBMED_EMAIL
except ImportError:
    Entrez.email = "user@example.com"


class SupplementaryDownloader:
    """從 PMC 下載並解析 Supplementary Materials"""
    
    def __init__(self, output_dir: Optional[Path] = None):
        """
        初始化下載器
        
        Args:
            output_dir: 輸出目錄
        """
        self.output_dir = output_dir or Path("./outputs/supplementary")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def get_pmc_id_from_pmid(self, pmid: str) -> Optional[str]:
        """從 PMID 獲取 PMC ID"""
        try:
            handle = Entrez.elink(
                dbfrom="pubmed",
                db="pmc",
                id=pmid,
                linkname="pubmed_pmc"
            )
            result = Entrez.read(handle)
            handle.close()
            
            if not result[0]['LinkSetDb']:
                return None
            
            pmc_id = result[0]['LinkSetDb'][0]['Link'][0]['Id']
            return pmc_id
            
        except Exception as e:
            logger.debug(f"Error getting PMC ID for PMID {pmid}: {e}")
            return None
    
    def get_supplementary_from_entrez(self, pmid: str) -> List[Dict[str, str]]:
        """
        使用 Entrez 獲取 PMC 文章的 supplementary files
        
        Args:
            pmid: PubMed ID
            
        Returns:
            List of dicts with 'name', 'url', 'type' keys
        """
        pmc_id = self.get_pmc_id_from_pmid(pmid)
        if not pmc_id:
            return []
        
        try:
            # 使用 Entrez efetch 獲取 PMC 文章的完整 XML
            handle = Entrez.efetch(
                db="pmc",
                id=pmc_id,
                rettype="xml",
                retmode="xml"
            )
            
            xml_content = handle.read()
            handle.close()
            
            # 解析 XML 找到 supplementary materials
            soup = BeautifulSoup(xml_content, 'xml')
            
            supplementary_links = []
            
            # 方法 1: 查找 <supplementary-material> 標籤
            supp_materials = soup.find_all('supplementary-material')
            
            for supp in supp_materials:
                # 獲取文件描述
                label = supp.find('label')
                caption = supp.find('caption')
                
                name = ""
                if label:
                    name = label.get_text(strip=True)
                elif caption:
                    name = caption.get_text(strip=True)[:100]  # 限制長度
                
                # 查找實際的文件連結
                media = supp.find('media')
                if media and media.get('xlink:href'):
                    href = media.get('xlink:href')
                    
                    # 構建完整的 URL
                    # 嘗試多個可能的 URL 格式
                    if not href.startswith('http'):
                        # 優先使用 pmc.ncbi.nlm.nih.gov 域名（更可靠）
                        full_url = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/bin/{href}"
                    else:
                        full_url = href
                    
                    file_type = self._guess_file_type(full_url, name)
                    
                    supplementary_links.append({
                        'name': name or href,
                        'url': full_url,
                        'type': file_type,
                        'pmc_id': pmc_id
                    })
                    
                    logger.info(f"  找到 supplementary: {name or href} ({file_type})")
            
            # 方法 2: 如果沒有找到，嘗試從 OA (Open Access) 包中獲取
            if not supplementary_links:
                logger.debug(f"PMC{pmc_id}: XML 中未找到 supplementary-material 標籤，嘗試其他方法")
                
                # 查找所有帶有 supplementary 相關屬性的 media 標籤
                all_media = soup.find_all('media')
                for media in all_media:
                    href = media.get('xlink:href', '')
                    mimetype = media.get('mimetype', '')
                    
                    # 檢查是否是補充文件
                    if any(keyword in href.lower() for keyword in ['supp', 'supplement', 'additional']):
                        if not href.startswith('http'):
                            # 優先使用 pmc.ncbi.nlm.nih.gov 域名
                            full_url = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/bin/{href}"
                        else:
                            full_url = href
                        
                        file_type = self._guess_file_type(full_url, href)
                        
                        supplementary_links.append({
                            'name': href.split('/')[-1],
                            'url': full_url,
                            'type': file_type,
                            'pmc_id': pmc_id
                        })
            
            logger.info(f"PMC{pmc_id}: 從 Entrez XML 找到 {len(supplementary_links)} 個 supplementary files")
            return supplementary_links
            
        except Exception as e:
            logger.warning(f"Entrez efetch 失敗 for PMC{pmc_id}: {e}")
            return []
    
    def scrape_pmc_supplementary_links(self, pmid: str) -> List[Dict[str, str]]:
        """
        從 PMC 抓取 supplementary files 連結（優先使用 Entrez API）
        
        Args:
            pmid: PubMed ID
            
        Returns:
            List of dicts with 'name', 'url', 'type' keys
        """
        # 優先使用 Entrez API（更可靠）
        logger.info(f"PMID {pmid}: 嘗試使用 Entrez API 獲取 supplementary files...")
        entrez_links = self.get_supplementary_from_entrez(pmid)
        
        if entrez_links:
            logger.info(f"PMID {pmid}: Entrez API 成功找到 {len(entrez_links)} 個 supplementary files")
            return entrez_links
        
        # 如果 Entrez 失敗，回退到網頁抓取
        logger.info(f"PMID {pmid}: Entrez API 未找到，嘗試網頁抓取...")
        
        pmc_id = self.get_pmc_id_from_pmid(pmid)
        if not pmc_id:
            logger.debug(f"PMID {pmid} not in PMC")
            return []
        
        try:
            # 訪問 PMC 文章頁面
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/"
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            }
            
            response = requests.get(pmc_url, headers=headers, timeout=15)
            if response.status_code != 200:
                logger.debug(f"Failed to access PMC{pmc_id}: {response.status_code}")
                return []
            
            soup = BeautifulSoup(response.content, 'html.parser')
            supplementary_links = []
            
            # 方法 1: 查找 "Supplementary Materials" 區塊
            supp_section = soup.find('div', class_='supplementary-material') or \
                          soup.find('div', id='supplementary-material-sec') or \
                          soup.find('section', class_='sec supplementary-material')
            
            if supp_section:
                # 找到所有連結
                for link in supp_section.find_all('a', href=True):
                    href = link['href']
                    name = link.get_text(strip=True)
                    
                    # 轉換相對路徑為絕對路徑
                    if href.startswith('/'):
                        # 修正 /articles/instance/ID/bin/file.xlsx 路徑為正確的 PMC URL
                        if '/articles/instance/' in href and '/bin/' in href:
                            # 提取 /bin/ 之後的文件名
                            file_part = href.split('/bin/')[-1]
                            # 使用正確的 PMC 域名
                            full_url = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/bin/{file_part}"
                        else:
                            full_url = f"https://www.ncbi.nlm.nih.gov{href}"
                    elif href.startswith('http'):
                        full_url = href
                    else:
                        full_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/{href}"
                    
                    # 判斷文件類型
                    file_type = self._guess_file_type(full_url, name)
                    
                    supplementary_links.append({
                        'name': name or f"Supplementary {len(supplementary_links) + 1}",
                        'url': full_url,
                        'type': file_type,
                        'pmc_id': pmc_id
                    })
            
            # 方法 2: 尋找所有實際的補充文件下載連結
            if not supplementary_links:
                for link in soup.find_all('a', href=True):
                    href = link['href']
                    text = link.get_text(strip=True).lower()
                    
                    # 只抓取實際的文件下載連結（包含文件擴展名或 bin/ 路徑）
                    file_extensions = ['.pdf', '.xlsx', '.xls', '.docx', '.doc', '.csv', '.tsv', '.zip', '.rar', '.txt']
                    has_extension = any(ext in href.lower() for ext in file_extensions)
                    is_bin_path = '/bin/' in href  # PMC 文件通常在 /bin/ 路徑下
                    
                    # 檢查是否為補充文件相關的連結
                    is_supplement_text = any(keyword in text for keyword in ['supplement', 'additional', 'table s', 'figure s', 'data file'])
                    
                    if (has_extension or is_bin_path) and is_supplement_text:
                        # 跳過錨點連結
                        if href.startswith('#'):
                            continue
                            
                        if href.startswith('/'):
                            # 修正 /articles/instance/ID/bin/file.xlsx 路徑為正確的 PMC URL
                            if '/articles/instance/' in href and '/bin/' in href:
                                # 提取 /bin/ 之後的文件名
                                file_part = href.split('/bin/')[-1]
                                # 使用正確的 PMC 域名
                                full_url = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/bin/{file_part}"
                            else:
                                full_url = f"https://www.ncbi.nlm.nih.gov{href}"
                        elif href.startswith('http'):
                            full_url = href
                        else:
                            continue
                        
                        # 從 URL 中提取更好的文件名
                        url_parts = href.split('/')
                        if url_parts:
                            file_name_from_url = url_parts[-1].split('?')[0]  # 移除查詢參數
                            # 如果文件名有意義（不是通用文本），使用它
                            if file_name_from_url and len(file_name_from_url) > 3 and '.' in file_name_from_url:
                                name = file_name_from_url
                            else:
                                name = link.get_text(strip=True)
                        else:
                            name = link.get_text(strip=True)
                        
                        file_type = self._guess_file_type(full_url, name)
                        
                        supplementary_links.append({
                            'name': name or f"Supplementary {len(supplementary_links) + 1}",
                            'url': full_url,
                            'type': file_type,
                            'pmc_id': pmc_id
                        })
            
            # 去重
            seen_urls = set()
            unique_links = []
            for link in supplementary_links:
                if link['url'] not in seen_urls:
                    seen_urls.add(link['url'])
                    unique_links.append(link)
            
            if unique_links:
                logger.info(f"PMID {pmid} (PMC{pmc_id}): 找到 {len(unique_links)} 個 supplementary files")
            
            return unique_links
            
        except Exception as e:
            logger.warning(f"Error scraping supplementary links for PMID {pmid}: {e}")
            return []
    
    def _guess_file_type(self, url: str, name: str) -> str:
        """猜測文件類型"""
        text = (url + " " + name).lower()
        
        if any(ext in text for ext in ['.xlsx', '.xls', 'excel']):
            return 'excel'
        elif any(ext in text for ext in ['.csv', '.tsv']):
            return 'csv'
        elif '.pdf' in text:
            return 'pdf'
        elif any(ext in text for ext in ['.docx', '.doc']):
            return 'word'
        elif any(keyword in text for keyword in ['table', 'data']):
            return 'table'
        elif any(keyword in text for keyword in ['figure', 'image']):
            return 'figure'
        else:
            return 'other'
    
    def download_supplementary_file(
        self, 
        url: str, 
        filename: Optional[str] = None,
        pmid: Optional[str] = None
    ) -> Optional[Path]:
        """
        下載 supplementary file
        
        Args:
            url: 文件 URL
            filename: 保存的文件名（可選）
            pmid: PubMed ID（用於組織文件夾）
            
        Returns:
            下載文件的路徑，失敗返回 None
        """
        try:
            headers = {
                'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
            }
            
            # 如果是 PMC 文件，添加 Referer header（避免 403 錯誤）
            if 'pmc.ncbi.nlm.nih.gov' in url or '/pmc/' in url:
                # 從 URL 中提取 PMC ID
                pmc_match = re.search(r'PMC(\d+)', url)
                if pmc_match:
                    pmc_id = pmc_match.group(1)
                    headers['Referer'] = f'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/'
            
            logger.debug(f"Attempting to download: {url}")
            response = requests.get(url, headers=headers, timeout=30, allow_redirects=True)
            
            if response.status_code != 200:
                logger.warning(f"HTTP {response.status_code} when downloading {url}")
                return None
            
            # 檢查是否為 HTML 頁面（可能是錯誤頁面）
            content_type = response.headers.get('content-type', '').lower()
            if 'text/html' in content_type:
                logger.warning(f"下載的是 HTML 頁面而非文件: {url}")
                return None
            
            # 決定文件名
            if not filename:
                # 從 URL 或 Content-Disposition 獲取
                if 'content-disposition' in response.headers:
                    cd = response.headers['content-disposition']
                    filename_match = re.findall('filename=(.+)', cd)
                    if filename_match:
                        filename = filename_match[0].strip('"\'')
                
                if not filename:
                    # 從 URL 提取
                    filename = url.split('/')[-1].split('?')[0]
                    if not filename or len(filename) < 3:
                        filename = f"supplementary_{hash(url) % 10000}.dat"
            
            # 創建子目錄（如果有 PMID）
            if pmid:
                save_dir = self.output_dir / f"pmid_{pmid}"
                save_dir.mkdir(parents=True, exist_ok=True)
            else:
                save_dir = self.output_dir
            
            save_path = save_dir / filename
            
            # 保存文件
            with open(save_path, 'wb') as f:
                f.write(response.content)
            
            logger.info(f"Downloaded supplementary file: {save_path} ({len(response.content)} bytes)")
            return save_path
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"網路錯誤下載 {url}: {e}")
            return None
        except Exception as e:
            logger.warning(f"下載錯誤 {url}: {type(e).__name__}: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return None
    
    def download_all_supplementary_files(self, pmid: str) -> List[Dict[str, any]]:
        """
        下載某篇文章的所有 supplementary files
        
        Args:
            pmid: PubMed ID
            
        Returns:
            List of dicts with 'name', 'url', 'type', 'local_path' keys
        """
        links = self.scrape_pmc_supplementary_links(pmid)
        
        results = []
        for link in links:
            local_path = self.download_supplementary_file(
                link['url'],
                pmid=pmid
            )
            
            results.append({
                **link,
                'local_path': str(local_path) if local_path else None,
                'downloaded': local_path is not None
            })
        
        return results
    
    def extract_tables_from_excel(self, excel_path: Path) -> List[str]:
        """
        從 Excel 文件中提取表格（轉為 Markdown）
        
        Args:
            excel_path: Excel 文件路徑
            
        Returns:
            List of Markdown table strings
        """
        try:
            import pandas as pd
            
            # 讀取所有 sheets
            excel_file = pd.ExcelFile(excel_path)
            tables = []
            
            for sheet_name in excel_file.sheet_names:
                df = pd.read_excel(excel_file, sheet_name=sheet_name)
                
                # 轉換為 Markdown
                markdown_table = df.to_markdown(index=False)
                
                tables.append(f"### {sheet_name}\n\n{markdown_table}\n")
            
            logger.info(f"Extracted {len(tables)} tables from {excel_path}")
            return tables
            
        except Exception as e:
            logger.warning(f"Error extracting tables from {excel_path}: {e}")
            return []


# 便捷函數
def get_supplementary_info_for_pmid(pmid: str) -> Dict[str, any]:
    """
    獲取某篇文章的 supplementary 信息（不下載文件）
    
    Args:
        pmid: PubMed ID
        
    Returns:
        Dict with 'links' and 'summary' keys
    """
    downloader = SupplementaryDownloader()
    links = downloader.scrape_pmc_supplementary_links(pmid)
    
    # 統計
    type_counts = {}
    for link in links:
        file_type = link['type']
        type_counts[file_type] = type_counts.get(file_type, 0) + 1
    
    return {
        'pmid': pmid,
        'links': links,
        'total_count': len(links),
        'type_distribution': type_counts
    }


if __name__ == "__main__":
    # 測試用
    import sys
    import json
    
    logging.basicConfig(level=logging.INFO)
    
    if len(sys.argv) > 1:
        pmid = sys.argv[1]
        
        print(f"\n=== 抓取 PMID {pmid} 的 Supplementary Files ===\n")
        
        info = get_supplementary_info_for_pmid(pmid)
        
        print(f"找到 {info['total_count']} 個 supplementary files:")
        print(json.dumps(info, indent=2, ensure_ascii=False))
        
        # 可選：下載所有文件
        download = input("\n是否下載所有文件？ (y/n): ")
        if download.lower() == 'y':
            downloader = SupplementaryDownloader()
            results = downloader.download_all_supplementary_files(pmid)
            print(f"\n下載完成！{sum(1 for r in results if r['downloaded'])} / {len(results)} 成功")
    else:
        print("用法: python supplementary_downloader.py <PMID>")
        print("範例: python supplementary_downloader.py 29255378")

