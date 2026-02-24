"""
Docling PDF Processor Module
使用 Docling 將 PDF 文章轉換為 Markdown，並優先提取關鍵內容

功能：
1. 將 PDF 轉換為結構化 Markdown
2. 擷取 Supplementary Data 連結
3. 轉換文中的表格（包含圖片形式的表格）
4. 優先提取 Results 段落、Tables、Supplementary Data
"""

import logging
import re
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# 嘗試導入 docling
try:
    from docling.document_converter import DocumentConverter
    from docling.datamodel.base_models import InputFormat
    from docling.datamodel.pipeline_options import PdfPipelineOptions
    from docling.backend.pypdfium2_backend import PyPdfiumDocumentBackend
    DOCLING_AVAILABLE = True
except ImportError:
    DOCLING_AVAILABLE = False
    logger.warning("Docling not installed. Install with: pip install docling")


@dataclass
class ExtractedContent:
    """結構化的提取內容"""
    full_markdown: str
    tables: List[str]
    results_section: str
    supplementary_links: List[Dict[str, str]]
    figures: List[str]
    abstract: str
    methods_section: str
    discussion_section: str
    priority_content: str  # 優先處理的內容（tables + results + supplementary）


class DoclingPDFProcessor:
    """
    使用 Docling 處理 PDF 文件的處理器
    
    特點：
    - 轉換 PDF 為 Markdown 格式
    - 智能識別並提取表格（包含圖片形式）
    - 擷取 Supplementary Data 連結
    - 優先處理關鍵段落以節省 GPU 記憶體
    """
    
    # 常見的 section 標題模式
    SECTION_PATTERNS = {
        'abstract': r'(?:^|\n)#+?\s*(?:Abstract|ABSTRACT|摘要)\s*\n',
        'results': r'(?:^|\n)#+?\s*(?:Results?|RESULTS?|結果)\s*\n',
        'methods': r'(?:^|\n)#+?\s*(?:Methods?|METHODS?|Materials?\s+and\s+Methods?|方法)\s*\n',
        'discussion': r'(?:^|\n)#+?\s*(?:Discussion|DISCUSSION|討論)\s*\n',
        'supplementary': r'(?:^|\n)#+?\s*(?:Supplement(?:ary|al)?(?:\s+(?:Data|Materials?|Information))?|Supporting\s+Information|附錄)\s*\n',
        'references': r'(?:^|\n)#+?\s*(?:References?|REFERENCES?|Bibliography|參考文獻)\s*\n',
    }
    
    # Supplementary 連結模式
    SUPPLEMENTARY_LINK_PATTERNS = [
        r'\[([^\]]*[Ss]uppl[^\]]*)\]\(([^)]+)\)',  # Markdown 格式連結
        r'(?:Supplementary|Supporting)\s+(?:Data|Materials?|Information|Table|Figure)\s*(?:S?\d+)?[:\s]*(?:available\s+(?:at|from))?\s*(https?://[^\s\)]+)',
        r'(?:see|refer\s+to)\s+(?:Supplementary|Supporting)\s+(?:Data|Materials?|Information)\s*(?:\(([^)]+)\))?',
    ]
    
    def __init__(
        self,
        max_priority_content_length: int = 8000,
        enable_ocr: bool = True,
        enable_table_extraction: bool = True,
        output_dir: Optional[Path] = None
    ):
        """
        初始化 PDF 處理器
        
        Args:
            max_priority_content_length: 優先內容的最大長度（預設 8000，適合有限 GPU 記憶體）
            enable_ocr: 是否啟用 OCR（用於圖片中的表格）
            enable_table_extraction: 是否啟用表格提取
            output_dir: 輸出目錄（用於保存轉換後的 Markdown）
        """
        if not DOCLING_AVAILABLE:
            raise ImportError(
                "Docling is not installed. Please install with:\n"
                "pip install docling"
            )
        
        self.max_priority_content_length = max_priority_content_length
        self.enable_ocr = enable_ocr
        self.enable_table_extraction = enable_table_extraction
        self.output_dir = output_dir or Path("./outputs/markdown")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 初始化 DocumentConverter
        self._init_converter()
    
    def _init_converter(self):
        """初始化 Docling DocumentConverter"""
        try:
            # 使用簡單的初始化方式（docling 會自動處理 OCR 和表格提取）
            self.converter = DocumentConverter()
            logger.info("Docling DocumentConverter 初始化成功")
        except Exception as e:
            logger.error(f"Docling 初始化失敗: {e}")
            raise
    
    def convert_pdf_to_markdown(self, pdf_path: str) -> str:
        """
        將 PDF 轉換為 Markdown
        
        Args:
            pdf_path: PDF 文件路徑
            
        Returns:
            Markdown 格式的文本
        """
        pdf_name = Path(pdf_path).stem
        md_path = self.output_dir / f"{pdf_name}.md"
        
        # 檢查 Markdown 檔案是否已存在（快取機制）
        if md_path.exists():
            logger.info(f"Markdown 已存在，跳過轉換: {md_path}")
            try:
                with open(md_path, 'r', encoding='utf-8') as f:
                    return f.read()
            except Exception as e:
                logger.warning(f"讀取快取的 Markdown 失敗: {e}，重新轉換")
        
        logger.info(f"開始轉換 PDF: {pdf_path}")
        
        try:
            result = self.converter.convert(pdf_path)
            markdown_content = result.document.export_to_markdown()
            
            # 保存 Markdown 文件
            with open(md_path, 'w', encoding='utf-8') as f:
                f.write(markdown_content)
            
            logger.info(f"PDF 轉換成功，Markdown 已保存至: {md_path}")
            return markdown_content
            
        except Exception as e:
            logger.error(f"PDF 轉換失敗: {e}")
            raise
    
    def extract_tables(self, markdown_content: str) -> List[str]:
        """
        從 Markdown 中提取所有表格
        
        Args:
            markdown_content: Markdown 內容
            
        Returns:
            表格列表
        """
        tables = []
        
        # 匹配 Markdown 表格（以 | 開頭的連續行）
        table_pattern = r'(?:^|\n)(\|[^\n]+\|(?:\n\|[^\n]+\|)+)'
        matches = re.findall(table_pattern, markdown_content, re.MULTILINE)
        
        for match in matches:
            # 清理並添加表格
            table = match.strip()
            if table and '|' in table:
                tables.append(table)
        
        # 也嘗試匹配 HTML 表格
        html_table_pattern = r'<table[^>]*>.*?</table>'
        html_matches = re.findall(html_table_pattern, markdown_content, re.DOTALL | re.IGNORECASE)
        tables.extend(html_matches)
        
        logger.info(f"提取到 {len(tables)} 個表格")
        return tables
    
    def extract_section(self, markdown_content: str, section_name: str) -> str:
        """
        提取特定的章節
        
        Args:
            markdown_content: Markdown 內容
            section_name: 章節名稱（如 'results', 'abstract' 等）
            
        Returns:
            該章節的內容
        """
        pattern = self.SECTION_PATTERNS.get(section_name.lower())
        if not pattern:
            return ""
        
        # 找到章節開始位置
        match = re.search(pattern, markdown_content, re.IGNORECASE)
        if not match:
            return ""
        
        start_pos = match.end()
        
        # 找到下一個章節的開始位置
        next_section_pos = len(markdown_content)
        for sec_name, sec_pattern in self.SECTION_PATTERNS.items():
            if sec_name != section_name.lower():
                next_match = re.search(sec_pattern, markdown_content[start_pos:], re.IGNORECASE)
                if next_match:
                    next_section_pos = min(next_section_pos, start_pos + next_match.start())
        
        section_content = markdown_content[start_pos:next_section_pos].strip()
        logger.debug(f"提取 {section_name} 章節: {len(section_content)} 字元")
        return section_content
    
    def extract_supplementary_links(self, markdown_content: str, pmid: Optional[str] = None) -> List[Dict[str, str]]:
        """
        提取 Supplementary Data 連結
        
        Args:
            markdown_content: Markdown 內容
            pmid: PubMed ID（可選，用於從網頁抓取）
            
        Returns:
            Supplementary 連結列表 [{"text": "...", "url": "..."}]
        """
        links = []
        
        # 方法 1: 從 Markdown 內容提取
        for pattern in self.SUPPLEMENTARY_LINK_PATTERNS:
            matches = re.findall(pattern, markdown_content, re.IGNORECASE)
            for match in matches:
                if isinstance(match, tuple):
                    text = match[0] if len(match) > 0 else ""
                    url = match[1] if len(match) > 1 else match[0]
                else:
                    text = "Supplementary Data"
                    url = match
                
                if url and url.startswith('http'):
                    links.append({
                        "name": text.strip() if text else "Supplementary Data",  # 使用 'name' 統一欄位名稱
                        "url": url.strip(),
                        "type": "unknown"  # 添加 type 欄位
                    })
        
        # 方法 2: 如果提供 PMID，從 PMC 網頁抓取（更準確）
        if pmid and len(links) == 0:
            try:
                from utils.supplementary_downloader import SupplementaryDownloader
                downloader = SupplementaryDownloader()
                web_links = downloader.scrape_pmc_supplementary_links(pmid)
                
                for web_link in web_links:
                    links.append({
                        "name": web_link['name'],  # 使用 'name' 而不是 'text'
                        "url": web_link['url'],
                        "type": web_link.get('type', 'unknown')
                    })
                
                if web_links:
                    logger.info(f"從 PMC 網頁抓取到 {len(web_links)} 個 Supplementary 連結")
                    
            except Exception as e:
                logger.debug(f"無法從網頁抓取 supplementary links: {e}")
        
        # 移除重複
        seen_urls = set()
        unique_links = []
        for link in links:
            if link['url'] not in seen_urls:
                seen_urls.add(link['url'])
                unique_links.append(link)
        
        logger.info(f"提取到 {len(unique_links)} 個 Supplementary 連結")
        return unique_links
    
    def extract_figures(self, markdown_content: str) -> List[str]:
        """
        提取圖片描述和連結
        
        Args:
            markdown_content: Markdown 內容
            
        Returns:
            圖片描述列表
        """
        figures = []
        
        # Markdown 圖片格式
        img_pattern = r'!\[([^\]]*)\]\(([^)]+)\)'
        matches = re.findall(img_pattern, markdown_content)
        
        for alt_text, url in matches:
            figures.append(f"![{alt_text}]({url})")
        
        # 也匹配圖片說明文字
        caption_pattern = r'(?:Figure|Fig\.?)\s*(\d+)[:\.]?\s*([^\n]+)'
        caption_matches = re.findall(caption_pattern, markdown_content, re.IGNORECASE)
        
        for fig_num, caption in caption_matches:
            figures.append(f"Figure {fig_num}: {caption.strip()}")
        
        logger.info(f"提取到 {len(figures)} 個圖片/圖說")
        return figures
    
    def create_priority_content(
        self,
        tables: List[str],
        results_section: str,
        supplementary_links: List[Dict[str, str]],
        abstract: str,
        full_markdown: str = "",
        supplementary_content: Optional[List[str]] = None,
        variant_aliases: Optional[List[str]] = None
    ) -> str:
        """
        創建優先處理的內容（根據用戶需求排序）
        
        優先順序：
        1. Supplementary Data (實際內容，如 Excel 表格)
        2. Tables (PDF 內的表格) - 如果提供 variant_aliases，只保留包含目標變異的表格
        3. Results 段落 (關鍵研究結果)
        4. 其他內容跳過（節省 token）
        
        Args:
            tables: 表格列表
            results_section: Results 段落
            supplementary_links: Supplementary 連結
            abstract: 摘要
            full_markdown: 完整 Markdown (暫不使用)
            supplementary_content: Supplementary files 的實際內容（如 Excel 表格轉 Markdown）
            variant_aliases: 目標變異的別名列表（用於過濾表格）
            
        Returns:
            優先內容的組合文本
        """
        priority_parts = []
        current_length = 0
        
        # 1. 最優先：Supplementary Data 實際內容
        if supplementary_content:
            supp_data_text = "## Supplementary Data (Extracted Tables)\n\n"
            # 根據是否使用過濾添加不同的註記
            if variant_aliases:
                supp_data_text += "*(Filtered: only tables from supplementary files mentioning the target variant)*\n\n"
            else:
                supp_data_text += "*(These tables are extracted from supplementary Excel/CSV files)*\n\n"
            
            for i, content in enumerate(supplementary_content, 1):
                entry = f"### Supplementary File {i}\n{content}\n\n"
                if current_length + len(entry) < self.max_priority_content_length:
                    supp_data_text += entry
                    current_length += len(entry)
                else:
                    supp_data_text += f"\n[Remaining supplementary files truncated to save space...]\n"
                    break
            
            if len(supp_data_text) > 100:  # 確保有內容
                priority_parts.append(supp_data_text)
        
        # 如果有 supplementary links 但沒有內容，提供連結
        if supplementary_links and not supplementary_content:
            supp_links_text = "## Supplementary Data Links\n\n"
            supp_links_text += "*(Important: These files may contain critical data tables)*\n\n"
            for link in supplementary_links:
                link_type = link.get('type', 'file')
                link_name = link.get('name', 'Supplementary File')
                supp_links_text += f"- **{link_name}** ({link_type}): {link['url']}\n"
            
            if current_length + len(supp_links_text) < self.max_priority_content_length:
                priority_parts.append(supp_links_text)
                current_length += len(supp_links_text)
        
        # 2. 第二優先：PDF 內的 Tables（使用 variant_aliases 過濾）
        if tables:
            # 如果提供了 variant_aliases，過濾只包含目標變異的表格
            if variant_aliases:
                filtered_tables = []
                for table in tables:
                    table_lower = table.lower()
                    # 檢查表格是否包含任何變異別名
                    if any(alias.lower() in table_lower for alias in variant_aliases if alias):
                        filtered_tables.append(table)
                
                if filtered_tables:
                    logger.info(f"表格過濾: {len(tables)} -> {len(filtered_tables)} 個包含目標變異的表格")
                    tables_to_use = filtered_tables
                else:
                    logger.warning(f"表格過濾後沒有匹配的表格，保留所有 {len(tables)} 個表格")
                    tables_to_use = tables
            else:
                tables_to_use = tables
            
            tables_text = "## Tables from Main Article\n\n"
            if variant_aliases:
                tables_text += "*(Filtered: only tables mentioning the target variant)*\n\n"
            
            for i, table in enumerate(tables_to_use, 1):
                table_entry = f"### Table {i}\n{table}\n\n"
                if current_length + len(table_entry) < self.max_priority_content_length:
                    tables_text += table_entry
                    current_length += len(table_entry)
                else:
                    tables_text += f"\n[Remaining tables truncated...]\n"
                    break
            
            if len(tables_text) > 50:
                priority_parts.append(tables_text)
        
        # 3. 第三優先：Results 段落
        if results_section:
            remaining = self.max_priority_content_length - current_length
            if remaining > 500:  # 至少保留 500 字元給 results
                results_text = "## Results Section\n\n"
                truncated_results = results_section[:remaining - 100]
                if len(results_section) > remaining - 100:
                    truncated_results += "\n\n[Results section truncated...]"
                results_text += truncated_results
                priority_parts.append(results_text)
                current_length += len(results_text)
        
        # 4. 其他內容一律跳過（abstract, methods, discussion 等）
        # 這樣可以節省大量 token，讓 LLM 專注在數據上
        
        if not priority_parts:
            # 如果完全沒有優先內容，提供一個最小的 fallback
            return f"## Note\n\nNo tables, supplementary data, or results section found in this article.\n\nAbstract: {abstract[:500] if abstract else 'N/A'}"
        
        priority_content = "\n\n".join(priority_parts)
        logger.info(f"創建優先內容: {len(priority_content)} 字元 (優先級: Supplementary Data > Tables > Results，其他內容已跳過)")
        return priority_content
    
    def process_pdf(self, pdf_path: str, pmid: Optional[str] = None, download_supplementary: bool = True, variant_aliases: Optional[List[str]] = None) -> ExtractedContent:
        """
        完整處理 PDF 文件 (整合 Playwright 下載)
        
        Args:
            pdf_path: PDF 文件路徑
            pmid: PubMed ID
            download_supplementary: 是否下載 supplementary files
            variant_aliases: 目標變異的別名列表（用於過濾表格）
        """
        logger.info(f"開始處理 PDF: {pdf_path}")
        
        # 轉換為 Markdown
        markdown_content = self.convert_pdf_to_markdown(pdf_path)
        tables = self.extract_tables(markdown_content)
        results_section = self.extract_section(markdown_content, 'results')
        figures = self.extract_figures(markdown_content)
        abstract = self.extract_section(markdown_content, 'abstract')
        methods_section = self.extract_section(markdown_content, 'methods')
        discussion_section = self.extract_section(markdown_content, 'discussion')
        
        # 嘗試從 Markdown 提取連結 (僅作參考)
        supplementary_links = [] # 您的 Playwright 腳本會自己找連結，不需要從 PDF 這裡抓

        # *** 修改部分：使用簡化版下載器 ***
        supplementary_content = []
        if download_supplementary and pmid:
            try:
                from utils.supplementary_downloader import SupplementaryDownloader
                downloader = SupplementaryDownloader(output_dir=self.output_dir.parent / "supplementary")
                
                # 直接執行您的下載邏輯
                downloaded_files = downloader.download_all_supplementary_files(pmid)
                
                # 簡單讀取下載下來的 Excel 內容（使用 variant_aliases 過濾）
                for file_info in downloaded_files:
                    if file_info['type'] == 'excel':
                        # 傳遞 variant_aliases 進行表格過濾
                        excel_tables = downloader.extract_tables_from_excel(
                            Path(file_info['local_path']), 
                            variant_aliases=variant_aliases
                        )
                        if excel_tables:
                            combined = "\n\n".join(excel_tables)
                            supplementary_content.append(f"**Supplementary File: {file_info['name']}**\n\n{combined}")
                            
            except Exception as e:
                logger.warning(f"Supplementary download failed: {e}")
        # *** 修改部分結束 ***

        # 創建優先內容 (這會被送到 LLM)
        priority_content = self.create_priority_content(
            tables, results_section, supplementary_links, abstract,
            full_markdown=markdown_content,
            supplementary_content=supplementary_content,
            variant_aliases=variant_aliases
        )
        
        return ExtractedContent(
            full_markdown=markdown_content,
            tables=tables,
            results_section=results_section,
            supplementary_links=supplementary_links,
            figures=[],
            abstract=abstract,
            methods_section="",
            discussion_section="",
            priority_content=priority_content
        )
    
    def process_pdf_for_llm(
        self,
        pdf_path: str,
        include_full_text: bool = False,
        pmid: Optional[str] = None,
        download_supplementary: bool = True,
        variant_aliases: Optional[List[str]] = None
    ) -> Dict[str, any]:
        """
        為 LLM 處理 PDF，返回適合模型輸入的格式
        
        Args:
            pdf_path: PDF 文件路徑
            include_full_text: 是否包含完整文本（預設 False 以節省 GPU 記憶體）
            pmid: PubMed ID（可選，用於從網頁抓取 supplementary links）
            download_supplementary: 是否下載並解析 supplementary files（預設 True）
            variant_aliases: 目標變異的別名列表（用於過濾表格）
            
        Returns:
            適合 LLM 輸入的字典
        """
        extracted = self.process_pdf(pdf_path, pmid=pmid, download_supplementary=download_supplementary, variant_aliases=variant_aliases)
        
        result = {
            "source": pdf_path,
            "priority_content": extracted.priority_content,
            "tables_count": len(extracted.tables),
            "supplementary_links": extracted.supplementary_links,
            "has_results": bool(extracted.results_section),
            "has_abstract": bool(extracted.abstract),
        }
        
        if include_full_text:
            result["full_markdown"] = extracted.full_markdown
        
        return result


def download_pdf_from_pmc(pmid: str, output_dir: Optional[Path] = None) -> Optional[str]:
    """
    從 PubMed Central 下載 PDF（如果可用）
    
    Args:
        pmid: PubMed ID
        output_dir: 輸出目錄
        
    Returns:
        PDF 文件路徑或 None
    """
    import requests
    try:
        from Bio import Entrez
    except ImportError:
        logger.warning("Biopython not installed, cannot download PDF from PMC")
        return None
    
    # 設定 Entrez email
    try:
        from config import PUBMED_EMAIL
        Entrez.email = PUBMED_EMAIL
    except ImportError:
        Entrez.email = "user@example.com"
    
    output_dir = output_dir or Path("./outputs/pdfs")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # 先查詢 PMC ID
        handle = Entrez.elink(
            dbfrom="pubmed",
            db="pmc",
            id=pmid,
            linkname="pubmed_pmc"
        )
        result = Entrez.read(handle)
        handle.close()
        
        if not result[0]['LinkSetDb']:
            logger.debug(f"PMID {pmid} 不在 PMC 中")
            return None
        
        pmc_id = result[0]['LinkSetDb'][0]['Link'][0]['Id']
        
        # 檢查 PDF 是否已存在
        pdf_path = output_dir / f"PMC{pmc_id}.pdf"
        if pdf_path.exists():
            logger.info(f"PDF 已存在: {pdf_path}")
            return str(pdf_path)
        
        # 設定 headers 模擬瀏覽器
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept': 'application/pdf,*/*',
        }
        
        # 嘗試多種 PDF 下載 URL 格式
        pdf_urls = [
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/main.pdf",
            f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/pdf/",
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/",
        ]
        
        for pdf_url in pdf_urls:
            try:
                response = requests.get(pdf_url, headers=headers, timeout=30, allow_redirects=True)
                content_type = response.headers.get('content-type', '').lower()
                
                if response.status_code == 200 and ('pdf' in content_type or response.content[:4] == b'%PDF'):
                    pdf_path = output_dir / f"PMC{pmc_id}.pdf"
                    with open(pdf_path, 'wb') as f:
                        f.write(response.content)
                    logger.info(f"PDF 下載成功: {pdf_path}")
                    return str(pdf_path)
            except Exception as url_error:
                logger.debug(f"URL {pdf_url} 失敗: {url_error}")
                continue
        
        # 如果直接下載失敗，嘗試使用 Europe PMC
        try:
            eupmc_url = f"https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC{pmc_id}&blobtype=pdf"
            response = requests.get(eupmc_url, headers=headers, timeout=30)
            
            if response.status_code == 200 and (response.content[:4] == b'%PDF'):
                pdf_path = output_dir / f"PMC{pmc_id}.pdf"
                with open(pdf_path, 'wb') as f:
                    f.write(response.content)
                logger.info(f"PDF 從 Europe PMC 下載成功: {pdf_path}")
                return str(pdf_path)
        except Exception as eupmc_error:
            logger.debug(f"Europe PMC 下載失敗: {eupmc_error}")
        
        logger.debug(f"PMID {pmid} (PMC{pmc_id}): 所有 PDF 下載方式都失敗")
        return None
        
    except Exception as e:
        logger.debug(f"下載 PDF 失敗: {e}")
        return None


# 便捷函數
def process_pubmed_pdf(
    pdf_path: str,
    max_priority_length: int = 8000
) -> Dict[str, any]:
    """
    便捷函數：處理 PubMed PDF 文件
    
    Args:
        pdf_path: PDF 文件路徑
        max_priority_length: 優先內容最大長度
        
    Returns:
        處理結果字典
    """
    processor = DoclingPDFProcessor(
        max_priority_content_length=max_priority_length,
        enable_ocr=True,
        enable_table_extraction=True
    )
    return processor.process_pdf_for_llm(pdf_path, include_full_text=False)


if __name__ == "__main__":
    # 測試用
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    if len(sys.argv) > 1:
        pdf_path = sys.argv[1]
        result = process_pubmed_pdf(pdf_path)
        print("\n=== 處理結果 ===")
        print(f"表格數量: {result['tables_count']}")
        print(f"Supplementary 連結: {len(result['supplementary_links'])}")
        print(f"有 Results 段落: {result['has_results']}")
        print(f"有 Abstract: {result['has_abstract']}")
        print(f"\n優先內容長度: {len(result['priority_content'])} 字元")
        print("\n=== 優先內容預覽 ===")
        print(result['priority_content'][:2000])
    else:
        print("用法: python docling_pdf_processor.py <pdf_path>")

