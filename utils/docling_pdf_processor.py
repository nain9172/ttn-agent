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
        logger.info(f"開始轉換 PDF: {pdf_path}")
        
        try:
            result = self.converter.convert(pdf_path)
            markdown_content = result.document.export_to_markdown()
            
            # 保存 Markdown 文件
            pdf_name = Path(pdf_path).stem
            md_path = self.output_dir / f"{pdf_name}.md"
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
    
    def extract_supplementary_links(self, markdown_content: str) -> List[Dict[str, str]]:
        """
        提取 Supplementary Data 連結
        
        Args:
            markdown_content: Markdown 內容
            
        Returns:
            Supplementary 連結列表 [{"text": "...", "url": "..."}]
        """
        links = []
        
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
                        "text": text.strip() if text else "Supplementary Data",
                        "url": url.strip()
                    })
        
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
        abstract: str
    ) -> str:
        """
        創建優先處理的內容（用於有限 GPU 記憶體）
        
        優先順序：
        1. Tables（最重要的數據）
        2. Results 段落
        3. Supplementary Data 連結
        4. Abstract（如果還有空間）
        
        Args:
            tables: 表格列表
            results_section: Results 段落
            supplementary_links: Supplementary 連結
            abstract: 摘要
            
        Returns:
            優先內容的組合文本
        """
        priority_parts = []
        current_length = 0
        
        # 1. 添加表格
        if tables:
            tables_text = "## Tables\n\n"
            for i, table in enumerate(tables, 1):
                table_entry = f"### Table {i}\n{table}\n\n"
                if current_length + len(table_entry) < self.max_priority_content_length:
                    tables_text += table_entry
                    current_length += len(table_entry)
            priority_parts.append(tables_text)
        
        # 2. 添加 Results 段落
        if results_section:
            remaining = self.max_priority_content_length - current_length
            if remaining > 500:  # 至少保留 500 字元給 results
                results_text = "## Results\n\n"
                truncated_results = results_section[:remaining - 100]
                if len(results_section) > remaining - 100:
                    truncated_results += "\n\n[Results section truncated...]"
                results_text += truncated_results
                priority_parts.append(results_text)
                current_length += len(results_text)
        
        # 3. 添加 Supplementary 連結
        if supplementary_links:
            supp_text = "## Supplementary Data Links\n\n"
            for link in supplementary_links:
                supp_text += f"- [{link['text']}]({link['url']})\n"
            if current_length + len(supp_text) < self.max_priority_content_length:
                priority_parts.append(supp_text)
                current_length += len(supp_text)
        
        # 4. 添加 Abstract（如果還有空間）
        if abstract:
            remaining = self.max_priority_content_length - current_length
            if remaining > 300:
                abstract_text = "## Abstract\n\n"
                truncated_abstract = abstract[:remaining - 50]
                if len(abstract) > remaining - 50:
                    truncated_abstract += "..."
                abstract_text += truncated_abstract
                priority_parts.append(abstract_text)
        
        priority_content = "\n\n".join(priority_parts)
        logger.info(f"創建優先內容: {len(priority_content)} 字元")
        return priority_content
    
    def process_pdf(self, pdf_path: str) -> ExtractedContent:
        """
        完整處理 PDF 文件
        
        Args:
            pdf_path: PDF 文件路徑
            
        Returns:
            ExtractedContent 包含所有提取的內容
        """
        logger.info(f"開始處理 PDF: {pdf_path}")
        
        # 轉換為 Markdown
        markdown_content = self.convert_pdf_to_markdown(pdf_path)
        
        # 提取各部分
        tables = self.extract_tables(markdown_content)
        results_section = self.extract_section(markdown_content, 'results')
        supplementary_links = self.extract_supplementary_links(markdown_content)
        figures = self.extract_figures(markdown_content)
        abstract = self.extract_section(markdown_content, 'abstract')
        methods_section = self.extract_section(markdown_content, 'methods')
        discussion_section = self.extract_section(markdown_content, 'discussion')
        
        # 創建優先內容
        priority_content = self.create_priority_content(
            tables, results_section, supplementary_links, abstract
        )
        
        return ExtractedContent(
            full_markdown=markdown_content,
            tables=tables,
            results_section=results_section,
            supplementary_links=supplementary_links,
            figures=figures,
            abstract=abstract,
            methods_section=methods_section,
            discussion_section=discussion_section,
            priority_content=priority_content
        )
    
    def process_pdf_for_llm(
        self,
        pdf_path: str,
        include_full_text: bool = False
    ) -> Dict[str, any]:
        """
        為 LLM 處理 PDF，返回適合模型輸入的格式
        
        Args:
            pdf_path: PDF 文件路徑
            include_full_text: 是否包含完整文本（預設 False 以節省 GPU 記憶體）
            
        Returns:
            適合 LLM 輸入的字典
        """
        extracted = self.process_pdf(pdf_path)
        
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
    from Bio import Entrez
    
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
        
        # 嘗試下載 PDF
        pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/"
        response = requests.get(pdf_url, timeout=30, allow_redirects=True)
        
        if response.status_code == 200 and 'pdf' in response.headers.get('content-type', '').lower():
            pdf_path = output_dir / f"PMC{pmc_id}.pdf"
            with open(pdf_path, 'wb') as f:
                f.write(response.content)
            logger.info(f"PDF 下載成功: {pdf_path}")
            return str(pdf_path)
        
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

