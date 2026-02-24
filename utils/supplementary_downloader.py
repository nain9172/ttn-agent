#!/usr/bin/env python3
"""
Supplementary Files Downloader (Simple Playwright Version)
"""

import logging
import os
import time
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
from Bio import Entrez
from playwright.sync_api import sync_playwright

logger = logging.getLogger(__name__)

# 設定 Entrez email (僅用於 PMID -> PMC ID 轉換)
try:
    from config import PUBMED_EMAIL, ALLOWED_SUPPLEMENTARY_FORMATS, FORCE_REDOWNLOAD_SUPPLEMENTARY
    Entrez.email = PUBMED_EMAIL
except ImportError:
    Entrez.email = "user@example.com"
    ALLOWED_SUPPLEMENTARY_FORMATS = ['.doc', '.docx', '.xlsx', '.xls', '.csv', '.pdf']
    FORCE_REDOWNLOAD_SUPPLEMENTARY = False

class SupplementaryDownloader:
    def __init__(self, output_dir: Optional[Path] = None):
        self.output_dir = output_dir or Path("./outputs/supplementary")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def get_pmc_id_from_pmid(self, pmid: str) -> Optional[str]:
        """簡單轉換：從 PMID 查 PMC ID (這是進入網址的必要步驟)"""
        try:
            handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc")
            result = Entrez.read(handle)
            handle.close()
            if not result[0]['LinkSetDb']: return None
            return result[0]['LinkSetDb'][0]['Link'][0]['Id'] # 回傳純數字 ID
        except:
            return None

    def download_all_supplementary_files(self, pmid: str) -> List[Dict[str, any]]:
        """
        使用 Playwright 下載 (邏輯與您提供的一致)
        如果文件已經下載過，則直接返回已存在的文件信息
        """
        pmc_id_num = self.get_pmc_id_from_pmid(pmid)
        if not pmc_id_num:
            logger.warning(f"PMID {pmid} 無對應 PMC ID，跳過。")
            return []

        PMC_ID = f"PMC{pmc_id_num}"
        BASE = "https://pmc.ncbi.nlm.nih.gov"
        
        # 設定輸出目錄
        out_dir = self.output_dir / f"pmid_{pmid}"
        
        # 檢查目錄是否已存在且有文件（除非強制重新下載）
        if out_dir.exists() and not FORCE_REDOWNLOAD_SUPPLEMENTARY:
            existing_files = list(out_dir.glob('*'))
            # 過濾掉目錄，只保留文件
            existing_files = [f for f in existing_files if f.is_file()]
            
            if existing_files:
                logger.info(f"PMID {pmid}: 發現 {len(existing_files)} 個已下載的 supplementary files，跳過下載")
                results = []
                
                # 從已存在的文件創建結果列表
                for file_path in existing_files:
                    file_ext = file_path.suffix.lower()
                    
                    # 只包含允許的格式
                    if file_ext in ALLOWED_SUPPLEMENTARY_FORMATS:
                        # 根據副檔名判斷檔案類型
                        if file_ext in ['.xlsx', '.xls', '.csv']:
                            file_type = 'excel'
                        elif file_ext in ['.doc', '.docx']:
                            file_type = 'word'
                        elif file_ext == '.pdf':
                            file_type = 'pdf'
                        else:
                            file_type = 'other'
                        
                        results.append({
                            'name': file_path.name,
                            'local_path': str(file_path),
                            'type': file_type
                        })
                        logger.info(f"  - 使用已存在文件: {file_path.name}")
                
                if results:
                    return results
                else:
                    logger.info(f"  已存在的文件都不在允許格式列表中，重新下載")
        
        # 如果目錄不存在或沒有有效文件，創建目錄並進行下載
        os.makedirs(out_dir, exist_ok=True)
        
        results = []
        
        logger.info(f"啟動 Playwright 下載: {PMC_ID}")

        try:
            with sync_playwright() as p:
                browser = p.chromium.launch(headless=True)
                context = browser.new_context(accept_downloads=True)
                page = context.new_page()
                
                # 開文章頁面
                url = f"{BASE}/articles/{PMC_ID}/"
                logger.info(f"前往: {url}")
                page.goto(url, wait_until="networkidle")
                
                # 找 supplement 的連結 (使用您指定的 locator)
                links = page.locator('a[href*="/articles/instance/"]').all()
                logger.info(f"Found {len(links)} supplements")
                # print(f'links: {links}')
                
                for i, link in enumerate(links, 1):
                    try:
                        href = link.get_attribute("href")
                        filename = href.split("/")[-1]
                        
                        # 檢查檔案格式是否在允許列表中
                        file_ext = Path(filename).suffix.lower()
                        if file_ext not in ALLOWED_SUPPLEMENTARY_FORMATS:
                            logger.info(f"[{i}/{len(links)}] Skipping {filename} (format not allowed: {file_ext})")
                            continue
                        
                        # 處理重複檔名
                        save_path = out_dir / filename
                        
                        logger.info(f"[{i}/{len(links)}] Downloading: {filename}")
                        
                        # 點擊並等待下載
                        with page.expect_download(timeout=60000) as d:
                            # link.scroll_into_view_if_needed() # 確保看得到才能點
                            link.click()
                        
                        download = d.value
                        print(f'download: {download}')
                        download.save_as(str(save_path))
                        
                        # 記錄成功下載的檔案資訊 (給後續程式讀取用)
                        # 根據副檔名判斷檔案類型
                        if file_ext in ['.xlsx', '.xls', '.csv']:
                            file_type = 'excel'
                        elif file_ext in ['.doc', '.docx']:
                            file_type = 'word'
                        elif file_ext == '.pdf':
                            file_type = 'pdf'
                        else:
                            file_type = 'other'
                        
                        results.append({
                            'name': filename,
                            'local_path': str(save_path),
                            'type': file_type
                        })
                        
                        time.sleep(1)
                        
                    except Exception as e:
                        logger.error(f"Failed to download link {i}: {e}")
                        
                browser.close()
                
        except Exception as e:
            logger.error(f"Playwright Error: {e}")

        return results

    def extract_tables_from_excel(self, excel_path: Path, variant_aliases: Optional[List[str]] = None) -> List[str]:
        """
        簡單讀取 Excel 轉 Markdown (給 LLM 看內容用)
        
        Args:
            excel_path: Excel 文件路徑
            variant_aliases: 變異別名列表，用於過濾只包含目標變異的表格
        
        Returns:
            表格列表（Markdown 格式）
        """
        try:
            tables = []
            file_str = str(excel_path).lower()
            
            if file_str.endswith('.csv'):
                df = pd.read_csv(excel_path)
                table_md = f"### {excel_path.name}\n{df.to_markdown(index=False)}"
                
                # 如果提供了 variant_aliases，檢查表格是否包含目標變異
                if variant_aliases:
                    table_lower = table_md.lower()
                    if any(alias.lower() in table_lower for alias in variant_aliases if alias):
                        tables.append(table_md)
                    else:
                        logger.debug(f"Supplementary 表格過濾: {excel_path.name} 不包含目標變異，已跳過")
                else:
                    tables.append(table_md)
            else:
                xl = pd.ExcelFile(excel_path)
                for sheet in xl.sheet_names:
                    df = pd.read_excel(xl, sheet_name=sheet)
                    table_md = f"### {excel_path.name} - {sheet}\n{df.to_markdown(index=False)}"
                    
                    # 如果提供了 variant_aliases，檢查表格是否包含目標變異
                    if variant_aliases:
                        table_lower = table_md.lower()
                        if any(alias.lower() in table_lower for alias in variant_aliases if alias):
                            tables.append(table_md)
                        else:
                            logger.debug(f"Supplementary 表格過濾: {excel_path.name} - {sheet} 不包含目標變異，已跳過")
                    else:
                        tables.append(table_md)
            
            if variant_aliases and tables:
                logger.info(f"Supplementary 表格過濾: {excel_path.name} 保留 {len(tables)} 個包含目標變異的表格/工作表")
            
            return tables
        except Exception as e:
            logger.error(f"Error extracting tables from {excel_path}: {e}")
            return []