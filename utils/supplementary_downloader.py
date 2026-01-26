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
    from config import PUBMED_EMAIL
    Entrez.email = PUBMED_EMAIL
except ImportError:
    Entrez.email = "user@example.com"

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
        """
        pmc_id_num = self.get_pmc_id_from_pmid(pmid)
        if not pmc_id_num:
            logger.warning(f"PMID {pmid} 無對應 PMC ID，跳過。")
            return []

        PMC_ID = f"PMC{pmc_id_num}"
        BASE = "https://pmc.ncbi.nlm.nih.gov"
        
        # 設定輸出目錄
        out_dir = self.output_dir / f"pmid_{pmid}"
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
                        results.append({
                            'name': filename,
                            'local_path': str(save_path),
                            'type': 'excel' if 'xls' in filename or 'csv' in filename else 'other'
                        })
                        
                        time.sleep(1)
                        
                    except Exception as e:
                        logger.error(f"Failed to download link {i}: {e}")
                        
                browser.close()
                
        except Exception as e:
            logger.error(f"Playwright Error: {e}")

        return results

    def extract_tables_from_excel(self, excel_path: Path) -> List[str]:
        """簡單讀取 Excel 轉 Markdown (給 LLM 看內容用)"""
        try:
            tables = []
            file_str = str(excel_path).lower()
            
            if file_str.endswith('.csv'):
                df = pd.read_csv(excel_path)
                tables.append(f"### {excel_path.name}\n{df.to_markdown(index=False)}")
            else:
                xl = pd.ExcelFile(excel_path)
                for sheet in xl.sheet_names:
                    df = pd.read_excel(xl, sheet_name=sheet)
                    tables.append(f"### {excel_path.name} - {sheet}\n{df.to_markdown(index=False)}")
            return tables
        except:
            return []