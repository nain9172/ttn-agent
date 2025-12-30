#!/usr/bin/env python3
"""
診斷 Supplementary Files 抓取和下載
"""

import logging
import sys
from utils.supplementary_downloader import SupplementaryDownloader

# 設定詳細的日誌
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def debug_pmid(pmid: str, test_download: bool = False):
    """診斷某個 PMID 的 supplementary files 抓取"""
    
    print(f"\n{'='*80}")
    print(f"診斷 PMID {pmid} 的 Supplementary Files 抓取")
    print(f"{'='*80}\n")
    
    downloader = SupplementaryDownloader()
    
    # 步驟 1: 獲取 PMC ID
    print("步驟 1: 獲取 PMC ID...")
    pmc_id = downloader.get_pmc_id_from_pmid(pmid)
    if pmc_id:
        print(f"  ✓ PMC ID: PMC{pmc_id}")
    else:
        print(f"  ✗ 無法找到 PMC ID（文章可能不在 PMC 中）")
        return
    
    # 步驟 2: 使用 Entrez API 抓取連結
    print(f"\n步驟 2: 使用 Entrez API 獲取 supplementary files...")
    entrez_links = downloader.get_supplementary_from_entrez(pmid)
    
    if entrez_links:
        print(f"  ✓ Entrez API 找到 {len(entrez_links)} 個文件")
        links = entrez_links
    else:
        print(f"  ⚠ Entrez API 未找到，嘗試網頁抓取...")
        print(f"  URL: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/")
        links = downloader.scrape_pmc_supplementary_links(pmid)
    print(f"\n  ✓ 找到 {len(links)} 個 supplementary files\n")
    
    if not links:
        print("  ⚠ 未找到任何 supplementary files")
        return
    
    # 顯示連結詳情
    for i, link in enumerate(links, 1):
        print(f"  [{i}] {link['name']}")
        print(f"      類型: {link['type']}")
        print(f"      URL: {link['url'][:100]}{'...' if len(link['url']) > 100 else ''}")
        print()
    
    # 步驟 3: 測試下載（可選）
    if test_download and links:
        print(f"\n步驟 3: 測試下載前 {min(2, len(links))} 個文件...")
        
        for i, link in enumerate(links[:2], 1):
            print(f"\n  測試 [{i}]: {link['name']}")
            print(f"    URL: {link['url']}")
            
            local_path = downloader.download_supplementary_file(
                link['url'],
                pmid=pmid
            )
            
            if local_path:
                import os
                file_size = os.path.getsize(local_path)
                print(f"    ✓ 下載成功: {local_path}")
                print(f"    ✓ 文件大小: {file_size:,} bytes")
            else:
                print(f"    ✗ 下載失敗")
    
    print(f"\n{'='*80}")
    print("診斷完成")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python debug_supplementary.py <PMID> [--download]")
        print("\n範例:")
        print("  python debug_supplementary.py 26147798")
        print("  python debug_supplementary.py 26147798 --download")
        sys.exit(1)
    
    pmid = sys.argv[1]
    test_download = '--download' in sys.argv
    
    debug_pmid(pmid, test_download)

