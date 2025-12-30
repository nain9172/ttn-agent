#!/usr/bin/env python3
"""
測試 Supplementary Files 優先級讀取
驗證 LLM prompt 是否按照優先級排序：Supplementary > Tables > Results
"""

import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)

def test_supplementary_priority(pmid: str):
    """測試特定 PMID 的 supplementary 處理"""
    
    print(f"\n{'='*70}")
    print(f"測試 PMID {pmid} 的 Supplementary Files 優先級處理")
    print(f"{'='*70}\n")
    
    # Step 1: 檢查是否有 supplementary files
    print("Step 1: 檢查 Supplementary Files...")
    from utils.supplementary_downloader import get_supplementary_info_for_pmid
    
    info = get_supplementary_info_for_pmid(pmid)
    print(f"  ✓ 找到 {info['total_count']} 個 supplementary files")
    for i, link in enumerate(info['links'], 1):
        print(f"    {i}. {link['name']} ({link['type']})")
    
    # Step 2: 使用 Docling 處理（會自動下載並解析）
    print(f"\nStep 2: 使用 Docling 處理 PDF + Supplementary Files...")
    from utils.docling_pdf_processor import DoclingPDFProcessor, download_pdf_from_pmc
    from config import OUTPUT_DIR
    
    # 下載 PDF
    pdf_path = download_pdf_from_pmc(pmid, output_dir=OUTPUT_DIR / "pdfs")
    if not pdf_path:
        print(f"  ✗ 無法下載 PDF")
        return
    
    print(f"  ✓ PDF 路徑: {pdf_path}")
    
    # 處理 PDF（啟用 supplementary 下載）
    processor = DoclingPDFProcessor(
        max_priority_content_length=12000,  # 增加到 12000 以容納更多內容
        output_dir=OUTPUT_DIR / "markdown"
    )
    
    result = processor.process_pdf_for_llm(
        pdf_path,
        pmid=pmid,
        download_supplementary=True  # 啟用 supplementary 下載
    )
    
    # Step 3: 顯示優先內容
    print(f"\nStep 3: 優先內容摘要")
    print(f"  - Tables: {result['tables_count']} 個")
    print(f"  - Supplementary Links: {len(result['supplementary_links'])} 個")
    print(f"  - Has Results: {result['has_results']}")
    print(f"  - Priority Content Length: {len(result['priority_content'])} 字元")
    
    # Step 4: 顯示優先內容的結構
    print(f"\nStep 4: 優先內容結構預覽")
    print(f"{'='*70}")
    
    content = result['priority_content']
    
    # 顯示前 2000 字元
    preview = content[:2000]
    print(preview)
    print(f"\n... (showing first 2000 chars of {len(content)} total)")
    print(f"{'='*70}")
    
    # Step 5: 檢查是否包含 supplementary data
    has_supp_data = "Supplementary Data" in content or "Supplementary File" in content
    has_tables = "## Tables" in content
    has_results = "## Results" in content
    
    print(f"\nStep 5: 內容驗證")
    print(f"  ✓ 包含 Supplementary Data: {has_supp_data}")
    print(f"  ✓ 包含 Tables: {has_tables}")
    print(f"  ✓ 包含 Results: {has_results}")
    
    # 檢查順序
    if has_supp_data and has_tables and has_results:
        supp_pos = content.find("Supplementary")
        table_pos = content.find("## Tables")
        results_pos = content.find("## Results")
        
        print(f"\n  順序檢查:")
        print(f"    Supplementary 位置: {supp_pos}")
        print(f"    Tables 位置: {table_pos}")
        print(f"    Results 位置: {results_pos}")
        
        if supp_pos < table_pos < results_pos or (supp_pos < results_pos and table_pos < 0):
            print(f"  ✓ 順序正確：Supplementary -> Tables -> Results")
        else:
            print(f"  ⚠ 順序可能不符合預期")
    
    print(f"\n{'='*70}")
    print(f"測試完成！")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        pmid = sys.argv[1]
    else:
        # 預設測試 PMID（已知有 supplementary files）
        pmid = "29255378"
        print(f"未提供 PMID，使用預設值: {pmid}")
    
    test_supplementary_priority(pmid)

