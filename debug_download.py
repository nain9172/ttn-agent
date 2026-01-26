import logging
import sys
from pathlib import Path
from utils.docling_pdf_processor import DoclingPDFProcessor
from utils.playwright_downloader import PlaywrightDownloader
from config import OUTPUT_DIR

# 設定 Logging 到控制台
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def test_integration():
    # 測試用 PMID (這是你提供的範例文章 PMC7654353 的 PMID)
    # PMID: 33184654 -> PMC: PMC7654353
    test_pmid = "33170376" 
    
    print("="*50)
    print(f"測試目標: PMID {test_pmid}")
    print(f"預期輸出根目錄: {OUTPUT_DIR.absolute()}")
    
    # 1. 測試路徑計算
    processor = DoclingPDFProcessor(output_dir=OUTPUT_DIR / "markdown")
    expected_supp_dir = processor.output_dir.parent / "supplementary" / test_pmid
    print(f"預期下載資料夾: {expected_supp_dir.absolute()}")
    
    # 2. 測試 PMC ID 轉換
    pmc_id = processor._get_pmcid_from_pmid(test_pmid)
    print(f"PMC ID 轉換結果: {pmc_id}")
    
    if not pmc_id:
        print("❌ 錯誤: 無法獲取 PMC ID，下載中止")
        return

    # 3. 測試 Playwright 下載
    print("\n啟動 Playwright 下載測試...")
    try:
        with PlaywrightDownloader(headless=True) as downloader:
            # 強制印出實際訪問的 URL
            target_url = f"{downloader.base_url}/articles/{pmc_id}/"
            print(f"正在訪問: {target_url}")
            
            results = downloader.download_all_supplements(pmc_id, expected_supp_dir)
            
            print(f"\n下載結果報告 (共 {len(results)} 個檔案):")
            for res in results:
                status = "✅" if res['downloaded'] else "❌"
                print(f"{status} {res['name']}")
                print(f"   └─ 本地路徑: {res.get('local_path')}")
                
    except Exception as e:
        print(f"❌ Playwright 發生例外錯誤: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_integration()