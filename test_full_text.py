#!/usr/bin/env python3
"""
測試腳本：驗證全文獲取和 LLM 分析功能
"""

import logging
from utils.enhanced_pubmed_search import EnhancedPubMedSearcher

# 設置日誌
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """測試全文獲取功能"""
    
    # 測試 PMID（使用者提供的範例）
    test_pmid = '34135346'
    
    logger.info("=" * 70)
    logger.info("測試：增強版 PubMed 搜尋器 - 全文獲取功能")
    logger.info("=" * 70)
    
    # 初始化搜尋器
    searcher = EnhancedPubMedSearcher(try_full_text=True, max_text_length=8000)
    
    # 測試 variant info (dummy data)
    variant_info = {
        'chrom': '2',
        'pos': 178612477,
        'ref': 'T',
        'alt': 'A'
    }
    
    # 使用 PMID 列表搜尋
    logger.info(f"\n📚 正在獲取 PMID {test_pmid} 的文章...")
    articles = searcher.search(variant_info, pmid_list=[test_pmid])
    
    if not articles:
        logger.error("❌ 無法獲取文章")
        return
    
    article = articles[0]
    
    # 顯示結果
    logger.info("\n" + "=" * 70)
    logger.info("📊 文章資訊：")
    logger.info("=" * 70)
    logger.info(f"PMID: {article['pmid']}")
    logger.info(f"Title: {article['title']}")
    logger.info(f"DOI: {article.get('doi', 'N/A')}")
    logger.info(f"PubMed Link: {article['pubmed_link']}")
    if article.get('doi_link'):
        logger.info(f"DOI Link: {article['doi_link']}")
    
    logger.info("\n" + "=" * 70)
    logger.info("📄 文本資訊：")
    logger.info("=" * 70)
    logger.info(f"Has full text: {article.get('has_full_text', False)}")
    logger.info(f"Abstract length: {len(article.get('abstract', ''))} 字元")
    
    if article.get('full_text'):
        logger.info(f"Full text length: {len(article['full_text'])} 字元")
        logger.info(f"Text for LLM length: {len(article['text_for_llm'])} 字元")
        
        # 顯示前 500 字元
        logger.info("\n前 500 字元：")
        logger.info("-" * 70)
        logger.info(article['text_for_llm'][:500])
        logger.info("-" * 70)
    else:
        logger.info("⚠️  未能獲取全文，將使用摘要")
        logger.info(f"Abstract 前 300 字元：")
        logger.info("-" * 70)
        logger.info(article.get('abstract', '')[:300])
        logger.info("-" * 70)
    
    logger.info("\n" + "=" * 70)
    logger.info("✅ 測試完成")
    logger.info("=" * 70)
    
    # 提示下一步
    logger.info("\n💡 提示：")
    logger.info("1. 全文已成功獲取（如果顯示 has_full_text=True）")
    logger.info("2. LLM 分析時會使用 'text_for_llm' 欄位（優先使用全文）")
    logger.info("3. 修改後的 batch_extract 方法會自動使用全文")
    logger.info("4. max_model_len 已增加到 8192，可處理更長的文章")
    logger.info("5. 文本限制調整為 5000 字元（若超過會截斷）")

if __name__ == "__main__":
    main()

