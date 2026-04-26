#!/usr/bin/env python3
"""
RAG Retriever for TTN Variant Pipeline
========================================

目的：讓 per-article 提取更準確
- 把單篇文章切成 chunk，用 FAISS 找出含有 variant 資訊的段落
- 解決「全文太長被截斷，剛好截掉有 variant 的表格」問題

依賴（輕量，CPU 就夠）：
    pip install sentence-transformers faiss-cpu
"""

import logging
import re
from typing import Dict, List

import numpy as np

logger = logging.getLogger(__name__)

_embedder = None


def _get_embedder(model_name: str = "pritamdeka/S-PubMedBert-MS-MARCO"):
    """
    懶載入 embedding 模型。
    - 預設：PubMedBERT，生醫文獻效果最好
    - 備用輕量版：'sentence-transformers/all-MiniLM-L6-v2'（更快，精度稍低）
    """
    global _embedder
    if _embedder is None:
        from sentence_transformers import SentenceTransformer
        logger.info(f"載入 embedding 模型: {model_name}")
        _embedder = SentenceTransformer(model_name)
        logger.info("Embedding 模型載入完成")
    return _embedder


def chunk_text(
    text: str,
    chunk_size: int = 400,
    overlap: int = 80,
    min_chunk_len: int = 50,
) -> List[str]:
    """
    把文章切成有重疊的段落，優先在段落/句子邊界切割。
    chunk_size: 約幾個 token（1 token ≈ 4 英文字元）
    """
    char_size = chunk_size * 4
    char_overlap = overlap * 4

    paragraphs = re.split(r'\n{2,}', text)
    chunks: List[str] = []
    current = ""

    for para in paragraphs:
        para = para.strip()
        if not para:
            continue

        if len(para) > char_size:
            for sent in re.split(r'(?<=[.!?])\s+', para):
                if len(current) + len(sent) <= char_size:
                    current += (" " + sent) if current else sent
                else:
                    if len(current) >= min_chunk_len:
                        chunks.append(current.strip())
                    overlap_text = current[-char_overlap:] if len(current) > char_overlap else current
                    current = overlap_text + " " + sent
        else:
            if len(current) + len(para) <= char_size:
                current += ("\n\n" + para) if current else para
            else:
                if len(current) >= min_chunk_len:
                    chunks.append(current.strip())
                overlap_text = current[-char_overlap:] if len(current) > char_overlap else current
                current = overlap_text + "\n\n" + para

    if len(current) >= min_chunk_len:
        chunks.append(current.strip())

    return chunks


class ArticleRAG:
    """
    對單篇文章建立小型 FAISS index，根據查詢撈最相關段落。
    取代截斷邏輯（之前只取文章頭部，容易錯過文章後段的 variant 表格）。

    用法：
        rag = ArticleRAG()
        context = rag.retrieve(article, queries=aliases + ["disease phenotype"])
    """

    def __init__(
        self,
        embedder_model: str = "pritamdeka/S-PubMedBert-MS-MARCO",
        top_k: int = 8,
        chunk_size: int = 400,
        overlap: int = 80,
    ):
        self.embedder_model = embedder_model
        self.top_k = top_k
        self.chunk_size = chunk_size
        self.overlap = overlap

    def retrieve(
        self,
        article: Dict,
        queries: List[str],
        max_total_chars: int = 12_000,
    ) -> str:
        """
        從單篇文章中撈最相關段落，組成帶標籤的 context 字串。

        Args:
            article:         pubmed article dict（含 title, abstract, full_text 等欄位）
            queries:         查詢列表（variant aliases + 疾病關鍵字）
            max_total_chars: 回傳字串的字元上限

        Returns:
            結構化 context 字串，可直接放進 LLM prompt
        """
        title = (article.get("title") or "").strip()
        abstract = (article.get("abstract") or "").strip()
        full_text = (
            article.get("text_for_llm")
            or article.get("full_text")
            or ""
        ).strip()

        # 永遠保留 title + abstract
        header_parts = []
        if title:
            header_parts.append(f"[TITLE]\n{title}")
        if abstract:
            header_parts.append(f"[ABSTRACT]\n{abstract[:2_000]}")
        header = "\n\n".join(header_parts)

        remaining_budget = max_total_chars - len(header) - 100
        if not full_text or remaining_budget < 300:
            return header

        try:
            import faiss

            chunks = chunk_text(full_text, self.chunk_size, self.overlap)
            if not chunks:
                return header

            embedder = _get_embedder(self.embedder_model)

            chunk_embs = embedder.encode(chunks, show_progress_bar=False, batch_size=32)
            chunk_embs = np.array(chunk_embs, dtype=np.float32)
            faiss.normalize_L2(chunk_embs)

            index = faiss.IndexFlatIP(chunk_embs.shape[1])
            index.add(chunk_embs)

            query_embs = embedder.encode(queries, show_progress_bar=False)
            query_embs = np.array(query_embs, dtype=np.float32)
            faiss.normalize_L2(query_embs)

            k = min(self.top_k, len(chunks))
            _, indices = index.search(query_embs, k)

            # 聯集，按原文順序排（保持閱讀連貫性）
            seen = set()
            hit_indices = []
            for row in indices:
                for idx in row:
                    if idx >= 0 and idx not in seen:
                        hit_indices.append(idx)
                        seen.add(idx)
            hit_indices.sort()

            parts = []
            used = 0
            for idx in hit_indices:
                chunk = chunks[idx].strip()
                if used + len(chunk) > remaining_budget:
                    leftover = remaining_budget - used
                    if leftover > 200:
                        parts.append(chunk[:leftover] + "...")
                    break
                parts.append(chunk)
                used += len(chunk)

            if parts:
                body = "\n\n---\n\n".join(parts)
                return f"{header}\n\n[RELEVANT SECTIONS (RAG-retrieved from full text)]\n{body}"

            return header

        except Exception as e:
            logger.warning(f"ArticleRAG 失敗，使用截取 fallback: {e}")
            return f"{header}\n\n[FULL TEXT (truncated)]\n{full_text[:remaining_budget]}"
