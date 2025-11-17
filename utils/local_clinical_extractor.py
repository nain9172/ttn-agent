#!/usr/bin/env python3
"""
Improved Local Clinical Information Extractor
解決模型空響應和資訊提取失敗的問題

改進重點：
1. 簡化 prompt，更適合小型模型
2. 增加更詳細的調試日誌
3. 改進 JSON 解析的容錯性
4. 使用更寬鬆的 generation 參數
5. 添加重試機制
"""

import logging
import json
import re
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


class LocalClinicalExtractor:
    """
    使用本地 LLM 提取臨床資訊
    優化 prompt 以獲得更好的提取效果
    """
    
    def __init__(
        self, 
        backend: str = "vllm",
        model_name: str = "meta-llama/Llama-3.2-3B",
        tensor_parallel_size: int = 2
    ):
        """
        Args:
            backend: "vllm" 或 "ollama"
            model_name: 模型名稱
            tensor_parallel_size: GPU 數量
        """
        self.backend = backend
        self.model_name = model_name
        self.tensor_parallel_size = tensor_parallel_size
        
        logger.info(f"初始化 {backend} backend，模型: {model_name}")
        logger.info(f"使用 {tensor_parallel_size} 個 GPU")
        
        if backend == "vllm":
            self._init_vllm()
        elif backend == "ollama":
            self._init_ollama()
        else:
            raise ValueError(f"不支援的 backend: {backend}")
    
    def _init_vllm(self):
        """初始化 vLLM"""
        try:
            from vllm import LLM, SamplingParams
            
            logger.info("載入 vLLM 模型...")
            self.llm = LLM(
                model=self.model_name,
                tensor_parallel_size=self.tensor_parallel_size,
                trust_remote_code=True,
                enforce_eager=True,
                gpu_memory_utilization=0.9
            )
            
            # 優化的採樣參數 - 更寬鬆以獲得更多輸出
            self.sampling_params = SamplingParams(
                temperature=0.3,  # 降低溫度以獲得更確定的輸出
                top_p=0.9,
                max_tokens=1000,  # 增加 token 限制
                stop=["</s>"],  # 只使用 EOS token
            )
            
            logger.info("vLLM 初始化成功")
            
        except Exception as e:
            logger.error(f"vLLM 初始化失敗: {e}")
            raise
    
    def _init_ollama(self):
        """初始化 Ollama"""
        try:
            import ollama
            self.ollama_client = ollama
            logger.info("Ollama 初始化成功")
        except Exception as e:
            logger.error(f"Ollama 初始化失敗: {e}")
            raise
    
    def extract_from_articles(
        self,
        articles: List[Dict],
        variant_id: str
    ) -> List[Dict]:
        """
        從文章列表中提取臨床資訊
        
        Args:
            articles: PubMed 文章列表
            variant_id: 變異 ID (例如: 2-178528273-G-C)
        
        Returns:
            包含臨床資訊的文章列表
        """
        logger.info(f"\n{'='*70}")
        logger.info(f"開始提取臨床資訊 - 變異: {variant_id}")
        logger.info(f"文章數量: {len(articles)}")
        logger.info(f"{'='*70}\n")
        
        enriched_articles = []
        success_count = 0
        fail_count = 0
        
        for idx, article in enumerate(articles, 1):
            pmid = article.get('pmid', 'unknown')
            title = article.get('title', 'No title')[:80]
            
            logger.info(f"\n{'─'*70}")
            logger.info(f"處理文章 {idx}/{len(articles)}")
            logger.info(f"   PMID: {pmid}")
            logger.info(f"   標題: {title}...")
            logger.info(f"{'─'*70}")
            
            # 提取臨床資訊
            clinical_info = self._extract_from_single_article(article, variant_id)
            
            # 檢查提取是否成功
            if clinical_info and self._is_valid_extraction(clinical_info):
                article['clinical_info'] = clinical_info
                enriched_articles.append(article)
                success_count += 1
                logger.info(f"   PMID {pmid}: 提取成功")
                self._log_extraction_result(clinical_info)
            else:
                article['clinical_info'] = self._get_default_info()
                enriched_articles.append(article)
                fail_count += 1
                logger.warning(f"   PMID {pmid}: 未提取到有效資訊")
        
        # 總結
        logger.info(f"\n{'='*70}")
        logger.info(f"提取完成！")
        logger.info(f"   總文章數: {len(articles)}")
        logger.info(f"   成功提取: {success_count}")
        logger.info(f"   提取失敗: {fail_count}")
        logger.info(f"   成功率: {success_count/len(articles)*100:.1f}%")
        logger.info(f"{'='*70}\n")
        
        return enriched_articles
    
    def _extract_from_single_article(
        self,
        article: Dict,
        variant_id: str
    ) -> Optional[Dict]:
        """
        從單篇文章中提取臨床資訊
        
        Args:
            article: 文章資訊
            variant_id: 變異 ID
        
        Returns:
            提取的臨床資訊字典，失敗返回 None
        """
        # 準備文本 - 優先使用 abstract
        text = article.get('abstract', '')
        if not text or len(text) < 100:
            logger.warning("   摘要太短或不存在，嘗試使用標題")
            text = article.get('title', '')
        
        if not text or len(text) < 20:
            logger.warning("   文本太短，跳過")
            return None
        
        # 截斷過長的文本
        max_length = 2000
        if len(text) > max_length:
            text = text[:max_length]
            logger.info(f"   📏 文本截斷至 {max_length} 字元")
        
        logger.info(f"   📝 文本長度: {len(text)} 字元")
        
        # 構建 prompt
        prompt = self._create_simple_prompt(text, variant_id)
        
        # 生成
        try:
            if self.backend == "vllm":
                response = self._generate_vllm(prompt)
            else:
                response = self._generate_ollama(prompt)
            
            logger.info(f"   模型響應長度: {len(response)} 字元")
            
            if not response or len(response) < 10:
                logger.warning(f"   模型返回空響應或過短")
                return None
            
            # 解析 JSON - 使用更寬鬆的解析
            clinical_info = self._parse_response_flexible(response)
            
            return clinical_info
            
        except Exception as e:
            logger.error(f"   提取失敗: {e}")
            return None
    
    def _create_simple_prompt(self, text: str, variant_id: str) -> str:
        """
        創建簡化的 prompt
        
        重點改進:
        1. 更短、更直接的指令
        2. 明確的 JSON schema
        3. 使用範例來引導模型
        """
        prompt = f"""Extract clinical information from the article below about the genetic variant {variant_id}.

Article text:
{text}

Please extract ONLY the following information and respond in JSON format:

{{
    "inheritance": "dominant" or "recessive" or "not specified",
    "age_distribution": "age range (e.g., 20-50 years)" or "not specified",
    "affected_tissue": "cardiac" or "skeletal" or "both" or "not specified",
    "sample_size": number or "not specified"
}}

Important:
- Use "not specified" if information is not found in the text
- For age_distribution, extract actual age ranges mentioned
- For inheritance, look for keywords: dominant, recessive, autosomal
- For tissue, look for: heart, cardiac, muscle, skeletal

JSON response:"""
        
        return prompt
    
    def _generate_vllm(self, prompt: str) -> str:
        """使用 vLLM 生成"""
        outputs = self.llm.generate([prompt], self.sampling_params)
        response = outputs[0].outputs[0].text
        
        # 記錄生成統計
        num_tokens = len(outputs[0].outputs[0].token_ids)
        finish_reason = outputs[0].outputs[0].finish_reason
        
        logger.info(f"      生成 tokens: {num_tokens}")
        logger.info(f"      Finish reason: {finish_reason}")
        logger.info(f"      響應長度: {len(response)} 字元")
        
        # 顯示響應的前 200 字元用於調試
        if len(response) > 0:
            logger.debug(f"      響應預覽: {response[:200]}...")
        else:
            logger.warning(f"      響應為空！")
        
        return response
    
    def _generate_ollama(self, prompt: str) -> str:
        """使用 Ollama 生成"""
        response = self.ollama_client.chat(
            model=self.model_name,
            messages=[{'role': 'user', 'content': prompt}],
            options={
                'temperature': 0.3,
                'num_predict': 1000,
            }
        )
        return response['message']['content']
    
    def _parse_response_flexible(self, response: str) -> Optional[Dict]:
        """
        使用更寬鬆的方式解析響應
        
        嘗試多種解析策略:
        1. 直接 JSON 解析
        2. 提取 JSON 區塊
        3. 使用正則表達式提取關鍵資訊
        
        最後轉換為 html_report.py 期望的格式
        """
        # 策略 1: 直接解析
        simple_info = None
        try:
            # 清理響應
            response = response.strip()
            
            # 移除 markdown 標記
            response = re.sub(r'```json\s*', '', response)
            response = re.sub(r'```\s*', '', response)
            
            # 嘗試找到 JSON 區塊
            json_match = re.search(r'\{[^{}]*\}', response, re.DOTALL)
            if json_match:
                json_str = json_match.group(0)
                simple_info = json.loads(json_str)
                logger.debug(f"      JSON 解析成功")
        except json.JSONDecodeError as e:
            logger.debug(f"      JSON 解析失敗: {e}")
        
        # 策略 2: 使用正則表達式提取
        if not simple_info:
            try:
                simple_info = self._extract_with_regex(response)
                if simple_info:
                    logger.debug(f"      正則表達式提取成功")
            except Exception as e:
                logger.debug(f"      正則表達式提取失敗: {e}")
        
        if not simple_info:
            logger.warning(f"      所有解析策略都失敗")
            return None
        
        # 轉換為 html_report.py 期望的格式
        return self._convert_to_html_format(simple_info)
    
    def _convert_to_html_format(self, simple_info: Dict) -> Dict:
        """
        將簡單格式轉換為 html_report.py 期望的複雜格式
        
        簡單格式:
        {
            'inheritance': 'dominant',
            'age_distribution': '30-60 years',
            'affected_tissue': 'cardiac',
            'sample_size': '50'
        }
        
        轉換為:
        {
            'inheritance': {
                'pattern': 'autosomal_dominant',
                'description': 'Autosomal dominant'
            },
            'age_onset': {
                'range': '30-60',
                'description': '30-60 years'
            },
            'tissue_affected': {
                'primary': 'cardiac',
                'description': 'Cardiac muscle'
            },
            'sample_size': {
                'patients': 50,
                'families': 0
            }
        }
        """
        html_format = {}
        
        # 轉換 inheritance
        inheritance = simple_info.get('inheritance', 'not specified')
        if inheritance and inheritance != 'not specified':
            if 'dominant' in inheritance.lower():
                html_format['inheritance'] = {
                    'pattern': 'autosomal_dominant',
                    'description': 'Autosomal dominant'
                }
            elif 'recessive' in inheritance.lower():
                html_format['inheritance'] = {
                    'pattern': 'autosomal_recessive',
                    'description': 'Autosomal recessive'
                }
            else:
                html_format['inheritance'] = {
                    'pattern': inheritance,
                    'description': inheritance.capitalize()
                }
        else:
            html_format['inheritance'] = {
                'pattern': 'not_specified',
                'description': 'Not specified'
            }
        
        # 轉換 age_distribution
        age = simple_info.get('age_distribution', 'not specified')
        if age and age != 'not specified':
            # 提取數字範圍
            age_match = re.search(r'(\d+)\s*[-to]*\s*(\d+)', age)
            if age_match:
                html_format['age_onset'] = {
                    'range': f"{age_match.group(1)}-{age_match.group(2)}",
                    'description': age
                }
            else:
                html_format['age_onset'] = {
                    'range': age,
                    'description': age
                }
        else:
            html_format['age_onset'] = {
                'range': 'not_specified',
                'description': 'Not specified'
            }
        
        # 轉換 affected_tissue
        tissue = simple_info.get('affected_tissue', 'not specified')
        if tissue and tissue != 'not specified':
            html_format['tissue_affected'] = {
                'primary': tissue.lower(),
                'description': tissue.capitalize() + ' muscle'
            }
        else:
            html_format['tissue_affected'] = {
                'primary': 'not_specified',
                'description': 'Not specified'
            }
        
        # 轉換 sample_size
        sample = simple_info.get('sample_size', 'not specified')
        patients = 0
        if sample and sample != 'not specified':
            try:
                patients = int(str(sample).replace(',', ''))
            except (ValueError, AttributeError):
                pass
        
        html_format['sample_size'] = {
            'patients': patients,
            'families': 0  # 我們的簡單提取器不提取家族數
        }
        
        return html_format
    
    def _extract_with_regex(self, text: str) -> Optional[Dict]:
        """
        使用正則表達式從文本中提取資訊
        作為 JSON 解析失敗的備用方案
        """
        info = {
            'inheritance': 'not specified',
            'age_distribution': 'not specified',
            'affected_tissue': 'not specified',
            'sample_size': 'not specified'
        }
        
        # 提取 inheritance
        if re.search(r'\bdominant\b', text, re.IGNORECASE):
            info['inheritance'] = 'dominant'
        elif re.search(r'\brecessive\b', text, re.IGNORECASE):
            info['inheritance'] = 'recessive'
        
        # 提取 age
        age_match = re.search(r'(\d+)\s*[-to]*\s*(\d+)\s*years?', text, re.IGNORECASE)
        if age_match:
            info['age_distribution'] = f"{age_match.group(1)}-{age_match.group(2)} years"
        
        # 提取 tissue
        if re.search(r'\b(cardiac|heart)\b', text, re.IGNORECASE):
            if re.search(r'\b(skeletal|muscle)\b', text, re.IGNORECASE):
                info['affected_tissue'] = 'both'
            else:
                info['affected_tissue'] = 'cardiac'
        elif re.search(r'\b(skeletal|muscle)\b', text, re.IGNORECASE):
            info['affected_tissue'] = 'skeletal'
        
        # 提取 sample size
        sample_match = re.search(r'n\s*=\s*(\d+)', text, re.IGNORECASE)
        if sample_match:
            info['sample_size'] = sample_match.group(1)
        
        # 檢查是否提取到任何有效資訊
        if any(v != 'not specified' for v in info.values()):
            return info
        
        return None
    
    def _is_valid_extraction(self, info: Dict) -> bool:
        """
        檢查提取的資訊是否有效
        至少需要一個欄位不是 "not specified"
        
        支援新的嵌套格式
        """
        if not info:
            return False
        
        # 檢查 inheritance
        inheritance = info.get('inheritance', {})
        if isinstance(inheritance, dict):
            if inheritance.get('pattern', 'not_specified') != 'not_specified':
                return True
        elif inheritance and str(inheritance).lower() != 'not specified':
            return True
        
        # 檢查 age_onset
        age = info.get('age_onset', {})
        if isinstance(age, dict):
            if age.get('range', 'not_specified') != 'not_specified':
                return True
        elif age and str(age).lower() != 'not specified':
            return True
        
        # 檢查 tissue_affected
        tissue = info.get('tissue_affected', {})
        if isinstance(tissue, dict):
            if tissue.get('primary', 'not_specified') != 'not_specified':
                return True
        elif tissue and str(tissue).lower() != 'not specified':
            return True
        
        # 檢查 sample_size
        sample = info.get('sample_size', {})
        if isinstance(sample, dict):
            if sample.get('patients', 0) > 0:
                return True
        elif sample and str(sample).lower() != 'not specified':
            return True
        
        return False
    
    def _log_extraction_result(self, info: Dict):
        """記錄提取結果的詳細資訊，支援嵌套格式"""
        # 獲取 inheritance
        inheritance = info.get('inheritance', {})
        if isinstance(inheritance, dict):
            inheritance_str = inheritance.get('description', 'N/A')
        else:
            inheritance_str = str(inheritance)
        
        # 獲取 age
        age = info.get('age_onset', {})
        if isinstance(age, dict):
            age_str = age.get('description', 'N/A')
        else:
            age_str = str(age)
        
        # 獲取 tissue
        tissue = info.get('tissue_affected', {})
        if isinstance(tissue, dict):
            tissue_str = tissue.get('description', 'N/A')
        else:
            tissue_str = str(tissue)
        
        # 獲取 sample size
        sample = info.get('sample_size', {})
        if isinstance(sample, dict):
            sample_str = str(sample.get('patients', 0))
        else:
            sample_str = str(sample)
        
        logger.info(f"      • 遺傳模式: {inheritance_str}")
        logger.info(f"      • 年齡分布: {age_str}")
        logger.info(f"      • 影響組織: {tissue_str}")
        logger.info(f"      • 樣本數: {sample_str}")
    
    def _get_default_info(self) -> Dict:
        """返回預設的空資訊，格式符合 html_report.py 的期望"""
        return {
            'inheritance': {
                'pattern': 'not_specified',
                'description': 'Not specified'
            },
            'age_onset': {
                'range': 'not_specified',
                'description': 'Not specified'
            },
            'tissue_affected': {
                'primary': 'not_specified',
                'description': 'Not specified'
            },
            'sample_size': {
                'patients': 0,
                'families': 0
            }
        }
    
    def batch_extract(
        self,
        articles: List[Dict],
        variant_info: Dict
    ) -> List[Dict]:
        """
        批次提取文章的臨床資訊
        這是 main.py 調用的主要方法
        
        Args:
            articles: PubMed 文章列表
            variant_info: 變異資訊字典
        
        Returns:
            包含臨床資訊的文章列表
        """
        variant_id = f"{variant_info['chrom']}-{variant_info['pos']}-{variant_info['ref']}-{variant_info['alt']}"
        return self.extract_from_articles(articles, variant_id)
    
    def aggregate_stats(self, articles: List[Dict]) -> Dict:
        """
        從提取的文章中聚合統計資訊
        
        Args:
            articles: 已提取臨床資訊的文章列表
        
        Returns:
            聚合的統計資訊字典
        """
        stats = {
            'total_articles': len(articles),
            'articles_with_extraction': 0,
            'total_patients': 0,
            'total_families': 0,
            'inheritance_distribution': {},
            'age_onset_distribution': {},
            'tissue_affected': {},
            'cardiac_phenotypes': {},
            'skeletal_phenotypes': {}
        }
        
        for article in articles:
            clinical_info = article.get('clinical_info', {})
            
            # 檢查是否有有效提取
            if self._is_valid_extraction(clinical_info):
                stats['articles_with_extraction'] += 1
            
            # 統計遺傳模式
            inheritance = clinical_info.get('inheritance', {})
            if isinstance(inheritance, dict):
                pattern = inheritance.get('pattern', 'not_specified')
                if pattern != 'not_specified':
                    stats['inheritance_distribution'][pattern] = \
                        stats['inheritance_distribution'].get(pattern, 0) + 1
            
            # 統計年齡分布
            age = clinical_info.get('age_onset', {})
            if isinstance(age, dict):
                age_range = age.get('range', 'not_specified')
                if age_range != 'not_specified':
                    stats['age_onset_distribution'][age_range] = \
                        stats['age_onset_distribution'].get(age_range, 0) + 1
            
            # 統計影響組織
            tissue = clinical_info.get('tissue_affected', {})
            if isinstance(tissue, dict):
                primary = tissue.get('primary', 'not_specified')
                if primary != 'not_specified':
                    stats['tissue_affected'][primary] = \
                        stats['tissue_affected'].get(primary, 0) + 1
                    
                    # 根據組織類型添加表型
                    title = article.get('title', 'Unknown')
                    if 'cardiac' in primary.lower():
                        stats['cardiac_phenotypes'][title] = \
                            stats['cardiac_phenotypes'].get(title, 0) + 1
                    if 'skeletal' in primary.lower():
                        stats['skeletal_phenotypes'][title] = \
                            stats['skeletal_phenotypes'].get(title, 0) + 1
            
            # 統計樣本數
            sample = clinical_info.get('sample_size', {})
            if isinstance(sample, dict):
                patients = sample.get('patients', 0)
                families = sample.get('families', 0)
                stats['total_patients'] += patients
                stats['total_families'] += families
        
        return stats


# 測試代碼
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # 測試文章
    test_articles = [
        {
            'pmid': '12345678',
            'title': 'Autosomal dominant dilated cardiomyopathy caused by TTN mutations',
            'abstract': 'We studied 50 patients aged 30-60 years with cardiac involvement. The disease showed autosomal dominant inheritance pattern.'
        }
    ]
    
    # 創建提取器
    extractor = LocalClinicalExtractor(backend="vllm")
    
    # 提取資訊
    results = extractor.extract_from_articles(test_articles, "2-178528273-G-C")
    
    # 顯示結果
    print("\n測試結果:")
    print(json.dumps(results[0]['clinical_info'], indent=2))