# TTN Agent 改進總結

## 完成的改進項目

### 1. 全文獲取功能
**問題：** LLM 處理速度太快，可能沒有讀取完整文章
**解決方案：**
- 修改 `batch_extract` 方法，優先使用完整文本而不是只用摘要
- 增加 vLLM 的 `max_model_len` 從 4096 到 8192
- 調整文本截斷限制從 3000 到 5000 字元
- 全文獲取功能已啟用（從 PMC, Europe PMC, Unpaywall）

**修改文件：**
- `utils/local_clinical_extractor.py`
- `utils/enhanced_pubmed_search.py`

### 2. 報告顯示優化
**問題：** 需要顯示病患年齡分布和影響組織，不需要顯示文章文字
**解決方案：**
- 在每篇文章顯示：
  - Age Distribution (年齡分布)：包含平均年齡、中位數年齡、年齡範圍
  - Affected Tissue (影響組織)：顯示心肌/骨骼肌及具體表型
  - Inheritance (遺傳模式)：顯性/隱性等
- 移除文章摘要文字的顯示

**修改文件：**
- `utils/html_report.py`

### 3. LLM Prompt 優化
**問題：** Llama-3.2-3B 模型只複製範例數據，所有文章提取結果相同
**解決方案：**
- 重新設計 prompt，移除具體數字範例
- 強調從文章提取實際數據
- 使用佔位符而不是具體範例值
- 建議使用 Llama-3.1-8B-Instruct 以獲得更好的提取效果

**修改文件：**
- `utils/local_clinical_extractor.py`
- `config.py` (建議模型：Llama-3.1-8B-Instruct)

### 4. 移除所有表情符號
**問題：** 用戶要求移除報告和日誌中的表情符號
**解決方案：**
- 移除 HTML 報告中的所有表情符號
- 移除日誌輸出中的所有表情符號

**修改文件：**
- `utils/html_report.py`
- `main.py`
- `utils/local_clinical_extractor.py`
- `utils/enhanced_pubmed_search.py`
- `utils/clinvar_parser.py`

## 配置建議

### 推薦配置
```python
# config.py
ENABLE_FULL_TEXT_FETCH = True
MAX_TEXT_LENGTH = 8000
LOCAL_LLM_MODEL = "meta-llama/Llama-3.1-8B-Instruct"
LOCAL_LLM_TENSOR_PARALLEL = 2
```

### 模型選擇
- **Llama-3.2-3B**: 快速但提取準確度較低（不推薦用於臨床資訊提取）
- **Llama-3.1-8B-Instruct**: 推薦！較好的指令理解和資訊提取能力
- **Llama-3.1-70B-Instruct**: 最佳效果但需要更多記憶體

## 使用說明

### 測試全文獲取
```bash
python3 test_full_text.py
```

### 運行完整分析
```bash
python3 main.py 2-178612477-T-A
```

## 報告顯示範例

每篇文章將顯示：

```
1. [文章標題]
   作者 | 期刊 (年份) | PMID: xxxxx

   Age Distribution (年齡分布): Adult (成人) (平均 45 歲, 範圍 20-65)
   Affected Tissue (影響組織): Cardiac (心肌): DCM, HCM
   Inheritance (遺傳模式): Autosomal Dominant (顯性)
   
   Sample Size: 120 patients, 35 families
   Severity: Moderate
```

## 注意事項

1. **全文獲取**：並非所有文章都能獲取全文，系統會自動回退到摘要
2. **模型性能**：Llama-3.1-8B-Instruct 需要約 16GB GPU 記憶體（雙 4090 可順利運行）
3. **提取準確度**：AI 提取的資訊僅供參考，建議閱讀原文確認
4. **處理時間**：完整文章分析需要更多時間，但提取質量更好

## 問題排查

### 如果提取結果仍然相同
1. 確認使用 Llama-3.1-8B-Instruct 而不是 3.2-3B
2. 檢查文章是否成功獲取全文（日誌中會顯示）
3. 查看日誌中的 "模型響應長度" 是否合理

### 如果記憶體不足
1. 降低 `gpu_memory_utilization` (目前 0.85)
2. 減少 `max_model_len` (目前 8192)
3. 使用單 GPU (`LOCAL_LLM_TENSOR_PARALLEL = 1`)

## 後續建議

1. 定期驗證 AI 提取的準確度
2. 考慮為不同類型文章（case report, cohort study 等）使用不同 prompt
3. 收集用戶反饋持續優化 prompt

