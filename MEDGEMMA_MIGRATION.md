# MedGemma 遷移指南

## 概述
已成功將系統從 Llama 3.2 遷移到 Google MedGemma-27B-IT，這是一個專為醫療領域優化的 27B 參數模型。

## 主要變更

### 1. 模型配置 (config.py)

```python
# 模型選擇
LOCAL_LLM_MODEL = "google/medgemma-27b-it"  # 從 Llama 3.2-3B 升級到 MedGemma-27B

# 上下文長度大幅提升
LOCAL_LLM_MAX_MODEL_LEN = 65536  # 從 32K 提升到 65K tokens
LOCAL_LLM_MAX_CONTEXT_LENGTH = 50000  # 從 24K 提升到 50K 字元
MAX_TEXT_LENGTH = 80000  # 從 50K 提升到 80K 字元
DOCLING_MAX_PRIORITY_LENGTH = 120000  # 從 80K 提升到 120K 字元
```

### 2. Prompt 格式更新 (local_clinical_extractor.py)

**Llama 3.2 格式（舊）：**
```
<|begin_of_text|><|start_header_id|>system<|end_header_id|>
...
<|eot_id|><|start_header_id|>user<|end_header_id|>
...
<|eot_id|><|start_header_id|>assistant<|end_header_id|>
```

**MedGemma 格式（新）：**
```
<start_of_turn>user
{system_message}

{user_message}<end_of_turn>
<start_of_turn>model
```

系統會自動檢測模型類型並使用正確的格式。

### 3. vLLM 參數優化

**MedGemma 特定優化：**
- GPU 記憶體利用率：90%（Llama 為 85%）
- Temperature：0.3（更保守，適合醫療任務，Llama 為 0.7）
- Stop tokens：`["<end_of_turn>", "\n\n\n"]`（Gemma 專用）

### 4. 自動模型檢測

系統現在會自動檢測模型類型：
- 檢測 "medgemma" 或 "gemma" → 使用 Gemma chat template
- 檢測 "llama" → 使用 Llama chat template
- 日誌會顯示檢測結果

## MedGemma 優勢

### 1. 醫療領域專業化
- 在 22+ 個醫療數據集上訓練
- 涵蓋放射學、病理學、皮膚科、眼科等
- 在 MedQA 上達到 89.8% 準確率（best-of-5）

### 2. 更長的上下文窗口
- 支持最高 128K tokens（Llama 3.2 為 128K，但實際使用上 MedGemma 更穩定）
- 可處理更長的醫學文獻和病歷

### 3. 更好的醫學推理
- 使用 test-time scaling 優化
- 在醫學推理任務上表現優異

### 4. 多模態能力
- 支持醫學影像（X-ray、病理切片等）
- 雖然目前只使用文本功能，但未來可擴展

## 性能對比

| 指標 | Llama 3.2-3B | MedGemma-27B |
|------|--------------|--------------|
| 參數量 | 3B | 27B |
| Context 長度 | 32K tokens | 65K+ tokens |
| MedQA 準確率 | ~50% | 89.8% |
| 醫療專業化 | 通用模型 | 醫療專用 |
| GPU 記憶體 | ~12GB | ~54GB (FP16) |

## 系統需求

### 硬體需求
- **GPU 記憶體**：建議 80GB+（A100 或 H100）
- **系統記憶體**：128GB（已滿足）
- **儲存空間**：~54GB（模型大小）

### 軟體需求
- vLLM 0.15.0+
- transformers 4.50.0+
- PyTorch 2.0+
- CUDA 12.0+

## 使用方式

### 運行模型
```bash
python main.py 2-178527121-T-C
```

系統會自動：
1. 檢測到 MedGemma 模型
2. 使用正確的 prompt 格式
3. 應用 MedGemma 優化參數
4. 記錄檢測結果到日誌

### 查看日誌
```
INFO - Initializing vllm backend with model: google/medgemma-27b-it
INFO - Max model length: 65536, Max context length: 50000
INFO - Detected MedGemma/Gemma model - using Gemma chat template
INFO - Using MedGemma-optimized vLLM settings
```

## 切換回 Llama（如需）

只需修改 `config.py`：

```python
LOCAL_LLM_MODEL = "meta-llama/Llama-3.2-3B-Instruct"
LOCAL_LLM_MAX_MODEL_LEN = 32768
LOCAL_LLM_MAX_CONTEXT_LENGTH = 24000
```

系統會自動檢測並切換回 Llama 格式。

## 注意事項

1. **首次運行**：第一次會下載 ~54GB 的模型文件
2. **記憶體使用**：確保 GPU 有足夠空閒記憶體（建議至少 60GB 可用）
3. **推理速度**：27B 模型比 3B 模型慢，但準確率更高
4. **Temperature**：MedGemma 使用較低的 temperature (0.3) 以提高醫療任務的準確性

## 效能優化建議

### 如果 GPU 記憶體不足
```python
# 減少 max_model_len
LOCAL_LLM_MAX_MODEL_LEN = 32768  # 從 65536 降低

# 或使用 4B 版本
LOCAL_LLM_MODEL = "google/medgemma-4b-it"
```

### 如果需要更快的推理
```python
# 使用多 GPU
LOCAL_LLM_TENSOR_PARALLEL = 2  # 或 4

# 或使用量化版本（需要額外設定）
```

## 參考資源

- [MedGemma 官方文檔](https://developers.google.com/health-ai-developer-foundations/medgemma)
- [Hugging Face 模型頁面](https://huggingface.co/google/medgemma-27b-it)
- [MedGemma Technical Report](https://arxiv.org/abs/2507.05201)
- [GitHub Repository](https://github.com/google-health/medgemma)

---
*更新日期：2026-02-04*
