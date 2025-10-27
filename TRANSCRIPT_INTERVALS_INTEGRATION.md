# Transcript Intervals Integration 整合說明

## 概要 Overview

成功將 transcript intervals 繪製功能整合到 TTN Variant AI Agent 系統中。

## 新增功能 New Features

### 1. ImageGenerator 類新增方法

在 `utils/image_generator.py` 中新增以下方法：

- `generate_transcript_intervals()`: 從 Excel 檔案讀取並繪製 transcript intervals
- `_parse_interval()`: 解析區間字串
- `_normalize_genomic_position()`: 標準化基因組位置
- `_plot_transcript_with_intervals()`: 繪製單個 transcript 的區間

### 2. Main Pipeline 更新

在 `main.py` 中新增 **Step 4.5**：
- 自動檢測 `transcript_interval.xlsx` 檔案
- 如果檔案存在，生成 transcript intervals 圖
- 將變異位置標記在圖上

### 3. HTML Report 增強

在 `utils/html_report.py` 中：
- 新增 `_generate_transcript_intervals_section()` 方法
- 在報告中顯示 transcript intervals 圖（位於 Protein Domain Localization 之後）
- 提供詳細的說明文字

## 使用方式 Usage

### 獨立使用 transcript intervals 繪圖

```python
from utils.image_generator import ImageGenerator

# 不標記變異位置
ig = ImageGenerator()
output_path = ig.generate_transcript_intervals('transcript_interval.xlsx')

# 標記變異位置
variant_info = {
    'variant_id': '2-178612477-T-A',
    'pos': 178612477,
    # ... 其他欄位
}
output_path = ig.generate_transcript_intervals('transcript_interval.xlsx', variant_info)
```

### 在完整 pipeline 中使用

```bash
# 一般使用（會自動包含 transcript intervals）
python main.py 2-178612477-T-A

# 產生的報告會包含：
# 1. Protein Domain Localization 圖
# 2. Transcript Intervals 圖（新增）
# 3. 其他分析結果
```

## 檔案需求 File Requirements

需要 `transcript_interval.xlsx` 檔案，格式如下：

| transcript | interval1 | interval2 | ... | interval9 |
|------------|-----------|-----------|-----|-----------|
| Meta | 178807423-178753124 | 178741921-178525989 | ... | NaN |
| N2BA | 178807423-178758984 | 178753180-178753124 | ... | 178646030-178525989 |
| N2B | 178807423-178792072 | 178790115-178758984 | ... | NaN |
| ... | ... | ... | ... | ... |

## 輸出檔案 Output Files

1. **圖片檔案**：`outputs/images/transcript_intervals.png` 或 `transcript_intervals_{variant_id}.png`
2. **HTML 報告**：包含 transcript intervals 圖的完整報告

## 技術細節 Technical Details

### 座標系統
- 使用基因組座標（chr2）
- TTN 基因位於負鏈（negative strand）
- 座標範圍：178,807,423 - 178,525,989

### 顏色配置
使用 9 種不同顏色標記不同的 intervals：
- #FF6B6B (紅)
- #4ECDC4 (青)
- #45B7D1 (藍)
- #96CEB4 (綠)
- #FECA57 (黃)
- #FF9FF3 (粉)
- #54A0FF (天藍)
- #48DBFB (水藍)
- #00D2D3 (青綠)

### 變異標記
- 紅色箭頭標記變異位置
- 自動計算在基因組座標上的位置
- 在所有 transcript 圖上顯示相同位置

## 測試 Testing

```bash
# 測試 transcript intervals 繪圖功能
python3 -c "from utils.image_generator import ImageGenerator; ig = ImageGenerator(); print(ig.generate_transcript_intervals('transcript_interval.xlsx'))"

# 測試完整 pipeline
python main.py 2-178612477-T-A --skip-evo2
```

## 相關檔案 Related Files

- `utils/image_generator.py` - 圖片生成器（已更新）
- `utils/html_report.py` - HTML 報告生成器（已更新）
- `main.py` - 主程式（已更新）
- `plot_transcript_intervals.py` - 獨立繪圖腳本（參考用）
- `transcript_interval.xlsx` - 資料檔案

## 更新日期 Update Date

2025-10-27

