# Integrated Visualization Update 整合視覺化更新

## 更新日期
2025-10-27

## 更新概要

成功將 Protein Domain Localization 與 Transcript Intervals 整合為單一視覺化圖片，所有元素完全對齊。

## 主要變更

### 1. **整合視覺化設計**

#### 之前：
- Protein Domain Localization（獨立圖片）
  - Domain Overview (Z-disk, I-band, A-band, M-band)
  - 各個 transcript 的蛋白質結構視圖
- Transcript Intervals（獨立圖片）
  - 各個 transcript 的基因組區間

#### 現在：
- **單一整合圖片**
  - **Top Panel**: Domain Overview (Z-disk, I-band, A-band, M-band) - 保留
  - **Subsequent Panels**: Transcript Intervals（取代原本的 transcript 視圖）
    - 每個 transcript 顯示其真實的基因組區間（來自 xlsx 檔案）
    - 彩色區塊代表不同的 intervals
    - 所有視圖使用相同的座標系統對齊

### 2. **座標對齊系統**

所有視圖現在使用統一的標準化座標系統：
- X 軸範圍：0 到 1（標準化基因組位置）
- 對應基因組座標：178,807,423 到 178,525,989 (chr2, 負鏈)
- 變異位置在所有視圖中垂直對齊

### 3. **視覺元素**

#### Domain Overview（頂部）
- 顯示 4 個主要 protein domains
- 顏色編碼：
  - Z-disk: 紅色 (#FF6B6B)
  - I-band: 青色 (#4ECDC4)
  - A-band: 藍色 (#45B7D1)
  - M-band: 綠色 (#96CEB4)

#### Transcript Intervals（各行）
- 左側：Transcript 名稱標籤
- 中間：彩色區塊代表基因組 intervals
- 右側：Interval 數量統計
- 紅色箭頭：變異位置標記（所有視圖對齊）

## 技術實作

### 修改的檔案

#### 1. `utils/image_generator.py`
新增/修改的方法：
- `generate_titin_schematic()`: 新增 `xlsx_file` 參數
- `_load_transcript_intervals()`: 從 Excel 讀取 transcript intervals
- `_plot_transcript_intervals_aligned()`: 繪製對齊的 transcript intervals

關鍵技術：
```python
# 統一的座標標準化方法
def _normalize_genomic_position(self, pos: int) -> float:
    gene_length = TTN_GENE_INFO['start'] - TTN_GENE_INFO['end']
    return (TTN_GENE_INFO['start'] - pos) / gene_length

# 所有視圖使用相同的 x 軸範圍
ax.set_xlim(-0.05, 1.05)
```

#### 2. `main.py`
- 自動檢測 `transcript_interval.xlsx` 檔案
- 將檔案路徑傳遞給 `generate_titin_schematic()`
- 移除獨立的 transcript intervals 生成步驟

#### 3. `utils/html_report.py`
- 更新圖片說明文字
- 移除獨立的 transcript intervals section
- 簡化報告結構

## 使用方式

### 標準使用（推薦）
```bash
# 確保 transcript_interval.xlsx 存在於專案根目錄
python main.py 2-178612477-T-A
```

### 不使用 transcript intervals（降級模式）
```bash
# 移除或重命名 transcript_interval.xlsx
# 系統會自動使用原始的 transcript 視圖
python main.py 2-178612477-T-A
```

### 程式化使用
```python
from utils.image_generator import ImageGenerator

variant_info = {
    'variant_id': '2-178612477-T-A',
    'chrom': '2',
    'pos': 178612477,
    'ref': 'T',
    'alt': 'A'
}

ig = ImageGenerator()

# 使用 transcript intervals
image_path = ig.generate_titin_schematic(
    variant_info, 
    'transcript_interval.xlsx'
)

# 不使用 transcript intervals
image_path = ig.generate_titin_schematic(variant_info)
```

## 輸出範例

### 生成的圖片包含：
1. **標題**：變異 ID 和基因組位置
2. **Domain Overview**：顯示 Z-disk, I-band, A-band, M-band
3. **Meta transcript**：2 個 intervals
4. **N2BA transcript**：9 個 intervals
5. **N2B transcript**：5 個 intervals
6. **N2A transcript**：9 個 intervals
7. **Nvx1 transcript**：5 個 intervals
8. **Nvx2 transcript**：5 個 intervals
9. **Nvx3 transcript**：2 個 intervals

每個 transcript 都顯示：
- 彩色編碼的 intervals
- Interval 編號（如果空間足夠）
- 變異位置的紅色箭頭標記
- Interval 總數

## 優點

### 1. **空間效率**
- 單一圖片取代兩個獨立圖片
- HTML 報告更簡潔

### 2. **視覺一致性**
- 所有元素完美對齊
- 變異位置在所有視圖中垂直對齊
- 統一的座標系統

### 3. **資訊整合**
- Domain 結構 + Transcript intervals 在同一視圖
- 更容易理解變異位置與結構的關係

### 4. **向下相容**
- 如果 xlsx 檔案不存在，自動降級到原始視圖
- 不會破壞現有功能

## 資料需求

需要 `transcript_interval.xlsx` 檔案，格式：

| transcript | interval1 | interval2 | interval3 | ... |
|------------|-----------|-----------|-----------|-----|
| Meta | 178807423-178753124 | 178741921-178525989 | NaN | ... |
| N2BA | 178807423-178758984 | 178753180-178753124 | 178741921-178677621 | ... |
| ... | ... | ... | ... | ... |

## 測試

```bash
# 測試整合視覺化
python3 -c "
from utils.image_generator import ImageGenerator

variant_info = {
    'variant_id': '2-178612477-T-A',
    'chrom': '2',
    'pos': 178612477,
    'ref': 'T',
    'alt': 'A'
}

ig = ImageGenerator()
result = ig.generate_titin_schematic(variant_info, 'transcript_interval.xlsx')
print(f'Success: {result}')
"

# 測試完整 pipeline
python main.py 2-178612477-T-A --skip-evo2
```

## 相關檔案

- `utils/image_generator.py` - 核心視覺化邏輯
- `utils/html_report.py` - HTML 報告生成
- `main.py` - 主程式 pipeline
- `transcript_interval.xlsx` - 資料來源
- `TRANSCRIPT_INTERVALS_INTEGRATION.md` - 先前的整合文件（已過時）

## 後續可能的增強

1. 在 domain overview 中標記每個 domain 的基因組座標
2. 添加 exon 編號標記
3. 支援更多 transcript isoforms
4. 添加可選的放大視圖（focus on variant region）
5. 支援批次處理多個變異

## 變更歷史

- 2025-10-27: 初始整合 - 將 protein domains 和 transcript intervals 整合為單一對齊視圖

