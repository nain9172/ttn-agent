# ✅ Supplementary Materials 已整合到 main.py

## 🎯 回答您的问题

**是的！** 当您执行 `main.py` 时，supplementary materials（补充材料）**会自动下载并输入到 LLM 的 prompt 中**。

## 📊 完整流程说明

### 1. 配置启用（已设置）

在 `config.py` 中：
```python
DOWNLOAD_SUPPLEMENTARY_FILES = True  # ✅ 已启用
ENABLE_DOCLING_PDF = True            # ✅ 需要启用才能处理 supplementary
```

### 2. 下载和处理流程

当您运行 `main.py` 时，整个流程如下：

```
main.py
  ↓
Step 4: EnhancedPubMedSearcher.search()
  ↓
  调用 docling_processor.process_pdf_for_llm()
    ↓
    1. 下载 PDF 文件
    2. 使用 Docling 提取内容
    3. 📥 下载 Supplementary Files (Excel, PDF 等)
       ↓
       使用 selenium_downloader (您的 sele.py 方法)
       ↓
       下载成功！
    4. 提取 Excel 表格内容转为 Markdown
    5. 创建 priority_content（优先级顺序）
       ↓
       最高优先级：Supplementary Data 📊
       第二优先：Tables from PDF
       第三优先：Results Section
  ↓
返回 article['text_for_llm'] = priority_content
  ↓
Step 4.5: LocalClinicalExtractor.batch_extract()
  ↓
  读取 article['text_for_llm']
  ↓
  构建 prompt（包含 supplementary 内容！）
  ↓
  发送到 LLM 分析
  ↓
  提取临床信息
```

### 3. 优先级排序

`docling_pdf_processor.py` 中的 `create_priority_content` 方法按以下优先级组织内容：

```python
priority_parts = [
    "## Supplementary Data (Extracted Tables)",  # 🥇 最优先！
    "  - Excel 表格内容",
    "  - CSV 数据",
    "",
    "## Tables from Main Article",              # 🥈 第二优先
    "  - PDF 内的表格",
    "",
    "## Results Section",                        # 🥉 第三优先
    "  - 研究结果段落"
]
```

### 4. 实际使用的代码位置

#### A. 下载 Supplementary Files
`utils/docling_pdf_processor.py` 第 428-475 行：

```python
if download_supplementary and pmid and supplementary_links:
    from utils.supplementary_downloader import SupplementaryDownloader
    downloader = SupplementaryDownloader(use_selenium=True)  # ← 使用您的 sele.py 方法
    
    # 下载 Excel 文件
    for link in excel_links:
        local_path = downloader.download_supplementary_file(file_url, pmid=pmid)
        
        # 提取 Excel 表格
        excel_tables = downloader.extract_tables_from_excel(local_path)
        supplementary_content.append(excel_tables)
```

#### B. 添加到 Priority Content
`utils/docling_pdf_processor.py` 第 332-347 行：

```python
# 1. 最优先：Supplementary Data 实际内容
if supplementary_content:
    supp_data_text = "## Supplementary Data (Extracted Tables)\n\n"
    for content in supplementary_content:
        supp_data_text += f"### Supplementary File\n{content}\n\n"
    priority_parts.append(supp_data_text)  # ← 添加到最前面！
```

#### C. LLM 使用此内容
`utils/local_clinical_extractor.py` 第 122 行：

```python
def _extract_single(self, article: Dict, variant_id: str, clinvar_info: Optional[Dict]) -> Dict:
    text = article.get('text_for_llm') or ...  # ← 优先使用 text_for_llm
    prompt = self._build_cot_prompt(text, aliases)  # ← 构建 prompt
    raw = self._generate_vllm(prompt)  # ← 发送到 LLM
```

## 📝 验证方法

### 方法 1：查看日志

运行 `main.py` 时，您会看到类似的日志：

```
INFO - PMID 30602777: 嘗試下載 9 個 supplementary files...
INFO -   → 正在下載: Supplementary Data 1 (excel)
INFO -   ✓ 下載成功: /path/to/file.xlsx
INFO -   → 正在提取 Excel 表格...
INFO -   ✓ 提取了 3 個表格來自 Supplementary Data 1
INFO - PMID 30602777: 成功解析 1 個 supplementary files
INFO - 創建優先內容: 15234 字元 (優先級: Supplementary Data > Tables > Results)
```

### 方法 2：查看 LLM 日志文件

LLM 的 prompt 和 response 会保存在：

```
llm_logs/
└── 2-178802273-G-A/          # 您的变异ID
    └── 30602777.txt          # 每篇文章的 PMID
```

打开任意一个文件，您会看到：

```
================================ PROMPT ================================

<文章内容，包含：>

## Supplementary Data (Extracted Tables)

*(These tables are extracted from supplementary Excel/CSV files)*

### Supplementary File 1
**Supplementary Data 1**

| Patient ID | Variant | Phenotype | Age Onset |
|------------|---------|-----------|-----------|
| P001       | ...     | ...       | ...       |
...

## Tables from Main Article
...

## Results Section
...
```

### 方法 3：运行测试

```bash
conda activate ai
python test_integrated_download.py
```

应该看到：
```
✅ 下载成功！
   文件: 41467_2018_7709_MOESM6_ESM.xlsx
   大小: 197,405 bytes (192.8 KB)
```

## 🔧 如何确保功能正常

### 检查清单

- [x] **config.py** 中 `DOWNLOAD_SUPPLEMENTARY_FILES = True`
- [x] **config.py** 中 `ENABLE_DOCLING_PDF = True`
- [x] **selenium** 和 **webdriver-manager** 已安装
- [x] **Chrome/Chromium** 浏览器已安装
- [x] **utils/selenium_downloader.py** 使用您的 sele.py 方法

### 运行完整测试

```bash
# 激活环境
conda activate ai

# 运行 main.py（使用一个有补充材料的变异）
python main.py 2-178802273-G-A

# 查看 LLM 日志确认 supplementary 内容
cat llm_logs/2-178802273-G-A/*.txt | grep -A 10 "Supplementary Data"
```

## 📊 实际效果

### 下载的文件类型

- ✅ **Excel 文件** (.xlsx, .xls) - 提取表格数据
- ✅ **PDF 文件** - 下载但不自动解析
- ✅ **CSV 文件** - 可以提取
- ⚠️ **其他格式** - 下载但可能不解析

### 性能

- **每篇文章**: 增加 10-30 秒（下载 + 处理时间）
- **最多下载**: 5 个 supplementary files per article
- **优先下载**: Excel 文件（包含数据表）

## 🎉 总结

**是的，supplementary materials 会自动输入到 prompt！**

整个流程已经完全整合：
1. ✅ 自动检测 supplementary files
2. ✅ 使用 Selenium 下载（您的 sele.py 方法）
3. ✅ 提取 Excel 表格内容
4. ✅ 作为**最高优先级**内容添加到 prompt
5. ✅ 发送到 LLM 进行分析

您不需要做任何额外配置，只需正常运行 `main.py` 即可！

---

*最后更新: 2026-01-06*
