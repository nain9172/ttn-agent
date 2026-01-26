# PubMed Supplementary Materials 下载功能

## 问题说明

之前的下载程序使用 `requests` 库直接下载 PMC 的补充材料，但 PMC 服务器会阻止这种请求，返回 HTML 页面而不是实际文件。

## 解决方案

使用 **Selenium WebDriver** 模拟真实浏览器行为来下载文件。

### 关键技术点

1. **Chrome DevTools Protocol**: 在无头模式下启用文件下载
2. **会话建立**: 先访问文章页面，然后再访问下载链接
3. **正确的标头**: 包括 Referer、User-Agent 等

## 测试结果

### 成功案例

```python
# 下载 PDF 文件 - ✓ 成功
URL: https://pmc.ncbi.nlm.nih.gov/articles/instance/5043576/bin/rsob160114supp1.pdf
文件大小: 3,805,330 bytes (3.7 MB)
```

### URL 格式说明

成功的下载 URL 格式：
```
https://pmc.ncbi.nlm.nih.gov/articles/instance/{instance_id}/bin/{filename}
```

失败的 URL 格式（从 Entrez API 获取）：
```
https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/bin/{filename}
```

**注意**: Entrez API 返回的 URL 格式可能不是实际的下载链接，需要从网页抓取获取正确的 `instance` URL。

## 使用方法

### 基本用法

```python
from utils.supplementary_downloader import SupplementaryDownloader
from pathlib import Path

# 创建下载器（启用 Selenium）
with SupplementaryDownloader(
    output_dir=Path("./downloads"),
    use_selenium=True  # 推荐！
) as downloader:
    
    # 查找补充材料
    links = downloader.scrape_pmc_supplementary_links("30602777")
    
    # 下载文件
    for link in links:
        local_path = downloader.download_supplementary_file(
            url=link['url'],
            pmid="30602777"
        )
        
        if local_path:
            print(f"✓ 下载成功: {local_path}")
```

### 下载特定文件

```python
with SupplementaryDownloader(use_selenium=True) as downloader:
    local_path = downloader.download_supplementary_file(
        url='https://pmc.ncbi.nlm.nih.gov/articles/instance/5043576/bin/rsob160114supp1.pdf',
        filename='my_supplement.pdf',
        pmid='27493940'
    )
```

## 依赖要求

```bash
# 安装 Selenium
conda activate ai
pip install selenium

# 确保系统已安装 Chrome 和 ChromeDriver
# ChromeDriver 会自动由 Selenium 管理
```

## 文件结构

```
utils/
├── supplementary_downloader.py  # 主下载器（支持 Selenium）
└── selenium_downloader.py       # Selenium 下载器封装

测试文件:
├── test_selenium_supplement.py  # Selenium 基础测试
├── test_chrome_download.py      # Chrome DevTools 测试
└── test_final_download.py       # 完整功能测试
```

## 已知问题

1. **URL 格式**: 
   - Entrez API 返回的某些 URL 可能无效
   - 建议使用网页抓取获取实际下载链接

2. **下载速度**:
   - Selenium 比 requests 慢（需要启动浏览器）
   - 每个文件需要 2-5 秒

3. **成功率**:
   - PDF 文件: 高成功率 ✓
   - Excel 文件: 取决于 URL 格式
   - 需要测试更多文章

## 建议优化

1. **改进 URL 获取**:
   - 优先使用网页抓取而不是 Entrez API
   - 解析页面中的实际下载链接

2. **批量下载**:
   - 复用同一个 WebDriver 实例
   - 减少启动开销

3. **错误处理**:
   - 自动重试失败的下载
   - 记录失败的 URL 供后续分析

## 测试命令

```bash
# 测试基本功能
python test_final_download.py

# 测试特定 PMID
python -c "
from utils.supplementary_downloader import SupplementaryDownloader
from pathlib import Path

with SupplementaryDownloader(use_selenium=True) as dl:
    links = dl.scrape_pmc_supplementary_links('YOUR_PMID')
    print(f'找到 {len(links)} 个补充材料')
"
```

## 总结

✓ **Selenium 可以成功下载 PMC 补充材料**
✓ **需要正确的 URL 格式（instance/ 而不是 PMC/）**
✓ **推荐在实际使用中启用 Selenium 模式**

---

*最后更新: 2026-01-05*
