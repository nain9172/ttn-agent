#!/usr/bin/env python3
"""
检查 PMC 页面的实际下载链接
"""

import requests
from bs4 import BeautifulSoup

pmc_url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6315030/"

headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
}

print(f"访问: {pmc_url}\n")
response = requests.get(pmc_url, headers=headers)

soup = BeautifulSoup(response.content, 'html.parser')

print("查找所有包含 'bin/' 的链接：\n")
found_count = 0

for link in soup.find_all('a', href=True):
    href = link['href']
    if '/bin/' in href and any(ext in href for ext in ['.pdf', '.xlsx', '.xls']):
        text = link.get_text(strip=True)
        print(f"文本: {text}")
        print(f"链接: {href}")
        print()
        found_count += 1
        
        if found_count >= 5:  # 只显示前5个
            break

print(f"总共找到 {found_count} 个下载链接")
