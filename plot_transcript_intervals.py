#!/usr/bin/env python3
"""
繪製 Transcript Interval 圖
根據 transcript_interval.xlsx 繪製各個 transcript 的基因組區間
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np
from pathlib import Path

# 設定字型支援（使用英文避免字型問題）
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# TTN 基因資訊
TTN_GENE_START = 178807423
TTN_GENE_END = 178525989
GENE_LENGTH = TTN_GENE_START - TTN_GENE_END  # 因為是負鏈

# 顏色配置
INTERVAL_COLORS = [
    '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', 
    '#FECA57', '#FF9FF3', '#54A0FF', '#48DBFB',
    '#00D2D3'
]

def parse_interval(interval_str):
    """解析區間字串 '178807423-178758984' 為 (start, end)"""
    if pd.isna(interval_str) or interval_str == '':
        return None
    parts = interval_str.strip().split('-')
    if len(parts) != 2:
        return None
    try:
        return (int(parts[0]), int(parts[1]))
    except:
        return None

def normalize_position(pos):
    """將基因組位置標準化到 0-1 範圍"""
    return (TTN_GENE_START - pos) / GENE_LENGTH

def plot_transcript_intervals(xlsx_file, output_file=None):
    """
    繪製 transcript intervals
    
    Args:
        xlsx_file: xlsx 檔案路徑
        output_file: 輸出圖片檔案路徑 (預設: transcript_intervals.png)
    """
    # 讀取資料
    df = pd.read_excel(xlsx_file)
    print(f"讀取到 {len(df)} 個 transcripts")
    
    # 創建圖形
    num_transcripts = len(df)
    fig_height = 3 + (num_transcripts * 1.5)
    fig = plt.figure(figsize=(20, fig_height))
    
    # 設定子圖布局
    gs = fig.add_gridspec(
        num_transcripts, 
        1, 
        hspace=0.6,
        top=0.95,
        bottom=0.05,
        left=0.08,
        right=0.95
    )
    
    # 主標題
    fig.suptitle(
        'TTN Transcript Intervals',
        fontsize=18,
        fontweight='bold',
        y=0.98
    )
    
    # 繪製每個 transcript
    for idx, row in df.iterrows():
        ax = fig.add_subplot(gs[idx])
        transcript_name = row['transcript']
        
        # 解析所有 intervals
        intervals = []
        for col in df.columns:
            if col.startswith('interval'):
                interval = parse_interval(row[col])
                if interval:
                    intervals.append(interval)
        
        print(f"{transcript_name}: {len(intervals)} intervals")
        
        # 繪製這個 transcript
        plot_single_transcript(ax, transcript_name, intervals)
    
    # 儲存圖片
    if output_file is None:
        output_file = 'outputs/transcript_intervals.png'
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        pad_inches=0.3
    )
    plt.close()
    
    print(f"\n✅ 圖片已儲存至: {output_path}")
    return output_path

def plot_single_transcript(ax, transcript_name, intervals):
    """
    繪製單個 transcript 的 intervals
    
    Args:
        ax: matplotlib axis
        transcript_name: transcript 名稱
        intervals: [(start, end), ...] 區間列表
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    # Transcript 名稱
    ax.text(
        0.02,
        0.75,
        f'{transcript_name}',
        fontsize=12,
        fontweight='bold',
        va='center',
        bbox=dict(
            boxstyle='round,pad=0.5', 
            facecolor='lightblue', 
            edgecolor='navy', 
            alpha=0.8,
            linewidth=2
        )
    )
    
    # 繪製基因組座標軸
    y_pos = 0.45
    height = 0.12
    
    # 背景參考線
    ax.plot([0.12, 0.95], [y_pos, y_pos], 'k-', linewidth=1, alpha=0.3, zorder=1)
    
    # 繪製每個 interval
    for i, (start, end) in enumerate(intervals):
        # 標準化位置 (注意 TTN 在負鏈上)
        start_norm = normalize_position(start)
        end_norm = normalize_position(end)
        
        # 映射到圖上的 x 座標
        x_start = 0.12 + start_norm * 0.83
        x_end = 0.12 + end_norm * 0.83
        width = x_end - x_start
        
        if width > 0:
            # 使用不同顏色
            color = INTERVAL_COLORS[i % len(INTERVAL_COLORS)]
            
            # 繪製區間矩形
            rect = FancyBboxPatch(
                (x_start, y_pos - height/2),
                width,
                height,
                boxstyle="round,pad=0.002",
                edgecolor='black',
                facecolor=color,
                alpha=0.8,
                linewidth=1.5,
                zorder=3
            )
            ax.add_patch(rect)
            
            # 添加區間標籤 (如果區間夠大)
            if width > 0.05:
                ax.text(
                    x_start + width/2,
                    y_pos,
                    f'{i+1}',
                    ha='center',
                    va='center',
                    fontsize=8,
                    fontweight='bold',
                    color='white',
                    zorder=4
                )
    
    # 添加座標軸標記
    # 起始位置
    ax.text(
        0.12,
        y_pos - height/2 - 0.12,
        f'{TTN_GENE_START:,}',
        ha='center',
        va='top',
        fontsize=8,
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.9)
    )
    
    # 結束位置
    ax.text(
        0.95,
        y_pos - height/2 - 0.12,
        f'{TTN_GENE_END:,}',
        ha='center',
        va='top',
        fontsize=8,
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.9)
    )
    
    # 添加中間參考點
    mid_pos = (TTN_GENE_START + TTN_GENE_END) // 2
    ax.text(
        0.535,
        y_pos - height/2 - 0.12,
        f'{mid_pos:,}',
        ha='center',
        va='top',
        fontsize=7,
        color='gray',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='gray', alpha=0.7)
    )
    
    # 添加 interval 數量資訊
    ax.text(
        0.98,
        0.75,
        f'{len(intervals)} intervals',
        ha='right',
        va='center',
        fontsize=9,
        color='darkgreen',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='lightgreen', edgecolor='green', alpha=0.7)
    )

def main():
    """主程式"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='繪製 Transcript Intervals 圖'
    )
    parser.add_argument(
        '--input',
        type=str,
        default='transcript_interval.xlsx',
        help='輸入 xlsx 檔案路徑 (預設: transcript_interval.xlsx)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='輸出圖片檔案路徑 (預設: outputs/transcript_intervals.png)'
    )
    
    args = parser.parse_args()
    
    # 檢查檔案是否存在
    input_file = Path(args.input)
    if not input_file.exists():
        print(f"❌ 錯誤: 找不到檔案 {input_file}")
        return 1
    
    # 繪製圖形
    try:
        plot_transcript_intervals(args.input, args.output)
        return 0
    except Exception as e:
        print(f"❌ 錯誤: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())

