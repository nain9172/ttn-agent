"""
Image Generator Module - FIXED VERSION
Generates protein schematic images showing variant location
Fixed: overlap/clipping issues and ensures variant is only marked on transcripts that actually contain the variant position.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import to_rgba
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

from config import (
    OUTPUT_DIR,
    TTN_DOMAINS,
    TTN_TRANSCRIPTS,
    TTN_GENE_INFO
)

logger = logging.getLogger(__name__)

# 各結構 band 的基因組範圍（hg38，負股：start > end），取自 data/cardiodb.xlsx 的 364 個外顯子；
# 相鄰 band 的分界放在中間的 intron
BAND_GENOMIC = {
    "Z-disk": (178807423, 178775280),
    "I-band": (178775279, 178618538),
    "A-band": (178618537, 178535916),
    "M-band": (178535915, 178525989),
}
# exon 49（LRG / meta-transcript 編號，hg38 178,741,921–178,739,141）：N2BA、N2B 含有，N2A 跳過
EXON49 = (178741921, 178739141)
EXON49_ROWS = ("N2BA", "N2B", "N2A")
EXON49_ALPHA = 0.35  # 心肌 isoform 上 exon 49 方框的不透明度
ISOFORM_GROUPS = (
    ("Reference (meta-transcript)", "#7f8c8d", ("Meta",)),
    ("Cardiac isoforms", "#e74c3c", ("N2BA", "N2B")),
    ("Skeletal muscle isoform", "#e67e22", ("N2A",)),
    ("Other isoforms (Novex)", "#95a5a6", ("Nvx1", "Nvx2", "Nvx3")),
)
ISOFORM_LABELS = {"Meta": "Meta", "Nvx1": "Novex-1", "Nvx2": "Novex-2", "Nvx3": "Novex-3"}


class ImageGenerator:
    """Generate protein domain schematics"""
    
    def __init__(self):
        self.output_dir = OUTPUT_DIR / "images"
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_titin_schematic(
        self,
        variant_info: Dict[str, str],
        xlsx_file: Optional[str] = None
    ) -> Path:
        """
        Generate titin protein schematic with variant location
        """
        logger.info("Generating titin protein schematic...")
        
        # Calculate variant position in protein coordinates
        genomic_pos = variant_info['pos']
        relative_pos = (genomic_pos - TTN_GENE_INFO['start']) / \
                      (TTN_GENE_INFO['end'] - TTN_GENE_INFO['start'])
        
        # Try to load transcript intervals if xlsx file provided or exists
        transcript_intervals_data = None
        if xlsx_file:
            try:
                df = pd.read_excel(xlsx_file)
                transcript_intervals_data = self._load_transcript_intervals(df)
                logger.info(f"Loaded transcript intervals from {xlsx_file}")
            except Exception as e:
                logger.warning(f"Failed to load transcript intervals: {e}")
        
        if transcript_intervals_data:
            fig = self._plot_isoform_panel(variant_info, transcript_intervals_data)
        else:
            fig = self._plot_simple_panel(variant_info, relative_pos)

        # Save figure with tight layout
        output_path = self.output_dir / f"titin_schematic_{variant_info['variant_id']}.png"
        fig.savefig(
            output_path,
            dpi=200,
            bbox_inches='tight',
            facecolor='white',
            pad_inches=0.15
        )
        plt.close(fig)
        
        logger.info(f"Schematic saved to: {output_path}")
        return output_path
    
    def _plot_simple_panel(self, variant_info: Dict, relative_pos: float):
        """沒有 transcript_interval.xlsx 時的簡易版（舊版畫法）。"""
        num_transcripts = len(TTN_TRANSCRIPTS)
        fig = plt.figure(figsize=(20, 5 + num_transcripts * 1.8))
        gs = fig.add_gridspec(num_transcripts + 1, 1, hspace=0.01, top=0.95, bottom=0.05, left=0.05, right=0.95)
        fig.suptitle(self._title(variant_info), fontsize=22, fontweight='bold', y=0.97)
        self._plot_domain_overview(fig.add_subplot(gs[0]), relative_pos)
        for idx, (name, info) in enumerate(TTN_TRANSCRIPTS.items()):
            self._plot_transcript(fig.add_subplot(gs[idx + 1]), name, info, relative_pos, variant_info)
        return fig

    @staticmethod
    def _title(variant_info: Dict) -> str:
        from utils import variant_ids

        d = variant_ids.add_display_fields(variant_info)
        return (
            f'TTN Variant Location: {d["display_id"]}\n'
            f'chr{variant_info["chrom"]}:{int(variant_info["pos"]):,} {d["plus_ref"]}>{d["plus_alt"]} (GRCh38 + strand)'
            + (f'  ·  {d["hgvs_c_display"]}' if d.get("hgvs_c_display") else '')
        )

    def _plot_isoform_panel(self, variant_info: Dict, intervals_by_tx: Dict[str, List[Tuple[int, int]]]):
        """一張圖、同一條基因組座標軸：上方是結構 band 總覽，下方是各 isoform 的外顯子區段。

        - 區段依所在 band（Z-disk / I-band / A-band / M-band）上色，與總覽同色
        - 心肌 isoform（N2BA、N2B）與骨骼肌 isoform（N2A）各用不同顏色的框框起來
        - 標出 exon 49（N2B 特有區段的外顯子；N2BA、N2B 含有，N2A 跳過）
        """
        pos = int(variant_info['pos'])
        x_var = self._normalize_genomic_position(pos)

        groups = []  # (標籤, 顏色, [transcript...])；沒有列在任何組的放到 Other
        listed = set()
        for label, color, names in ISOFORM_GROUPS:
            present = [n for n in names if n in intervals_by_tx]
            if present:
                groups.append((label, color, present))
                listed.update(present)
        rest = [n for n in intervals_by_tx if n not in listed]
        if rest:
            groups.append(('Other isoforms', '#95a5a6', rest))

        ROW = 0.62        # 每列間距（座標單位）
        BAR_H = 0.30
        GROUP_GAP = 0.55  # 組與組之間（含組名）
        LABEL_PAD = 0.24  # 組內第一列上方有 exon 49 文字時，外框往上多留的空間（文字不壓到框線）

        # exon 49 文字：同一組只寫一次（心肌組的 N2BA、N2B 共用一個「exon 49」），跳過的列寫「exon 49 skipped」
        ex49 = (EXON49[0], EXON49[1])

        def has_exon49(name: str) -> bool:
            return any(min(a, b) <= ex49[1] and ex49[0] <= max(a, b) for a, b in intervals_by_tx[name])

        label_rows = set()
        for _, _, names in groups:
            shown = False
            for name in names:
                if name not in EXON49_ROWS:
                    continue
                if not has_exon49(name):
                    label_rows.add(name)
                elif not shown:
                    label_rows.add(name)
                    shown = True
        padded_groups = sum(1 for g in groups if g[2][0] in label_rows)

        n_rows = sum(len(g[2]) for g in groups)
        total_h = 2.0 + n_rows * ROW + len(groups) * GROUP_GAP + padded_groups * LABEL_PAD + 0.3

        fig, ax = plt.subplots(figsize=(18, 1.6 + total_h * 1.05))
        fig.subplots_adjust(left=0.10, right=0.98, top=1 - 1.25 / (1.6 + total_h * 1.05), bottom=0.02)
        fig.suptitle(self._title(variant_info), fontsize=19, fontweight='bold', y=0.995)
        ax.set_xlim(-0.115, 1.035)
        ax.set_ylim(total_h, -0.15)
        ax.axis('off')

        # ── 總覽：band 長條（基因組比例）
        y = 0.35
        for band, (b_start, b_end) in BAND_GENOMIC.items():
            x0 = self._normalize_genomic_position(b_start)
            x1 = self._normalize_genomic_position(b_end)
            ax.add_patch(patches.Rectangle((x0, y - 0.2), x1 - x0, 0.4, facecolor=TTN_DOMAINS[band]['color'],
                                           edgecolor='black', linewidth=1.2, alpha=0.85))
            if x1 - x0 > 0.06:
                ax.text((x0 + x1) / 2, y, band, ha='center', va='center', fontsize=14, fontweight='bold')
        ax.text(-0.02, y, 'Band', ha='right', va='center', fontsize=15, fontweight='bold')
        ax.annotate('', xy=(x_var, y - 0.22), xytext=(x_var, y - 0.55),
                    arrowprops=dict(arrowstyle='-|>', color='red', linewidth=2.5, mutation_scale=22))
        # 圖例
        ly = 0.95
        for i, (band, info) in enumerate(TTN_DOMAINS.items()):
            lx = 0.0 + i * 0.25
            ax.add_patch(patches.Rectangle((lx, ly - 0.09), 0.025, 0.18, facecolor=info['color'],
                                           edgecolor='black', linewidth=1))
            ax.text(lx + 0.032, ly, f"{band}: exon {info['start']}–{info['end']}", va='center', fontsize=13)
        ax.text(0.0, 1.3, f'5′  chr2:{TTN_GENE_INFO["start"]:,}', ha='left', va='center', fontsize=11, color='#555')
        ax.text(0.5, 1.3, 'TTN is on the minus strand: drawn 5′ → 3′ in transcript order (genomic coordinates decrease left → right)',
                ha='center', va='center', fontsize=11, color='#555')
        ax.text(1.0, 1.3, f'chr2:{TTN_GENE_INFO["end"]:,}  3′', ha='right', va='center', fontsize=11, color='#555')

        # ── 各 isoform
        y = 2.0
        ex49_x0 = self._normalize_genomic_position(ex49[0])
        ex49_x1 = self._normalize_genomic_position(ex49[1])
        first_row_y = None
        for label, color, names in groups:
            pad = LABEL_PAD if names[0] in label_rows else 0.0
            y += GROUP_GAP + pad
            top = y - ROW / 2 + 0.02 - pad
            for name in names:
                y_row = y
                if first_row_y is None:
                    first_row_y = y_row
                ax.plot([0, 1], [y_row, y_row], color='#cccccc', linewidth=0.8, zorder=1)
                ax.text(-0.02, y_row, ISOFORM_LABELS.get(name, name), ha='right', va='center', fontsize=15,
                        fontweight='bold')
                contains = False
                has49 = False
                for a, b in intervals_by_tx[name]:
                    hi, lo = max(a, b), min(a, b)
                    contains |= lo <= pos <= hi
                    has49 |= lo <= ex49[1] and ex49[0] <= hi
                    self._draw_banded_interval(ax, hi, lo, y_row, BAR_H)
                if name in EXON49_ROWS:
                    # 文字標在方框正上方（y 軸是反向的：值越小越上面）。變異剛好落在 exon 49 時，
                    # 方框正上方是紅色箭頭，文字改放在箭頭左上方，避免重疊
                    label_y, label_va = y_row - BAR_H / 2 - 0.08, 'bottom'
                    if ex49[1] <= pos <= ex49[0]:
                        label_x, label_ha = ex49_x0 - 0.008, 'right'
                    else:
                        label_x, label_ha = (ex49_x0 + ex49_x1) / 2, 'center'
                    if has49:
                        # 半透明，看得到底下外顯子的 band 顏色
                        ax.add_patch(patches.Rectangle((ex49_x0, y_row - BAR_H / 2 - 0.05), ex49_x1 - ex49_x0,
                                                       BAR_H + 0.1, facecolor=to_rgba('#8e44ad', EXON49_ALPHA),
                                                       edgecolor='#4a235a', linewidth=1.5, zorder=6))
                        if name in label_rows:
                            ax.text(label_x, label_y, 'exon 49', ha=label_ha, va=label_va,
                                    fontsize=12, color='#4a235a', fontweight='bold', zorder=10)
                    else:
                        ax.add_patch(patches.Rectangle((ex49_x0, y_row - BAR_H / 2 - 0.05), ex49_x1 - ex49_x0,
                                                       BAR_H + 0.1, facecolor='none', edgecolor='#8e44ad',
                                                       linestyle='--', linewidth=1.5, zorder=6))
                        ax.text(label_x, label_y, 'exon 49 skipped', ha=label_ha,
                                va=label_va, fontsize=12, color='#8e44ad', fontstyle='italic', zorder=10)
                if contains:
                    ax.plot([x_var, x_var], [y_row - BAR_H / 2 - 0.1, y_row + BAR_H / 2 + 0.1], color='red',
                            linewidth=2.5, zorder=8)
                    ax.scatter(x_var, y_row - BAR_H / 2 - 0.14, marker='v', s=110, color='red', zorder=9,
                               edgecolors='black', linewidth=1.2)
                y += ROW
            bottom = y - ROW / 2 - 0.02
            ax.add_patch(FancyBboxPatch((-0.11, top), 1.14, bottom - top, boxstyle='round,pad=0,rounding_size=0.012',
                                        facecolor=to_rgba(color, 0.07), edgecolor=color, linewidth=2.2, zorder=0,
                                        mutation_aspect=10))
            ax.text(-0.105, top - 0.03, label, ha='left', va='bottom', fontsize=13, fontweight='bold', color=color,
                    zorder=10)
        # 變異位置貫穿所有列的參考線
        if first_row_y is not None:
            ax.plot([x_var, x_var], [0.55, y - ROW / 2], color='red', linewidth=1, linestyle=':', alpha=0.6,
                    zorder=2)
        return fig

    def _draw_banded_interval(self, ax, hi: int, lo: int, y: float, h: float) -> None:
        """把一段連續的外顯子區段依 band 切開上色。"""
        for band, (b_start, b_end) in BAND_GENOMIC.items():
            s_hi, s_lo = min(hi, b_start), max(lo, b_end)
            if s_hi < s_lo:
                continue
            x0 = self._normalize_genomic_position(s_hi)
            x1 = self._normalize_genomic_position(s_lo)
            ax.add_patch(patches.Rectangle((x0, y - h / 2), max(x1 - x0, 0.0015), h,
                                           facecolor=TTN_DOMAINS[band]['color'], edgecolor='black',
                                           linewidth=0.9, alpha=0.85, zorder=3))

    def _load_transcript_intervals(self, df: pd.DataFrame) -> Dict[str, List[Tuple[int, int]]]:
        """Load transcript intervals from DataFrame"""
        transcript_data = {}
        for _, row in df.iterrows():
            transcript_name = row['transcript']
            intervals = []
            for col in df.columns:
                if col.startswith('interval'):
                    interval = self._parse_interval(row[col])
                    if interval:
                        intervals.append(interval)
            if intervals:
                transcript_data[transcript_name] = intervals
        return transcript_data
    
    def _plot_domain_overview(self, ax, variant_pos):
        """Plot overview of TTN domains"""
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(0, 1.3)
        ax.axis('off')
        
        # Draw domains
        y_center = 0.5
        height = 0.22
        
        # Calculate domain positions (normalized)
        total_length = TTN_DOMAINS['M-band']['end'] - TTN_DOMAINS['Z-disk']['start'] + 1
        first_exon = TTN_DOMAINS['Z-disk']['start']
        
        for domain_name, domain_info in TTN_DOMAINS.items():
            start_norm = (domain_info['start'] - first_exon) / total_length
            end_norm = (domain_info['end'] - first_exon + 1) / total_length
            width = end_norm - start_norm
            
            # Draw domain rectangle without internal label
            rect = FancyBboxPatch(
                (start_norm, y_center - height/2),
                width,
                height,
                boxstyle="square,pad=0",
                edgecolor='black',
                facecolor=domain_info['color'],
                alpha=0.8,
                linewidth=1.5
            )
            ax.add_patch(rect)
            
            # 移除內部文字標籤
            # label_x = start_norm + width/2
            # ax.text(
            #     label_x,
            #     y_center,
            #     domain_name,
            #     ha='center',
            #     va='center',
            #     fontsize=13,
            #     fontweight='bold',
            #     color='white' if domain_name != 'I-band' else 'black'
            # )
        
        # Mark variant position (only line and arrow marker, no text label)
        marker_height = 0.18
        ax.plot(
            [variant_pos, variant_pos],
            [y_center - height/2 - marker_height, y_center + height/2 + marker_height],
            'r-',
            linewidth=3,
            zorder=5
        )
        ax.scatter(
            variant_pos, 
            y_center + height/2 + marker_height + 0.05, 
            marker='v', 
            s=250, 
            color='red', 
            zorder=10, 
            edgecolors='black', 
            linewidth=2
        )
        
        # Legend (增加顏色塊大小)
        legend_y = 0.08
        legend_x_start = 0.05
        legend_spacing = 0.24
        
        for idx, (domain_name, domain_info) in enumerate(TTN_DOMAINS.items()):
            legend_x = legend_x_start + (idx * legend_spacing)
            # 增加顏色塊大小（從 0.025x0.05 增加到 0.04x0.08）
            rect = patches.Rectangle(
                (legend_x, legend_y),
                0.04,  # 增加寬度從 0.025 到 0.04
                0.08,  # 增加高度從 0.05 到 0.08
                facecolor=domain_info['color'],
                edgecolor='black',
                linewidth=1.5  # 增加邊框粗細
            )
            ax.add_patch(rect)
            ax.text(
                legend_x + 0.045,  # 調整文字位置以配合新的顏色塊大小
                legend_y + 0.04,   # 垂直置中
                f"{domain_name}: Exon {domain_info['start']}-{domain_info['end']}",
                va='center',
                fontsize=20
            )
    
    def _plot_transcript(self, ax, name, info, variant_pos, variant_info):
        """Plot individual transcript (Simple bar fallback)"""
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.1, 1.0)
        ax.axis('off')
        
        ax.text(
            0.5,
            0.88,
            f'{name} - {info["description"]}\n'
            f'Length: {info["length"]} AA | Transcript: {info["id"]}',
            ha='center',
            va='top',
            fontsize=13,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', edgecolor='gray', alpha=0.8)
        )
        
        y_center = 0.42
        height = 0.15
        bar_start = 0.08
        bar_end = 0.92
        bar_width = bar_end - bar_start
        
        rect = FancyBboxPatch(
            (bar_start, y_center - height/2),
            bar_width,
            height,
            boxstyle="round,pad=0.01",
            edgecolor='black',
            facecolor='lightgray',
            alpha=0.3,
            linewidth=2
        )
        ax.add_patch(rect)
        
        # Colors based on domain proportions
        total_exons = TTN_DOMAINS['M-band']['end'] - TTN_DOMAINS['Z-disk']['start'] + 1
        first_exon = TTN_DOMAINS['Z-disk']['start']
        
        for domain_name, domain_info in TTN_DOMAINS.items():
            domain_start_norm = (domain_info['start'] - first_exon) / total_exons
            domain_end_norm = (domain_info['end'] - first_exon + 1) / total_exons
            
            start_x = bar_start + domain_start_norm * bar_width
            end_x = bar_start + domain_end_norm * bar_width
            width = end_x - start_x
            
            if width > 0.001:
                rect = patches.Rectangle(
                    (start_x, y_center - height/2),
                    width,
                    height,
                    edgecolor='black',
                    facecolor=domain_info['color'],
                    alpha=0.7,
                    linewidth=0.8
                )
                ax.add_patch(rect)
        
        # Mark variant position
        variant_x = bar_start + variant_pos * bar_width
        marker_height = 0.12
        
        ax.plot(
            [variant_x, variant_x],
            [y_center - height/2 - marker_height, y_center + height/2 + marker_height],
            'r-',
            linewidth=2.5,
            zorder=5
        )
        ax.scatter(
            variant_x,
            y_center + height/2 + marker_height + 0.03,
            marker='v',
            s=150,
            color='red',
            zorder=10,
            edgecolors='black',
            linewidth=1.5
        )
        
        # Labels
        label_y = y_center - height/2 - marker_height - 0.05
        ax.text(
            bar_start, label_y, '1', ha='center', va='top', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='black', alpha=0.9)
        )
        ax.text(
            bar_end, label_y, str(info['length']), ha='center', va='top', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='black', alpha=0.9)
        )
    
    def _parse_interval(self, interval_str) -> Optional[Tuple[int, int]]:
        """Parse interval string '178807423-178758984' to (start, end)"""
        if pd.isna(interval_str) or interval_str == '':
            return None
        parts = str(interval_str).strip().split('-')
        if len(parts) != 2:
            return None
        try:
            return (int(parts[0]), int(parts[1]))
        except:
            return None
    
    def _normalize_genomic_position(self, pos: int) -> float:
        """Normalize genomic position to 0-1 range"""
        gene_length = TTN_GENE_INFO['start'] - TTN_GENE_INFO['end']
        return (TTN_GENE_INFO['start'] - pos) / gene_length