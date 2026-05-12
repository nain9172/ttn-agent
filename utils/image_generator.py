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
        
        # Determine number of transcripts to plot
        if transcript_intervals_data:
            num_transcripts = len(transcript_intervals_data)
        else:
            num_transcripts = len(TTN_TRANSCRIPTS)
        
        # Create figure with proper spacing
        fig_height = 5 + (num_transcripts * 1.8)
        fig = plt.figure(figsize=(20, fig_height))
        
        # Adjust spacing to prevent overlap (reduced hspace for tighter layout)
        gs = fig.add_gridspec(
            num_transcripts + 1, 
            1, 
            hspace=0.1,  # 減少間距從 0.8 到 0.3
            top=0.95,
            bottom=0.05,
            left=0.05,
            right=0.95
        )
        
        # Main title
        fig.suptitle(
            f'TTN Variant Location: {variant_info["variant_id"]}\n'
            f'Genomic Position: chr{variant_info["chrom"]}:{variant_info["pos"]} '
            f'({variant_info["ref"]}>{variant_info["alt"]})',
            fontsize=22,
            fontweight='bold',
            y=0.97
        )
        
        # Plot domain overview (Top row - Protein Structure)
        ax_main = fig.add_subplot(gs[0])
        self._plot_domain_overview(ax_main, relative_pos)
        
        # Plot transcripts
        if transcript_intervals_data:
            for idx, (transcript_name, intervals) in enumerate(transcript_intervals_data.items()):
                ax = fig.add_subplot(gs[idx + 1])
                self._plot_transcript_intervals_aligned(
                    ax,
                    transcript_name,
                    intervals,
                    genomic_pos
                )
        else:
            # Fallback if no interval file is provided
            for idx, (transcript_name, transcript_info) in enumerate(TTN_TRANSCRIPTS.items()):
                ax = fig.add_subplot(gs[idx + 1])
                self._plot_transcript(
                    ax,
                    transcript_name,
                    transcript_info,
                    relative_pos,
                    variant_info
                )
        
        # Save figure with tight layout
        output_path = self.output_dir / f"titin_schematic_{variant_info['variant_id']}.png"
        plt.savefig(
            output_path, 
            dpi=300, 
            bbox_inches='tight',
            facecolor='white',
            pad_inches=0.2
        )
        plt.close()
        
        logger.info(f"Schematic saved to: {output_path}")
        return output_path
    
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
    
    def _plot_transcript_intervals_aligned(
        self,
        ax,
        transcript_name: str,
        intervals: List[Tuple[int, int]],
        variant_pos: int
    ):
        """
        Plot transcript intervals aligned with domain overview.
        Checks if variant_pos lies within any interval before plotting the marker.
        """
        # Colors for intervals
        interval_colors = [
            '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4',
            '#FECA57', '#FF9FF3', '#54A0FF', '#48DBFB',
            '#00D2D3'
        ]
        
        # Use same x-axis range as domain overview (0-1 normalized)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(0, 1.0)
        ax.axis('off')
        
        # Transcript name label (Left side)
        ax.text(
            -0.03,
            0.5,
            f'{transcript_name}',
            fontsize=20,
            fontweight='bold',
            va='center',
            ha='right',
            bbox=dict(
                boxstyle='round,pad=0.4',
                facecolor='lightblue',
                edgecolor='navy',
                alpha=0.8,
                linewidth=1.5
            )
        )
        
        # Draw genomic coordinate reference
        y_pos = 0.5
        height = 0.15
        
        # Background reference line
        ax.plot([0, 1], [y_pos, y_pos], 'k-', linewidth=0.8, alpha=0.2, zorder=1)
        
        # Check if variant is inside any interval of this transcript
        variant_in_transcript = False
        
        # Draw each interval
        for i, (start, end) in enumerate(intervals):
            # Normalize positions using same method as domain overview
            start_norm = self._normalize_genomic_position(start)
            end_norm = self._normalize_genomic_position(end)
            
            # Check overlap logic: intervals can be start>end or end>start depending on strand notation
            # We use min/max to be safe
            interval_min = min(start, end)
            interval_max = max(start, end)
            
            if interval_min <= variant_pos <= interval_max:
                variant_in_transcript = True
            
            # Calculate width
            width = end_norm - start_norm
            
            # Use different colors
            color = interval_colors[i % len(interval_colors)]
            
            # Draw interval rectangle
            rect = FancyBboxPatch(
                (start_norm, y_pos - height/2),
                width,
                height,
                boxstyle="round,pad=0.002",
                edgecolor='black',
                facecolor=color,
                alpha=0.8,
                linewidth=1.2,
                zorder=3
            )
            ax.add_patch(rect)
        
        # Mark variant position ONLY if it falls within an interval
        if variant_in_transcript:
            variant_norm = self._normalize_genomic_position(variant_pos)
            marker_height = 0.12
            
            ax.plot(
                [variant_norm, variant_norm],
                [y_pos - height/2 - marker_height, y_pos + height/2 + marker_height],
                'r-',
                linewidth=2.5,
                zorder=5
            )
            ax.scatter(
                variant_norm,
                y_pos + height/2 + marker_height + 0.03,
                marker='v',
                s=120,
                color='red',
                zorder=10,
                edgecolors='black',
                linewidth=1.5
            )
    
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