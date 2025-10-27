"""
Image Generator Module - FIXED VERSION
Generates protein schematic images showing variant location
Fixed: overlap and clipping issues
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
        self.output_dir.mkdir(exist_ok=True)
    
    def generate_titin_schematic(
        self,
        variant_info: Dict[str, str],
        xlsx_file: Optional[str] = None
    ) -> Path:
        """
        Generate titin protein schematic with variant location
        
        Args:
            variant_info: Variant information dictionary
            xlsx_file: Optional path to transcript intervals Excel file
        
        Returns:
            Path to generated image
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
        
        # Adjust spacing to prevent overlap
        gs = fig.add_gridspec(
            num_transcripts + 1, 
            1, 
            hspace=0.8,
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
            fontsize=16,
            fontweight='bold',
            y=0.97
        )
        
        # Plot domain overview (保留 Z-disk, I-band 等)
        ax_main = fig.add_subplot(gs[0])
        self._plot_domain_overview(ax_main, relative_pos)
        
        # Plot transcripts: use intervals if available, otherwise use original method
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
        """Plot transcript intervals aligned with domain overview"""
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
        
        # Transcript name label (左側)
        ax.text(
            -0.03,
            0.5,
            f'{transcript_name}',
            fontsize=11,
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
        
        # Draw each interval
        for i, (start, end) in enumerate(intervals):
            # Normalize positions using same method as domain overview
            start_norm = self._normalize_genomic_position(start)
            end_norm = self._normalize_genomic_position(end)
            
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
        
        # Mark variant position
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
        """Plot overview of TTN domains - FIXED"""
        ax.set_xlim(-0.05, 1.05)  # Extended limits
        ax.set_ylim(0, 1.3)  # Extended to give more space for labels
        ax.axis('off')
        # ax.set_title(
        #     'TTN Protein Domain Structure (Exon-based)', 
        #     fontsize=13, 
        #     fontweight='bold', 
        #     pad=25
        # )
        
        # Draw domains
        y_center = 0.5
        height = 0.22  # Reduced height for better spacing
        
        # Calculate domain positions (normalized) - use exon 1 as start
        total_length = TTN_DOMAINS['M-band']['end'] - TTN_DOMAINS['Z-disk']['start'] + 1
        first_exon = TTN_DOMAINS['Z-disk']['start']
        
        for domain_name, domain_info in TTN_DOMAINS.items():
            # Normalize from first exon (1) to last exon (364)
            start_norm = (domain_info['start'] - first_exon) / total_length
            end_norm = (domain_info['end'] - first_exon + 1) / total_length  # +1 to include end exon
            width = end_norm - start_norm
            
            # Draw domain rectangle with no overlap
            rect = FancyBboxPatch(
                (start_norm, y_center - height/2),
                width,
                height,
                boxstyle="square,pad=0",  # Use square to prevent gaps
                edgecolor='black',
                facecolor=domain_info['color'],
                alpha=0.8,
                linewidth=1.5
            )
            ax.add_patch(rect)
            
            # Add domain label
            label_x = start_norm + width/2
            ax.text(
                label_x,
                y_center,
                domain_name,
                ha='center',
                va='center',
                fontsize=10,
                fontweight='bold',
                color='white' if domain_name != 'I-band' else 'black'
            )
        
        # Mark variant position - IMPROVED
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
        
        # Add variant annotation - REPOSITIONED
        ax.text(
            variant_pos,
            y_center + height/2 + marker_height + 0.18,
            'Variant\nLocation',
            ha='center',
            va='bottom',
            fontsize=9,
            fontweight='bold',
            color='red',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='red', alpha=0.8)
        )
        
        # Add legend for domains - REPOSITIONED
        legend_y = 0.08  # Higher position to avoid overlap
        legend_x_start = 0.05
        legend_spacing = 0.24
        
        for idx, (domain_name, domain_info) in enumerate(TTN_DOMAINS.items()):
            legend_x = legend_x_start + (idx * legend_spacing)
            rect = patches.Rectangle(
                (legend_x, legend_y),
                0.025,
                0.05,
                facecolor=domain_info['color'],
                edgecolor='black',
                linewidth=1
            )
            ax.add_patch(rect)
            ax.text(
                legend_x + 0.03,
                legend_y + 0.025,
                f"{domain_name}: Exon {domain_info['start']}-{domain_info['end']}",
                va='center',
                fontsize=7.5
            )
    
    def _plot_transcript(self, ax, name, info, variant_pos, variant_info):
        """Plot individual transcript with variant location - FIXED"""
        ax.set_xlim(-0.05, 1.05)  # Extended limits
        ax.set_ylim(-0.1, 1.0)  # Adjusted range for better spacing
        ax.axis('off')
        
        # Transcript title - REPOSITIONED with more space
        ax.text(
            0.5,
            0.88,  # Adjusted position
            f'{name} - {info["description"]}\n'
            f'Length: {info["length"]} AA | Transcript: {info["id"]}',
            ha='center',
            va='top',
            fontsize=9.5,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', edgecolor='gray', alpha=0.8)
        )
        
        # Draw transcript bar
        y_center = 0.42  # Adjusted for better spacing
        height = 0.15
        bar_start = 0.08
        bar_end = 0.92
        bar_width = bar_end - bar_start
        
        # Transcript background
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
        
        # Get transcript length for calculations
        total_length = info['length']
        
        # Color-code domains proportionally (based on exon ranges)
        total_exons = TTN_DOMAINS['M-band']['end'] - TTN_DOMAINS['Z-disk']['start'] + 1
        first_exon = TTN_DOMAINS['Z-disk']['start']
        
        for domain_name, domain_info in TTN_DOMAINS.items():
            # Calculate normalized domain position based on exon ranges
            domain_start_norm = (domain_info['start'] - first_exon) / total_exons
            domain_end_norm = (domain_info['end'] - first_exon + 1) / total_exons  # +1 to include end exon
            
            # Map to bar coordinates
            start_x = bar_start + domain_start_norm * bar_width
            end_x = bar_start + domain_end_norm * bar_width
            width = end_x - start_x
            
            if width > 0.001:  # Only draw if visible
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
        
        # Add position labels - positioned below the bar
        label_y = y_center - height/2 - marker_height - 0.05
        
        ax.text(
            bar_start,
            label_y,
            '1',
            ha='center',
            va='top',
            fontsize=7.5,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='black', alpha=0.9)
        )
        ax.text(
            bar_end,
            label_y,
            str(total_length),
            ha='center',
            va='top',
            fontsize=7.5,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='black', alpha=0.9)
        )
        
        # Add variant position estimate
        estimated_aa = int(variant_pos * total_length)
        ax.text(
            variant_x,
            label_y,
            f'~{estimated_aa}',
            ha='center',
            va='top',
            fontsize=7.5,
            color='red',
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='red', alpha=0.95)
        )
    
    def generate_transcript_intervals(
        self,
        xlsx_file: str,
        variant_info: Optional[Dict[str, str]] = None
    ) -> Path:
        """
        Generate transcript intervals diagram from Excel file
        
        Args:
            xlsx_file: Path to Excel file containing transcript intervals
            variant_info: Optional variant information to mark on the diagram
        
        Returns:
            Path to generated image
        """
        logger.info("Generating transcript intervals diagram...")
        
        # Read data
        df = pd.read_excel(xlsx_file)
        logger.info(f"Loaded {len(df)} transcripts from {xlsx_file}")
        
        # Create figure
        num_transcripts = len(df)
        fig_height = 3 + (num_transcripts * 1.5)
        fig = plt.figure(figsize=(20, fig_height))
        
        # Setup grid layout
        gs = fig.add_gridspec(
            num_transcripts,
            1,
            hspace=0.6,
            top=0.95,
            bottom=0.05,
            left=0.08,
            right=0.95
        )
        
        # Main title
        title = 'TTN Transcript Intervals'
        if variant_info:
            title += f' - Variant: {variant_info["variant_id"]}'
        
        fig.suptitle(
            title,
            fontsize=18,
            fontweight='bold',
            y=0.98
        )
        
        # Calculate variant position if provided
        variant_pos = None
        if variant_info:
            genomic_pos = variant_info['pos']
            variant_pos = genomic_pos
        
        # Plot each transcript
        for idx, row in df.iterrows():
            ax = fig.add_subplot(gs[idx])
            transcript_name = row['transcript']
            
            # Parse all intervals
            intervals = []
            for col in df.columns:
                if col.startswith('interval'):
                    interval = self._parse_interval(row[col])
                    if interval:
                        intervals.append(interval)
            
            logger.info(f"{transcript_name}: {len(intervals)} intervals")
            
            # Plot this transcript
            self._plot_transcript_with_intervals(
                ax,
                transcript_name,
                intervals,
                variant_pos
            )
        
        # Save figure
        output_filename = 'transcript_intervals.png'
        if variant_info:
            output_filename = f"transcript_intervals_{variant_info['variant_id']}.png"
        
        output_path = self.output_dir / output_filename
        plt.savefig(
            output_path,
            dpi=300,
            bbox_inches='tight',
            facecolor='white',
            pad_inches=0.3
        )
        plt.close()
        
        logger.info(f"Transcript intervals diagram saved to: {output_path}")
        return output_path
    
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
    
    def _plot_transcript_with_intervals(
        self,
        ax,
        transcript_name: str,
        intervals: List[Tuple[int, int]],
        variant_pos: Optional[int] = None
    ):
        """Plot a single transcript with its intervals"""
        # Colors for intervals
        interval_colors = [
            '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4',
            '#FECA57', '#FF9FF3', '#54A0FF', '#48DBFB',
            '#00D2D3'
        ]
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        
        # Transcript name
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
        
        # Draw genomic coordinate axis
        y_pos = 0.45
        height = 0.12
        
        # Background reference line
        ax.plot([0.12, 0.95], [y_pos, y_pos], 'k-', linewidth=1, alpha=0.3, zorder=1)
        
        # Draw each interval
        for i, (start, end) in enumerate(intervals):
            # Normalize positions (TTN is on negative strand)
            start_norm = self._normalize_genomic_position(start)
            end_norm = self._normalize_genomic_position(end)
            
            # Map to x coordinates
            x_start = 0.12 + start_norm * 0.83
            x_end = 0.12 + end_norm * 0.83
            width = x_end - x_start
            
            if width > 0:
                # Use different colors
                color = interval_colors[i % len(interval_colors)]
                
                # Draw interval rectangle
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
                
                # Add interval label if wide enough
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
        
        # Mark variant position if provided
        if variant_pos is not None:
            variant_norm = self._normalize_genomic_position(variant_pos)
            variant_x = 0.12 + variant_norm * 0.83
            marker_height = 0.15
            
            ax.plot(
                [variant_x, variant_x],
                [y_pos - height/2 - marker_height, y_pos + height/2 + marker_height],
                'r-',
                linewidth=2.5,
                zorder=5
            )
            ax.scatter(
                variant_x,
                y_pos + height/2 + marker_height + 0.03,
                marker='v',
                s=150,
                color='red',
                zorder=10,
                edgecolors='black',
                linewidth=1.5
            )
            
            # Add variant label
            ax.text(
                variant_x,
                y_pos + height/2 + marker_height + 0.08,
                'Variant',
                ha='center',
                va='bottom',
                fontsize=8,
                fontweight='bold',
                color='red',
                bbox=dict(
                    boxstyle='round,pad=0.3',
                    facecolor='white',
                    edgecolor='red',
                    alpha=0.9
                )
            )
        
        # Add coordinate axis markers
        # Start position
        ax.text(
            0.12,
            y_pos - height/2 - 0.12,
            f'{TTN_GENE_INFO["start"]:,}',
            ha='center',
            va='top',
            fontsize=8,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.9)
        )
        
        # End position
        ax.text(
            0.95,
            y_pos - height/2 - 0.12,
            f'{TTN_GENE_INFO["end"]:,}',
            ha='center',
            va='top',
            fontsize=8,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.9)
        )
        
        # Middle reference point
        mid_pos = (TTN_GENE_INFO['start'] + TTN_GENE_INFO['end']) // 2
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
        
        # Add interval count info
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