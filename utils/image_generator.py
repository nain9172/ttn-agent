"""
Image Generator Module - FIXED VERSION
Generates protein schematic images showing variant location
Fixed: overlap and clipping issues
"""

import logging
from pathlib import Path
from typing import Dict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np

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
        variant_info: Dict[str, str]
    ) -> Path:
        """
        Generate titin protein schematic with variant location
        
        Args:
            variant_info: Variant information dictionary
        
        Returns:
            Path to generated image
        """
        logger.info("Generating titin protein schematic...")
        
        # Calculate variant position in protein coordinates
        genomic_pos = variant_info['pos']
        relative_pos = (genomic_pos - TTN_GENE_INFO['start']) / \
                      (TTN_GENE_INFO['end'] - TTN_GENE_INFO['start'])
        
        # Create figure with proper spacing - FIXED
        num_transcripts = len(TTN_TRANSCRIPTS)
        fig_height = 5 + (num_transcripts * 3)  # Much larger spacing
        fig = plt.figure(figsize=(18, fig_height))
        
        # Adjust spacing to prevent overlap - FIXED
        gs = fig.add_gridspec(
            num_transcripts + 1, 
            1, 
            hspace=1.2,  # Much larger vertical spacing to prevent overlap
            top=0.92,    # Add margins
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
            y=0.96  # Position at top with more space
        )
        
        # Plot domain overview
        ax_main = fig.add_subplot(gs[0])
        self._plot_domain_overview(ax_main, relative_pos)
        
        # Plot each transcript with better spacing
        for idx, (transcript_name, transcript_info) in enumerate(TTN_TRANSCRIPTS.items()):
            ax = fig.add_subplot(gs[idx + 1])
            self._plot_transcript(
                ax,
                transcript_name,
                transcript_info,
                relative_pos,
                variant_info
            )
        
        # Save figure with tight layout - FIXED
        output_path = self.output_dir / f"titin_schematic_{variant_info['variant_id']}.png"
        plt.savefig(
            output_path, 
            dpi=300, 
            bbox_inches='tight',  # Prevent clipping
            facecolor='white',
            pad_inches=0.2  # Add padding
        )
        plt.close()
        
        logger.info(f"Schematic saved to: {output_path}")
        return output_path
    
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