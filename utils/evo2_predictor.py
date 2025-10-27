"""
Evo2 Predictor Module
Uses Evo2 model to predict variant pathogenicity
Based on the notebook: TTN_predict.ipynb
"""

import logging
from typing import Dict, Optional, Tuple
from pathlib import Path

from config import (
    EVO2_MODEL,
    EVO2_WINDOW_SIZE,
    PATHOGENIC_THRESHOLD,
    REFERENCE_GENOME_PATH
)

logger = logging.getLogger(__name__)


class Evo2Predictor:
    """Evo2-based variant pathogenicity predictor"""
    
    def __init__(self):
        self.model = None
        self.seq_chr2 = None
        self.window_size = EVO2_WINDOW_SIZE
        
    def _load_model(self):
        """Load Evo2 model (lazy loading)"""
        if self.model is not None:
            return
        
        from evo2.models import Evo2
        logger.info(f"Loading Evo2 model: {EVO2_MODEL}")
        self.model = Evo2(EVO2_MODEL)
        logger.info("Evo2 model loaded successfully")
    
    def _load_reference_sequence(self):
        """Load chromosome 2 reference sequence"""
        if self.seq_chr2 is not None:
            return
        
        if not REFERENCE_GENOME_PATH.exists():
            raise FileNotFoundError(
                f"Reference genome not found at {REFERENCE_GENOME_PATH}. "
                "Please download GRCh38 reference genome from:\n"
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
                "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
            )
        
        try:
            from Bio import SeqIO
            logger.info("Loading reference sequence for chromosome 2...")
            
            with open(REFERENCE_GENOME_PATH, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    header = record.description.lower()
                    if 'chromosome 2' in header or record.id.startswith('NC_000002'):
                        self.seq_chr2 = str(record.seq)
                        logger.info(f"Reference sequence loaded: {len(self.seq_chr2)} bp")
                        break
            
            if self.seq_chr2 is None:
                raise ValueError("Chromosome 2 sequence not found in reference genome")
                
        except Exception as e:
            logger.error(f"Failed to load reference sequence: {e}")
            raise
    
    def _parse_sequences(
        self,
        pos: int,
        ref: str,
        alt: str
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Parse reference and variant sequences from genome
        
        Args:
            pos: Genomic position (1-based)
            ref: Reference base
            alt: Alternate base
        
        Returns:
            Tuple of (ref_seq, var_seq) or (None, None) if invalid
        """
        p = pos - 1  # Convert to 0-indexed
        
        ref_seq_start = max(0, p - self.window_size // 2)
        ref_seq_end = min(len(self.seq_chr2), p + self.window_size // 2)
        ref_seq = self.seq_chr2[ref_seq_start:ref_seq_end]
        
        snv_pos_in_ref = p - ref_seq_start
        
        if snv_pos_in_ref < 0 or snv_pos_in_ref >= len(ref_seq):
            logger.error(f"SNV position {p} out of bounds")
            return None, None
        
        # Validate reference base
        if ref_seq[snv_pos_in_ref] != ref:
            logger.warning(
                f"Reference mismatch at position {pos}: "
                f"expected {ref}, got {ref_seq[snv_pos_in_ref]}"
            )
            return None, None
        
        # Create variant sequence
        var_seq = ref_seq[:snv_pos_in_ref] + alt + ref_seq[snv_pos_in_ref + 1:]
        
        if len(var_seq) != len(ref_seq):
            logger.error("Variant and reference sequences have different lengths")
            return None, None
        
        return ref_seq, var_seq
    
    def predict(self, variant_info: Dict[str, str]) -> Dict:
        """
        Predict variant pathogenicity using Evo2
        
        Args:
            variant_info: Dictionary with variant information
        
        Returns:
            Dictionary with prediction results
        """
        logger.info(f"Predicting pathogenicity for {variant_info['variant_id']}")
        
        try:
            # Load model and reference sequence
            self._load_model()
            self._load_reference_sequence()
            
            # Parse sequences
            ref_seq, var_seq = self._parse_sequences(
                variant_info['pos'],
                variant_info['ref'],
                variant_info['alt']
            )
            
            if ref_seq is None or var_seq is None:
                return {
                    'success': False,
                    'error': 'Reference sequence mismatch or parsing error'
                }
            
            # Score sequences
            logger.info("Scoring reference sequence...")
            ref_scores = self.model.score_sequences([ref_seq])
            ref_score = float(ref_scores[0])
            
            logger.info("Scoring variant sequence...")
            var_scores = self.model.score_sequences([var_seq])
            var_score = float(var_scores[0])
            
            # Calculate delta score
            delta_score = var_score - ref_score
            
            # Determine prediction
            prediction = "pathogenic" if delta_score < PATHOGENIC_THRESHOLD else "benign"
            
            result = {
                'success': True,
                'ref_score': ref_score,
                'var_score': var_score,
                'delta_score': delta_score,
                'prediction': prediction,
                'confidence': abs(delta_score),
                'threshold': PATHOGENIC_THRESHOLD
            }
            
            logger.info(f"Prediction: {prediction} (delta_score: {delta_score:.6f})")
            return result
            
        except Exception as e:
            logger.error(f"Error during prediction: {e}", exc_info=True)
            return {
                'success': False,
                'error': str(e)
            }