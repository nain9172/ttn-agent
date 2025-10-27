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
    REFERENCE_GENOME_PATH,
    TTN_SEQUENCE_START,
    TTN_SEQUENCE_END
)

logger = logging.getLogger(__name__)


class Evo2Predictor:
    """Evo2-based variant pathogenicity predictor"""
    
    def __init__(self):
        self.model = None
        self.seq_ttn = None  # TTN reference sequence
        self.window_size = EVO2_WINDOW_SIZE
        self.ttn_start = TTN_SEQUENCE_START
        self.ttn_end = TTN_SEQUENCE_END
        
    def _load_model(self):
        """Load Evo2 model (lazy loading)"""
        if self.model is not None:
            return
        
        from evo2.models import Evo2
        logger.info(f"Loading Evo2 model: {EVO2_MODEL}")
        self.model = Evo2(EVO2_MODEL)
        logger.info("Evo2 model loaded successfully")
    
    def _load_reference_sequence(self):
        """Load TTN reference sequence"""
        if self.seq_ttn is not None:
            return
        
        if not REFERENCE_GENOME_PATH.exists():
            raise FileNotFoundError(
                f"Reference sequence not found at {REFERENCE_GENOME_PATH}. "
                "Please ensure sequence.fasta is present in the data directory."
            )
        
        try:
            from Bio import SeqIO
            logger.info("Loading TTN reference sequence...")
            
            with open(REFERENCE_GENOME_PATH, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    self.seq_ttn = str(record.seq)
                    logger.info(
                        f"TTN reference sequence loaded: {len(self.seq_ttn)} bp "
                        f"(chr2:{self.ttn_start}-{self.ttn_end})"
                    )
                    break
            
            if self.seq_ttn is None:
                raise ValueError("TTN sequence not found in reference file")
            
            # Verify sequence length matches expected TTN region
            expected_length = self.ttn_start - self.ttn_end + 1
            if len(self.seq_ttn) != expected_length:
                logger.warning(
                    f"Sequence length mismatch: got {len(self.seq_ttn)} bp, "
                    f"expected {expected_length} bp"
                )
                
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
        Parse reference and variant sequences from TTN sequence
        
        Args:
            pos: Genomic position (1-based, chr2 coordinates)
            ref: Reference base
            alt: Alternate base
        
        Returns:
            Tuple of (ref_seq, var_seq) or (None, None) if invalid
        """
        # Check if position is within TTN region
        if pos > self.ttn_start or pos < self.ttn_end:
            logger.error(
                f"Position {pos} is outside TTN region "
                f"({self.ttn_start}-{self.ttn_end})"
            )
            return None, None
        
        # Convert genomic position to sequence index
        # TTN is on negative strand: position 178807423 = index 0
        # sequence.fasta is already reverse complemented (indicated by 'c' in header)
        seq_index = self.ttn_start - pos
        
        if seq_index < 0 or seq_index >= len(self.seq_ttn):
            logger.error(
                f"Calculated sequence index {seq_index} out of bounds "
                f"(sequence length: {len(self.seq_ttn)})"
            )
            return None, None
        
        # Extract window around the variant position
        window_half = self.window_size // 2
        ref_seq_start = max(0, seq_index - window_half)
        ref_seq_end = min(len(self.seq_ttn), seq_index + window_half)
        ref_seq = self.seq_ttn[ref_seq_start:ref_seq_end]
        
        snv_pos_in_ref = seq_index - ref_seq_start
        
        if snv_pos_in_ref < 0 or snv_pos_in_ref >= len(ref_seq):
            logger.error(f"SNV position out of bounds in extracted sequence")
            return None, None
        
        # Validate reference base
        if ref_seq[snv_pos_in_ref].upper() != ref.upper():
            logger.warning(
                f"Reference mismatch at position {pos} (index {seq_index}): "
                f"expected {ref}, got {ref_seq[snv_pos_in_ref]}"
            )
            return None, None
        
        # Create variant sequence
        var_seq = ref_seq[:snv_pos_in_ref] + alt + ref_seq[snv_pos_in_ref + 1:]
        
        if len(var_seq) != len(ref_seq):
            logger.error("Variant and reference sequences have different lengths")
            return None, None
        
        logger.info(
            f"Extracted sequence window: genomic pos {pos} -> "
            f"seq index {seq_index}, window size {len(ref_seq)} bp"
        )
        
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