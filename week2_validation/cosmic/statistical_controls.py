"""
COSMIC Statistical Controls - Negative Control Validation
Week 2: Data Integrity & Statistical Validation

This module provides negative control tests to validate that observed COSMIC
agreement is not due to chance. The negative control shuffles recurrence counts
and recomputes Spearman correlation to verify that real agreement is significantly
stronger than random.

WHAT DOES THIS MODULE DO?
This module implements a negative control test by:
1. Shuffling recurrence counts in COSMIC data
2. Recomputing Spearman correlation with shuffled data
3. Comparing shuffled correlation to real correlation

If real correlation >> shuffled correlation, this provides evidence that
the observed agreement is biologically meaningful, not random.

DESIGN PRINCIPLE:
Negative controls verify that statistical signals are real, not artifacts
of the analysis method or random chance.
"""

import logging
from typing import Dict, Optional, Tuple

import numpy as np

_logger = logging.getLogger(__name__)

# Try to import scipy for Spearman correlation
try:
    from scipy.stats import spearmanr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    _logger.warning("scipy not available; negative control test will not be computed")


def compute_negative_control_correlation(
    fusion_pairs: Dict[Tuple[str, str], float],
    cosmic_pairs: Dict[Tuple[str, str], float],
    overlap: set,
    n_shuffles: int = 1000,
    random_seed: int = 42,
) -> Dict[str, Optional[float]]:
    """
    Compute negative control Spearman correlation using shuffled COSMIC recurrence counts.
    
    This function shuffles the recurrence counts in COSMIC data (breaking the
    biological relationship) and recomputes Spearman correlation. If the real
    correlation is much stronger than shuffled correlations, this provides evidence
    that the observed agreement is biologically meaningful.
    
    Args:
        fusion_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        cosmic_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        overlap: Set of overlapping fusion pairs
        n_shuffles: Number of shuffle iterations (default: 1000)
        random_seed: Random seed for reproducibility (default: 42)
    
    Returns:
        Dictionary with:
        - negative_control_rho: Mean Spearman rho from shuffled data (or None)
        - negative_control_p_value: p-value comparing real vs shuffled (or None)
        - shuffled_rho_mean: Mean of shuffled correlations
        - shuffled_rho_std: Standard deviation of shuffled correlations
    """
    if not SCIPY_AVAILABLE:
        return {
            "negative_control_rho": None,
            "negative_control_p_value": None,
            "shuffled_rho_mean": None,
            "shuffled_rho_std": None,
        }
    
    if len(overlap) < 3:
        # Need at least 3 pairs for meaningful correlation
        return {
            "negative_control_rho": None,
            "negative_control_p_value": None,
            "shuffled_rho_mean": None,
            "shuffled_rho_std": None,
        }
    
    # Extract recurrence counts for overlapping pairs (real data)
    fusion_counts = []
    cosmic_counts = []
    overlap_list = list(overlap)
    
    for pair in overlap_list:
        fusion_counts.append(fusion_pairs[pair])
        cosmic_counts.append(cosmic_pairs[pair])
    
    fusion_counts = np.array(fusion_counts)
    cosmic_counts = np.array(cosmic_counts)
    
    # Compute real Spearman correlation
    try:
        real_rho, _ = spearmanr(fusion_counts, cosmic_counts)
        if np.isnan(real_rho):
            return {
                "negative_control_rho": None,
                "negative_control_p_value": None,
                "shuffled_rho_mean": None,
                "shuffled_rho_std": None,
            }
    except Exception as e:
        _logger.warning(f"Error computing real Spearman correlation: {e}")
        return {
            "negative_control_rho": None,
            "negative_control_p_value": None,
            "shuffled_rho_mean": None,
            "shuffled_rho_std": None,
        }
    
    # Shuffle COSMIC counts and recompute correlation
    rng = np.random.default_rng(random_seed)
    shuffled_rhos = []
    
    for _ in range(n_shuffles):
        # Shuffle COSMIC counts (preserves distribution, breaks biological relationship)
        shuffled_cosmic = rng.permutation(cosmic_counts)
        
        try:
            rho_shuffled, _ = spearmanr(fusion_counts, shuffled_cosmic)
            if not np.isnan(rho_shuffled):
                shuffled_rhos.append(rho_shuffled)
        except Exception:
            continue
    
    if len(shuffled_rhos) == 0:
        return {
            "negative_control_rho": None,
            "negative_control_p_value": None,
            "shuffled_rho_mean": None,
            "shuffled_rho_std": None,
        }
    
    shuffled_rhos = np.array(shuffled_rhos)
    shuffled_mean = float(np.mean(shuffled_rhos))
    shuffled_std = float(np.std(shuffled_rhos))
    
    # Compute p-value: proportion of shuffled correlations >= real correlation
    # This tests if real correlation is significantly higher than random
    if shuffled_std > 0:
        # Use z-test approximation
        z_score = (real_rho - shuffled_mean) / shuffled_std
        # One-tailed p-value: P(shuffled >= real)
        try:
            from scipy.stats import norm
            p_value = float(1 - norm.cdf(z_score))
        except ImportError:
            # Fallback if norm not available
            p_value = 0.0 if real_rho > shuffled_mean else 1.0
    else:
        # All shuffled correlations are identical
        p_value = 0.0 if real_rho > shuffled_mean else 1.0
    
    return {
        "negative_control_rho": float(shuffled_mean),
        "negative_control_p_value": float(p_value),
        "shuffled_rho_mean": shuffled_mean,
        "shuffled_rho_std": shuffled_std,
    }
