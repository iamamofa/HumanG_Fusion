"""
Week 2: Data Integrity & Statistical Validation — External Validity Stability Test.

This module tests distribution generalization by checking if recurrence signal
is stable across different recurrence strata (high, middle, low).

SAFETY: This module does NOT modify real data or core correlation formulas.
It uses ONLY existing pipeline correlation functions.
"""

import sys
from typing import Dict, Optional

import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# Import existing correlation function to ensure consistency
from week2_validation.cosmic.diagnostics import compute_spearman_correlation


def assess_distribution_generalization(
    fusion_df: pd.DataFrame,
    cosmic_df: pd.DataFrame,
) -> Dict[str, Optional[float]]:
    """
    Test if recurrence signal is stable across recurrence strata.
    
    PURPOSE:
    Check if correlation is stable across:
    - High recurrence fusions (top 25%)
    - Middle recurrence fusions (middle 50%)
    - Low recurrence fusions (bottom 25%)
    
    METHOD:
    Split fusion dataset into strata by recurrence_count, then compute
    correlation separately for each stratum.
    
    Args:
        fusion_df: DataFrame with columns: gene_1, gene_2, recurrence_count
        cosmic_df: DataFrame with columns: gene_1, gene_2, recurrence_count
    
    Returns:
        Dictionary with:
        - high_recurrence_rho: Correlation for top 25% recurrence fusions
        - mid_recurrence_rho: Correlation for middle 50% recurrence fusions
        - low_recurrence_rho: Correlation for bottom 25% recurrence fusions
        - stability_index: Standard deviation of [high, mid, low] correlations
        - n_high: Number of fusions in high stratum
        - n_mid: Number of fusions in middle stratum
        - n_low: Number of fusions in low stratum
    """
    base_result = {
        "high_recurrence_rho": None,
        "mid_recurrence_rho": None,
        "low_recurrence_rho": None,
        "stability_index": None,
        "n_high": 0,
        "n_mid": 0,
        "n_low": 0,
    }
    
    if not SCIPY_AVAILABLE:
        print("Warning: SciPy not available. Distribution generalization test skipped.", file=sys.stderr)
        return base_result
    
    try:
        # Build fusion pairs dictionary (using existing pipeline pattern)
        fusion_pairs = {}
        if "gene_1" in fusion_df.columns and "gene_2" in fusion_df.columns:
            for _, row in fusion_df.iterrows():
                gene_1 = str(row["gene_1"]).strip().upper()
                gene_2 = str(row["gene_2"]).strip().upper()
                recurrence = float(row.get("recurrence_count", 0))
                pair = tuple(sorted([gene_1, gene_2]))
                # If duplicate pairs, keep the one with higher count
                if pair not in fusion_pairs or recurrence > fusion_pairs[pair]:
                    fusion_pairs[pair] = recurrence
        
        # Build cosmic pairs dictionary
        cosmic_pairs = {}
        if "gene_1" in cosmic_df.columns and "gene_2" in cosmic_df.columns:
            for _, row in cosmic_df.iterrows():
                gene_1 = str(row["gene_1"]).strip().upper()
                gene_2 = str(row["gene_2"]).strip().upper()
                recurrence = float(row.get("recurrence_count", 0))
                pair = tuple(sorted([gene_1, gene_2]))
                # If duplicate pairs, keep the one with higher count
                if pair not in cosmic_pairs or recurrence > cosmic_pairs[pair]:
                    cosmic_pairs[pair] = recurrence
        
        # Sort fusion pairs by recurrence_count (descending)
        fusion_sorted = sorted(fusion_pairs.items(), key=lambda x: x[1], reverse=True)
        n_total = len(fusion_sorted)
        
        if n_total < 4:
            print("Warning: Insufficient fusions for stratification.", file=sys.stderr)
            return base_result
        
        # Define strata boundaries
        # Top 25%: highest recurrence
        # Middle 50%: middle recurrence
        # Bottom 25%: lowest recurrence
        n_high = max(1, int(n_total * 0.25))
        n_low = max(1, int(n_total * 0.25))
        n_mid = n_total - n_high - n_low
        
        if n_mid < 1:
            # If too few, adjust: high = top 1/3, mid = middle 1/3, low = bottom 1/3
            n_high = max(1, n_total // 3)
            n_low = max(1, n_total // 3)
            n_mid = n_total - n_high - n_low
        
        # Split into strata
        high_pairs = set(pair for pair, _ in fusion_sorted[:n_high])
        mid_pairs = set(pair for pair, _ in fusion_sorted[n_high:n_high + n_mid])
        low_pairs = set(pair for pair, _ in fusion_sorted[n_high + n_mid:])
        
        # Find overlap for each stratum
        cosmic_set = set(cosmic_pairs.keys())
        high_overlap = high_pairs & cosmic_set
        mid_overlap = mid_pairs & cosmic_set
        low_overlap = low_pairs & cosmic_set
        
        # Compute correlation for each stratum using existing function
        rhos = []
        
        # High recurrence stratum
        if len(high_overlap) >= 3:
            high_result = compute_spearman_correlation(
                fusion_pairs=fusion_pairs,
                cosmic_pairs=cosmic_pairs,
                overlap=high_overlap,
                compute_bootstrap_ci=False,
            )
            high_rho = high_result.get("spearman_rho")
            if high_rho is not None and not np.isnan(high_rho):
                rhos.append(high_rho)
                base_result["high_recurrence_rho"] = float(high_rho)
                base_result["n_high"] = len(high_overlap)
        
        # Middle recurrence stratum
        if len(mid_overlap) >= 3:
            mid_result = compute_spearman_correlation(
                fusion_pairs=fusion_pairs,
                cosmic_pairs=cosmic_pairs,
                overlap=mid_overlap,
                compute_bootstrap_ci=False,
            )
            mid_rho = mid_result.get("spearman_rho")
            if mid_rho is not None and not np.isnan(mid_rho):
                rhos.append(mid_rho)
                base_result["mid_recurrence_rho"] = float(mid_rho)
                base_result["n_mid"] = len(mid_overlap)
        
        # Low recurrence stratum
        if len(low_overlap) >= 3:
            low_result = compute_spearman_correlation(
                fusion_pairs=fusion_pairs,
                cosmic_pairs=cosmic_pairs,
                overlap=low_overlap,
                compute_bootstrap_ci=False,
            )
            low_rho = low_result.get("spearman_rho")
            if low_rho is not None and not np.isnan(low_rho):
                rhos.append(low_rho)
                base_result["low_recurrence_rho"] = float(low_rho)
                base_result["n_low"] = len(low_overlap)
        
        # Compute stability index: std([high, mid, low])
        if len(rhos) >= 2:
            stability_index = float(np.std(rhos))
            base_result["stability_index"] = stability_index
        elif len(rhos) == 1:
            # Only one stratum has sufficient overlap
            base_result["stability_index"] = None
        
        return base_result
        
    except Exception as e:
        print(f"Warning: Distribution generalization test failed: {e}", file=sys.stderr)
        return base_result
