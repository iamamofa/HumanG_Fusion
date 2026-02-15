"""
Week 2: Data Integrity & Statistical Validation — Null Model Simulation (Statistical Falsification Test).

This module implements a null model simulation to test whether observed correlation
could happen by chance. It shuffles COSMIC recurrence counts and recomputes Spearman
correlation to build a null distribution.

SAFETY: This module does NOT modify real data or core correlation formulas.
It uses ONLY existing pipeline correlation functions.
"""

import sys
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# Import existing correlation function to ensure consistency
from week2_validation.cosmic.diagnostics import compute_spearman_correlation


def run_cosmic_null_model_test(
    fusion_df: pd.DataFrame,
    cosmic_df: pd.DataFrame,
    n_simulations: int = 1000,
) -> Dict[str, Optional[float]]:
    """
    Run null model simulation test for COSMIC correlation.
    
    Tests the hypothesis: "Could observed correlation happen by chance?"
    
    METHOD (SAFE):
    - Does NOT change real data
    - For each simulation:
      1. Shuffle COSMIC recurrence counts
      2. Recompute Spearman using existing function
      3. Store rho
    
    Args:
        fusion_df: DataFrame with columns: gene_1, gene_2, recurrence_count
        cosmic_df: DataFrame with columns: gene_1, gene_2, recurrence_count
        n_simulations: Number of null simulations (default: 1000)
    
    Returns:
        Dictionary with:
        - observed_rho: Observed Spearman correlation
        - null_mean_rho: Mean of null distribution
        - null_std_rho: Standard deviation of null distribution
        - empirical_p_value: Empirical p-value (fraction of null >= observed)
        - z_score_vs_null: Z-score of observed rho vs null distribution
        - n_simulations: Number of simulations run
    """
    base_result = {
        "observed_rho": None,
        "null_mean_rho": None,
        "null_std_rho": None,
        "empirical_p_value": None,
        "z_score_vs_null": None,
        "n_simulations": n_simulations,
    }
    
    if not SCIPY_AVAILABLE:
        print("Warning: SciPy not available. Null model test skipped.", file=sys.stderr)
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
        
        # Find actual overlap
        fusion_set = set(fusion_pairs.keys())
        cosmic_set = set(cosmic_pairs.keys())
        overlap = fusion_set & cosmic_set
        
        if len(overlap) < 3:
            print("Warning: Insufficient overlap for null model test.", file=sys.stderr)
            return base_result
        
        # Compute observed correlation using existing function
        observed_result = compute_spearman_correlation(
            fusion_pairs=fusion_pairs,
            cosmic_pairs=cosmic_pairs,
            overlap=overlap,
            compute_bootstrap_ci=False,  # Don't need CI for null model
        )
        
        observed_rho = observed_result.get("spearman_rho")
        if observed_rho is None or np.isnan(observed_rho):
            return base_result
        
        # Extract recurrence counts for overlapping pairs (in consistent order)
        overlap_list = sorted(list(overlap))  # Sort for consistency
        fusion_counts_original = np.array([fusion_pairs[pair] for pair in overlap_list])
        cosmic_counts_original = np.array([cosmic_pairs[pair] for pair in overlap_list])
        
        # Run null simulations: shuffle COSMIC counts
        null_rhos = []
        rng = np.random.default_rng(42)  # Fixed seed for reproducibility
        
        for _ in range(n_simulations):
            # Shuffle COSMIC recurrence counts (permutation test)
            cosmic_counts_shuffled = rng.permutation(cosmic_counts_original)
            
            # Recompute correlation with shuffled COSMIC
            try:
                rho_shuffled, _ = spearmanr(fusion_counts_original, cosmic_counts_shuffled)
                if not np.isnan(rho_shuffled) and np.isfinite(rho_shuffled):
                    null_rhos.append(float(rho_shuffled))
            except Exception:
                continue
        
        if len(null_rhos) < 10:
            print("Warning: Insufficient valid null simulations.", file=sys.stderr)
            return base_result
        
        null_rhos = np.array(null_rhos)
        
        # Compute null distribution statistics
        null_mean = float(np.mean(null_rhos))
        null_std = float(np.std(null_rhos))
        
        # Compute empirical p-value: fraction of null >= |observed|
        # Two-tailed test: how extreme is observed correlation?
        abs_observed = abs(observed_rho)
        abs_null = np.abs(null_rhos)
        empirical_p_value = float(np.mean(abs_null >= abs_observed))
        
        # Compute Z-score: (observed - null_mean) / null_std
        if null_std > 0:
            z_score = float((observed_rho - null_mean) / null_std)
        else:
            z_score = None
        
        return {
            "observed_rho": float(observed_rho),
            "null_mean_rho": null_mean,
            "null_std_rho": null_std,
            "empirical_p_value": empirical_p_value,
            "z_score_vs_null": z_score,
            "n_simulations": len(null_rhos),
        }
        
    except Exception as e:
        print(f"Warning: Null model test failed: {e}", file=sys.stderr)
        return base_result
