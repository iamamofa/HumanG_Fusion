"""
Week 2: Data Integrity & Statistical Validation — Correlation Stability Analysis.

Simulates correlation stability vs overlap size to assess inference robustness.
Uses ONLY existing pipeline correlation functions.
"""

import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def simulate_correlation_stability(
    fusion_df: pd.DataFrame,
    cosmic_df: pd.DataFrame,
    overlap_sizes: List[int] = [5, 10, 20, 50, 100],
    bootstrap_n: int = 200,
) -> Optional[pd.DataFrame]:
    """
    Estimate how correlation CI width changes with overlap size.
    
    This function uses ONLY existing pipeline correlation computation methods.
    It does NOT modify core statistical logic.
    
    Args:
        fusion_df: DataFrame with gene_1, gene_2, recurrence_count columns.
        cosmic_df: DataFrame with gene_1, gene_2, recurrence_count columns.
        overlap_sizes: List of overlap sizes to simulate (default: [5, 10, 20, 50, 100]).
        bootstrap_n: Number of bootstrap iterations per overlap size (default: 200).
    
    Returns:
        DataFrame with columns: overlap_size, mean_rho, rho_std, ci_width, bootstrap_samples
        Returns None if computation fails or insufficient data.
    """
    if not SCIPY_AVAILABLE:
        print("Warning: SciPy not available. Correlation stability simulation skipped.", file=sys.stderr)
        return None
    
    try:
        # Build fusion pairs dictionary (using existing pipeline pattern)
        fusion_pairs = {}
        if "gene_1" in fusion_df.columns and "gene_2" in fusion_df.columns:
            for _, row in fusion_df.iterrows():
                gene_1 = str(row["gene_1"]).strip().upper()
                gene_2 = str(row["gene_2"]).strip().upper()
                recurrence = float(row.get("recurrence_count", 0))
                fusion_pairs[(gene_1, gene_2)] = recurrence
        
        # Build cosmic pairs dictionary
        cosmic_pairs = {}
        if "gene_1" in cosmic_df.columns and "gene_2" in cosmic_df.columns:
            for _, row in cosmic_df.iterrows():
                gene_1 = str(row["gene_1"]).strip().upper()
                gene_2 = str(row["gene_2"]).strip().upper()
                recurrence = float(row.get("recurrence_count", 0))
                cosmic_pairs[(gene_1, gene_2)] = recurrence
        
        # Find actual overlap
        fusion_set = set(fusion_pairs.keys())
        cosmic_set = set(cosmic_pairs.keys())
        actual_overlap = fusion_set & cosmic_set
        
        if len(actual_overlap) < 3:
            print("Warning: Insufficient overlap for correlation stability simulation.", file=sys.stderr)
            return None
        
        results = []
        
        # For each overlap size, simulate by sampling from actual overlap
        for target_size in overlap_sizes:
            if target_size > len(actual_overlap):
                # Cannot simulate larger than actual overlap
                continue
            
            if target_size < 3:
                # Need at least 3 for correlation
                continue
            
            # Bootstrap simulation: sample target_size pairs multiple times
            bootstrap_rhos = []
            
            rng = np.random.default_rng(42)  # Fixed seed for reproducibility
            overlap_list = list(actual_overlap)
            
            for _ in range(bootstrap_n):
                # Sample target_size pairs with replacement
                sampled_pairs = rng.choice(len(overlap_list), size=target_size, replace=True)
                sampled_overlap = {overlap_list[i] for i in sampled_pairs}
                
                # Extract counts for sampled pairs
                fusion_counts = []
                cosmic_counts = []
                
                for pair in sampled_overlap:
                    fusion_counts.append(fusion_pairs[pair])
                    cosmic_counts.append(cosmic_pairs[pair])
                
                fusion_counts = np.array(fusion_counts)
                cosmic_counts = np.array(cosmic_counts)
                
                # Compute correlation using existing method
                try:
                    rho, _ = spearmanr(fusion_counts, cosmic_counts)
                    if not np.isnan(rho) and np.isfinite(rho):
                        bootstrap_rhos.append(float(rho))
                except Exception:
                    continue
            
            if len(bootstrap_rhos) < 10:
                # Not enough valid samples
                continue
            
            bootstrap_rhos = np.array(bootstrap_rhos)
            
            # Compute statistics
            mean_rho = float(np.mean(bootstrap_rhos))
            rho_std = float(np.std(bootstrap_rhos))
            
            # Compute CI width (95% CI)
            ci_lower = float(np.percentile(bootstrap_rhos, 2.5))
            ci_upper = float(np.percentile(bootstrap_rhos, 97.5))
            ci_width = ci_upper - ci_lower
            
            results.append({
                "overlap_size": target_size,
                "mean_rho": mean_rho,
                "rho_std": rho_std,
                "ci_width": ci_width,
                "bootstrap_samples": len(bootstrap_rhos),
            })
        
        if not results:
            return None
        
        return pd.DataFrame(results)
        
    except Exception as e:
        print(f"Warning: Correlation stability simulation failed: {e}", file=sys.stderr)
        return None
