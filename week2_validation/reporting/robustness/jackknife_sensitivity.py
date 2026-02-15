"""
Week 2: Data Integrity & Statistical Validation — Jackknife Sensitivity Analysis.

Performs jackknife sensitivity analysis by removing each overlapping fusion pair
and recomputing correlation to assess sensitivity.
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


def compute_cosmic_jackknife_sensitivity(
    fusion_df: pd.DataFrame,
    cosmic_df: pd.DataFrame,
) -> Optional[pd.DataFrame]:
    """
    Compute jackknife sensitivity: remove each overlapping fusion and recompute correlation.
    
    This function uses ONLY existing pipeline correlation computation methods.
    It does NOT modify core statistical logic.
    
    Args:
        fusion_df: DataFrame with gene_1, gene_2, recurrence_count columns.
        cosmic_df: DataFrame with gene_1, gene_2, recurrence_count columns.
    
    Returns:
        DataFrame with columns: fusion_pair_removed, rho_after_removal, delta_from_full_rho
        Returns None if computation fails or insufficient overlap.
    """
    if not SCIPY_AVAILABLE:
        print("Warning: SciPy not available. Jackknife sensitivity analysis skipped.", file=sys.stderr)
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
        
        if len(actual_overlap) < 4:
            # Need at least 4 pairs: 3 for correlation after removal, 1 to remove
            print("Warning: Insufficient overlap for jackknife analysis (need at least 4 pairs).", file=sys.stderr)
            return None
        
        # Compute full correlation (baseline)
        fusion_counts_full = []
        cosmic_counts_full = []
        
        for pair in actual_overlap:
            fusion_counts_full.append(fusion_pairs[pair])
            cosmic_counts_full.append(cosmic_pairs[pair])
        
        fusion_counts_full = np.array(fusion_counts_full)
        cosmic_counts_full = np.array(cosmic_counts_full)
        
        try:
            rho_full, _ = spearmanr(fusion_counts_full, cosmic_counts_full)
            if np.isnan(rho_full) or not np.isfinite(rho_full):
                return None
            rho_full = float(rho_full)
        except Exception:
            return None
        
        # Jackknife: remove each pair and recompute
        results = []
        
        for pair_to_remove in actual_overlap:
            # Create overlap without this pair
            overlap_minus_one = actual_overlap - {pair_to_remove}
            
            if len(overlap_minus_one) < 3:
                # Need at least 3 pairs for correlation
                continue
            
            # Extract counts
            fusion_counts = []
            cosmic_counts = []
            
            for pair in overlap_minus_one:
                fusion_counts.append(fusion_pairs[pair])
                cosmic_counts.append(cosmic_pairs[pair])
            
            fusion_counts = np.array(fusion_counts)
            cosmic_counts = np.array(cosmic_counts)
            
            # Compute correlation
            try:
                rho_removed, _ = spearmanr(fusion_counts, cosmic_counts)
                if np.isnan(rho_removed) or not np.isfinite(rho_removed):
                    continue
                rho_removed = float(rho_removed)
            except Exception:
                continue
            
            # Compute delta
            delta = rho_removed - rho_full
            
            # Format pair name
            pair_name = f"{pair_to_remove[0]}-{pair_to_remove[1]}"
            
            results.append({
                "fusion_pair_removed": pair_name,
                "rho_after_removal": rho_removed,
                "delta_from_full_rho": delta,
            })
        
        if not results:
            return None
        
        return pd.DataFrame(results)
        
    except Exception as e:
        print(f"Warning: Jackknife sensitivity analysis failed: {e}", file=sys.stderr)
        return None
