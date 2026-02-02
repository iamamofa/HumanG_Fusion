"""
COSMIC Recurrence Diagnostic - Descriptive Rank-Order Comparison with Statistical Metrics.

This module provides diagnostic comparison between fusion recurrence data and COSMIC recurrence data,
including statistical metrics for biological consistency assessment.

WHAT DOES THIS MODULE DO?
This module compares the rank order of fusion pairs between two datasets:
1. Your fusion recurrence data
2. COSMIC fusion recurrence data

STATISTICAL METRICS:
- Spearman rank correlation (rho, p-value)
- Bootstrap confidence intervals for Spearman rho
- Top fusion enrichment metric
- Distribution comparison metrics (mean, variance, zero inflation)

IT DOES NOT:
- Conclude agreement or disagreement
- Validate or invalidate data
- Make biological claims

IT PROVIDES:
- Descriptive rank ordering comparison
- Statistical metrics for consistency assessment
- Bootstrap confidence intervals for uncertainty quantification
- Counts and mismatches
"""

import hashlib
import logging
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Tuple, Any

import numpy as np
import pandas as pd

_logger = logging.getLogger(__name__)

# Import alias normalization
from week2_validation.cosmic.gene_alias_map import normalize_gene_alias

# Import negative control test
from week2_validation.cosmic.statistical_controls import compute_negative_control_correlation

# Import quality gate
from week2_validation.cosmic.quality_gate import compute_cosmic_validation_score

# Try to import scipy for Spearman correlation and hypergeometric test (optional dependency)
try:
    from scipy.stats import spearmanr, hypergeom
    import scipy
    SCIPY_AVAILABLE = True
    SCIPY_VERSION = scipy.__version__
except ImportError:
    SCIPY_AVAILABLE = False
    SCIPY_VERSION = None
    _logger.warning("scipy not available; Spearman correlation and enrichment tests will not be computed")

# Version info for reproducibility lock
NUMPY_VERSION = np.__version__
PANDAS_VERSION = pd.__version__
PYTHON_VERSION = platform.python_version()
COSMIC_VALIDATION_CODE_VERSION = "1.1.0"  # Increment when diagnostics logic changes

# Threshold for GENE_NORMALIZATION_HIGH_IMPACT warning (fraction of gene values changed)
_GENE_NORM_WARN_THRESHOLD = 0.05


def _normalize_gene_columns(df, gene_cols=("gene_1", "gene_2")):
    """
    Apply alias normalization, then .astype(str).str.upper().str.strip() to gene columns.
    Returns (df_with_normalized_columns, n_changed, n_total).
    """
    n_changed = 0
    n_total = 0
    for col in gene_cols:
        if col not in df.columns:
            continue
        before = df[col].astype(str)
        # Apply alias normalization first, then uppercase and strip
        after = before.apply(lambda x: normalize_gene_alias(str(x).strip()))
        df[col] = after
        n_total += len(df)
        n_changed += (before.str.upper().str.strip() != after).sum()
    return df, n_changed, n_total


def compute_cosmic_provenance(
    cosmic_df,
    cosmic_reference_source: Optional[str] = None,
    cosmic_file_path: Optional[Path] = None,
    mock_generation_seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Compute COSMIC provenance metadata for reproducibility tracking.
    
    Args:
        cosmic_df: COSMIC DataFrame (may be None for mock)
        cosmic_reference_source: Source identifier ("user_provided", "mock", etc.)
        cosmic_file_path: Path to COSMIC file (if available)
        mock_generation_seed: Random seed used for mock COSMIC generation
    
    Returns:
        Dictionary with provenance metadata:
        - cosmic_reference_source: Source identifier
        - cosmic_reference_version: Version string
        - cosmic_reference_file_hash: SHA256 hash of file (if available)
        - cosmic_reference_load_timestamp: ISO timestamp of load
        - reproducibility_lock: Version information for full reproducibility
    """
    provenance = {
        "cosmic_reference_source": cosmic_reference_source or "unknown",
        "cosmic_reference_version": "unknown",
        "cosmic_reference_file_hash": None,
        "cosmic_reference_load_timestamp": datetime.now().isoformat(),
    }
    
    # Set version based on source
    if cosmic_reference_source == "mock" or cosmic_reference_source == "mock_fallback":
        provenance["cosmic_reference_version"] = "mock_v1"
    elif cosmic_reference_source == "user_provided":
        provenance["cosmic_reference_version"] = "user_provided_v1"
    
    # Compute file hash if path is available
    if cosmic_file_path and cosmic_file_path.exists():
        try:
            with open(cosmic_file_path, "rb") as f:
                file_hash = hashlib.sha256(f.read()).hexdigest()
                provenance["cosmic_reference_file_hash"] = file_hash
        except Exception as e:
            _logger.warning(f"Could not compute COSMIC file hash: {e}")
    
    # Add reproducibility lock with version information
    provenance["reproducibility_lock"] = {
        "python_version": PYTHON_VERSION,
        "numpy_version": NUMPY_VERSION,
        "scipy_version": SCIPY_VERSION,
        "pandas_version": PANDAS_VERSION,
        "cosmic_validation_code_version": COSMIC_VALIDATION_CODE_VERSION,
        "random_seed_mock_generation": mock_generation_seed,
    }
    
    return provenance


def compute_bootstrap_spearman_ci(
    fusion_counts: np.ndarray,
    cosmic_counts: np.ndarray,
    n_bootstrap: int = 1000,
    confidence_level: float = 0.95,
    random_seed: int = 42,
) -> Dict[str, Optional[float]]:
    """
    Compute bootstrap confidence interval for Spearman correlation coefficient.
    
    Uses resampling with replacement to estimate the sampling distribution
    of Spearman's rho and derive confidence interval bounds.
    
    Args:
        fusion_counts: Array of fusion recurrence counts for overlapping pairs
        cosmic_counts: Array of COSMIC recurrence counts for overlapping pairs
        n_bootstrap: Number of bootstrap iterations (default: 1000)
        confidence_level: Confidence level for CI (default: 0.95 for 95% CI)
        random_seed: Random seed for reproducibility
    
    Returns:
        Dictionary with:
        - rho_ci_lower: Lower bound of CI
        - rho_ci_upper: Upper bound of CI
        - bootstrap_iterations: Number of iterations used
        - bootstrap_seed: Random seed used
    """
    if not SCIPY_AVAILABLE:
        return {
            "rho_ci_lower": None,
            "rho_ci_upper": None,
            "bootstrap_iterations": n_bootstrap,
            "bootstrap_seed": random_seed,
        }
    
    n = len(fusion_counts)
    if n < 3:
        return {
            "rho_ci_lower": None,
            "rho_ci_upper": None,
            "bootstrap_iterations": n_bootstrap,
            "bootstrap_seed": random_seed,
        }
    
    rng = np.random.default_rng(random_seed)
    bootstrap_rhos = []
    
    for _ in range(n_bootstrap):
        # Resample with replacement
        indices = rng.choice(n, size=n, replace=True)
        fusion_sample = fusion_counts[indices]
        cosmic_sample = cosmic_counts[indices]
        
        try:
            rho, _ = spearmanr(fusion_sample, cosmic_sample)
            if not np.isnan(rho):
                bootstrap_rhos.append(rho)
        except Exception:
            continue
    
    if len(bootstrap_rhos) < 10:
        # Not enough valid bootstrap samples
        return {
            "rho_ci_lower": None,
            "rho_ci_upper": None,
            "bootstrap_iterations": n_bootstrap,
            "bootstrap_seed": random_seed,
        }
    
    bootstrap_rhos = np.array(bootstrap_rhos)
    
    # Compute percentile-based CI
    alpha = 1 - confidence_level
    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100
    
    ci_lower = float(np.percentile(bootstrap_rhos, lower_percentile))
    ci_upper = float(np.percentile(bootstrap_rhos, upper_percentile))
    
    return {
        "rho_ci_lower": ci_lower,
        "rho_ci_upper": ci_upper,
        "bootstrap_iterations": n_bootstrap,
        "bootstrap_seed": random_seed,
    }


def run_cosmic_recurrence_diagnostic(
    *,
    fusion_df,
    cosmic_df,
    top_n: int = 10,
    cosmic_reference_source: Optional[str] = None,
    cosmic_file_path: Optional[Path] = None,
    mock_generation_seed: Optional[int] = 42,
    bootstrap_iterations: int = 1000,
    bootstrap_seed: int = 42,
) -> dict:
    """
    Run a descriptive rank-order comparison between fusion and COSMIC data.
    
    This function performs statistical comparison with:
    - Spearman rank correlation with bootstrap CI
    - Enrichment metrics
    - Distribution comparison
    - Quality gate scoring with component breakdown
    
    Args:
        fusion_df: DataFrame with columns: gene_1, gene_2, recurrence_count.
        cosmic_df: DataFrame with columns: gene_1, gene_2, recurrence_count.
        top_n: Number of top rank discrepancies to report (default: 10).
        cosmic_reference_source: Source identifier ("user_provided", "mock", etc.)
        cosmic_file_path: Path to COSMIC file (if available)
        mock_generation_seed: Random seed used for mock COSMIC generation
        bootstrap_iterations: Number of bootstrap iterations for CI (default: 1000)
        bootstrap_seed: Random seed for bootstrap reproducibility (default: 42)
    
    Returns:
        A dictionary containing:
        - total_fusions_ours: Count of fusion pairs in our data
        - total_fusions_cosmic: Count of fusion pairs in COSMIC data
        - overlap_count: Count of pairs present in both datasets
        - only_in_ours_count: Count of pairs only in our data
        - only_in_cosmic_count: Count of pairs only in COSMIC data
        - top_rank_discrepancies: List of dicts with top rank mismatches
        - Statistical metrics with bootstrap CI
        - Provenance metadata with reproducibility lock
    """
    # -------------------------------------------------------------------------
    # STEP 1: Validate required columns
    # -------------------------------------------------------------------------
    required_columns = {"gene_1", "gene_2", "recurrence_count"}
    
    fusion_columns = set(fusion_df.columns) if fusion_df is not None else set()
    cosmic_columns = set(cosmic_df.columns) if cosmic_df is not None else set()
    
    # Check fusion_df columns
    if fusion_df is None:
        return {
            "total_fusions_ours": 0,
            "total_fusions_cosmic": 0,
            "overlap_count": 0,
            "only_in_ours_count": 0,
            "only_in_cosmic_count": 0,
            "top_rank_discrepancies": [],
            "message": "Fusion DataFrame is None.",
        }
    
    missing_fusion_cols = required_columns - fusion_columns
    if missing_fusion_cols:
        return {
            "total_fusions_ours": 0,
            "total_fusions_cosmic": 0,
            "overlap_count": 0,
            "only_in_ours_count": 0,
            "only_in_cosmic_count": 0,
            "top_rank_discrepancies": [],
            "message": f"Fusion data missing required columns: {missing_fusion_cols}",
        }
    
    # Check cosmic_df columns
    if cosmic_df is None:
        return {
            "total_fusions_ours": len(fusion_df),
            "total_fusions_cosmic": 0,
            "overlap_count": 0,
            "only_in_ours_count": len(fusion_df),
            "only_in_cosmic_count": 0,
            "top_rank_discrepancies": [],
            "message": "COSMIC DataFrame is None.",
        }
    
    missing_cosmic_cols = required_columns - cosmic_columns
    if missing_cosmic_cols:
        return {
            "total_fusions_ours": len(fusion_df),
            "total_fusions_cosmic": 0,
            "overlap_count": 0,
            "only_in_ours_count": len(fusion_df),
            "only_in_cosmic_count": 0,
            "top_rank_discrepancies": [],
            "message": f"COSMIC data missing required columns: {missing_cosmic_cols}",
        }

    # -------------------------------------------------------------------------
    # STEP 1.5: Normalize gene_1, gene_2 (uppercase, strip) before pair creation
    # -------------------------------------------------------------------------
    fusion_df, n_fusion_changed, n_fusion_total = _normalize_gene_columns(fusion_df.copy())
    cosmic_df, n_cosmic_changed, n_cosmic_total = _normalize_gene_columns(cosmic_df.copy())
    n_changed_total = n_fusion_changed + n_cosmic_changed
    n_gene_values_total = n_fusion_total + n_cosmic_total
    if n_gene_values_total > 0 and n_changed_total / n_gene_values_total > _GENE_NORM_WARN_THRESHOLD:
        _logger.warning("GENE_NORMALIZATION_HIGH_IMPACT")

    # -------------------------------------------------------------------------
    # STEP 2: Normalize fusion pairs (uppercase, strip, sort alphabetically)
    # -------------------------------------------------------------------------
    def normalize_pair(gene_1, gene_2):
        """Normalize a gene pair: uppercase, strip, sort alphabetically."""
        g1 = str(gene_1).strip().upper()
        g2 = str(gene_2).strip().upper()
        # Return sorted tuple for consistent ordering
        return tuple(sorted([g1, g2]))
    
    # Build normalized lookup for fusion data
    fusion_pairs = {}
    for idx, row in fusion_df.iterrows():
        pair = normalize_pair(row["gene_1"], row["gene_2"])
        count = row["recurrence_count"]
        # If duplicate pairs, keep the one with higher count
        if pair not in fusion_pairs or count > fusion_pairs[pair]:
            fusion_pairs[pair] = count
    
    # Build normalized lookup for COSMIC data
    cosmic_pairs = {}
    for idx, row in cosmic_df.iterrows():
        pair = normalize_pair(row["gene_1"], row["gene_2"])
        count = row["recurrence_count"]
        # If duplicate pairs, keep the one with higher count
        if pair not in cosmic_pairs or count > cosmic_pairs[pair]:
            cosmic_pairs[pair] = count
    
    # -------------------------------------------------------------------------
    # STEP 3: Compute counts
    # -------------------------------------------------------------------------
    fusion_set = set(fusion_pairs.keys())
    cosmic_set = set(cosmic_pairs.keys())
    
    overlap = fusion_set & cosmic_set
    only_in_ours = fusion_set - cosmic_set
    only_in_cosmic = cosmic_set - fusion_set
    
    total_fusions_ours = len(fusion_set)
    total_fusions_cosmic = len(cosmic_set)
    overlap_count = len(overlap)
    only_in_ours_count = len(only_in_ours)
    only_in_cosmic_count = len(only_in_cosmic)
    
    # -------------------------------------------------------------------------
    # STEP 4: Compute ranks for overlapping pairs
    # -------------------------------------------------------------------------
    if overlap_count == 0:
        # No overlap - nothing to compare
        return {
            "total_fusions_ours": total_fusions_ours,
            "total_fusions_cosmic": total_fusions_cosmic,
            "overlap_count": overlap_count,
            "only_in_ours_count": only_in_ours_count,
            "only_in_cosmic_count": only_in_cosmic_count,
            "top_rank_discrepancies": [],
            "message": "No overlapping fusion pairs found between datasets.",
        }
    
    # Rank fusion pairs by recurrence_count (descending)
    # Rank 1 = highest recurrence count
    fusion_sorted = sorted(fusion_pairs.items(), key=lambda x: x[1], reverse=True)
    fusion_ranks = {pair: rank + 1 for rank, (pair, _) in enumerate(fusion_sorted)}
    
    cosmic_sorted = sorted(cosmic_pairs.items(), key=lambda x: x[1], reverse=True)
    cosmic_ranks = {pair: rank + 1 for rank, (pair, _) in enumerate(cosmic_sorted)}
    
    # -------------------------------------------------------------------------
    # STEP 5: Compute rank discrepancies for overlapping pairs
    # -------------------------------------------------------------------------
    discrepancies = []
    for pair in overlap:
        our_rank = fusion_ranks[pair]
        cosmic_rank = cosmic_ranks[pair]
        abs_diff = abs(our_rank - cosmic_rank)
        discrepancies.append({
            "gene_1": pair[0],
            "gene_2": pair[1],
            "our_rank": our_rank,
            "cosmic_rank": cosmic_rank,
            "absolute_rank_difference": abs_diff,
            "our_recurrence_count": fusion_pairs[pair],
            "cosmic_recurrence_count": cosmic_pairs[pair],
        })
    
    # Sort by absolute rank difference (descending) to get top mismatches
    discrepancies.sort(key=lambda x: x["absolute_rank_difference"], reverse=True)
    
    # Take top N
    top_discrepancies = discrepancies[:top_n]
    
    # -------------------------------------------------------------------------
    # STEP 6: Compute statistical metrics
    # -------------------------------------------------------------------------
    stats = compute_cosmic_statistical_metrics(
        fusion_pairs=fusion_pairs,
        cosmic_pairs=cosmic_pairs,
        overlap=overlap,
        bootstrap_iterations=bootstrap_iterations,
        bootstrap_seed=bootstrap_seed,
    )
    
    # -------------------------------------------------------------------------
    # STEP 7: Compute provenance metadata
    # -------------------------------------------------------------------------
    provenance = compute_cosmic_provenance(
        cosmic_df=cosmic_df,
        cosmic_reference_source=cosmic_reference_source,
        cosmic_file_path=cosmic_file_path,
        mock_generation_seed=mock_generation_seed,
    )
    
    # -------------------------------------------------------------------------
    # STEP 8: Return descriptive dictionary with statistical metrics and provenance
    # -------------------------------------------------------------------------
    result = {
        "total_fusions_ours": total_fusions_ours,
        "total_fusions_cosmic": total_fusions_cosmic,
        "overlap_count": overlap_count,
        "only_in_ours_count": only_in_ours_count,
        "only_in_cosmic_count": only_in_cosmic_count,
        "top_rank_discrepancies": top_discrepancies,
    }
    
    # Add statistical metrics
    result.update(stats)
    
    # Add provenance metadata
    result.update(provenance)
    
    return result


def compute_spearman_correlation(
    fusion_pairs: Dict[Tuple[str, str], float],
    cosmic_pairs: Dict[Tuple[str, str], float],
    overlap: set,
    compute_bootstrap_ci: bool = True,
    bootstrap_iterations: int = 1000,
    bootstrap_seed: int = 42,
) -> Dict[str, Optional[float]]:
    """
    Compute Spearman rank correlation between dataset and COSMIC recurrence counts.
    
    Args:
        fusion_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        cosmic_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        overlap: Set of overlapping fusion pairs (gene_1, gene_2) tuples
        compute_bootstrap_ci: Whether to compute bootstrap CI (default: True)
        bootstrap_iterations: Number of bootstrap iterations (default: 1000)
        bootstrap_seed: Random seed for bootstrap (default: 42)
    
    Returns:
        Dictionary with:
        - spearman_rho: Spearman correlation coefficient (or None if cannot compute)
        - spearman_p_value: p-value for Spearman test (or None if cannot compute)
        - rho_ci_lower: Lower bound of 95% bootstrap CI (or None)
        - rho_ci_upper: Upper bound of 95% bootstrap CI (or None)
        - bootstrap_iterations: Number of bootstrap iterations used
    """
    base_result = {
        "spearman_rho": None,
        "spearman_p_value": None,
        "rho_ci_lower": None,
        "rho_ci_upper": None,
        "bootstrap_iterations": bootstrap_iterations,
    }
    
    if not SCIPY_AVAILABLE:
        return base_result
    
    if len(overlap) < 3:
        # Need at least 3 pairs for meaningful correlation
        return base_result
    
    # Extract recurrence counts for overlapping pairs
    fusion_counts = []
    cosmic_counts = []
    
    for pair in overlap:
        fusion_counts.append(fusion_pairs[pair])
        cosmic_counts.append(cosmic_pairs[pair])
    
    fusion_counts = np.array(fusion_counts)
    cosmic_counts = np.array(cosmic_counts)
    
    try:
        rho, p_value = spearmanr(fusion_counts, cosmic_counts)
        
        # Handle edge cases
        if np.isnan(rho) or np.isnan(p_value):
            return base_result
        
        result = {
            "spearman_rho": float(rho),
            "spearman_p_value": float(p_value),
            "rho_ci_lower": None,
            "rho_ci_upper": None,
            "bootstrap_iterations": bootstrap_iterations,
        }
        
        # Compute bootstrap CI if requested
        if compute_bootstrap_ci:
            ci_result = compute_bootstrap_spearman_ci(
                fusion_counts=fusion_counts,
                cosmic_counts=cosmic_counts,
                n_bootstrap=bootstrap_iterations,
                confidence_level=0.95,
                random_seed=bootstrap_seed,
            )
            result["rho_ci_lower"] = ci_result["rho_ci_lower"]
            result["rho_ci_upper"] = ci_result["rho_ci_upper"]
        
        return result
    except Exception as e:
        _logger.warning(f"Error computing Spearman correlation: {e}")
        return base_result


def compute_enrichment_metric(
    fusion_pairs: Dict[Tuple[str, str], float],
    cosmic_pairs: Dict[Tuple[str, str], float],
    top_n: int = 10,
) -> Dict[str, int]:
    """
    Compute enrichment metric: overlap between top N fusions in dataset vs COSMIC.
    
    Args:
        fusion_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        cosmic_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        top_n: Number of top fusions to compare (default: 10)
    
    Returns:
        Dictionary with:
        - top_fusion_overlap: Number of overlapping fusions in top N
        - top_fusion_enrichment_ratio: Ratio of overlap to top_n (0.0 to 1.0)
    """
    # Get top N fusions from each dataset
    fusion_sorted = sorted(fusion_pairs.items(), key=lambda x: x[1], reverse=True)
    cosmic_sorted = sorted(cosmic_pairs.items(), key=lambda x: x[1], reverse=True)
    
    top_fusion = set(pair for pair, _ in fusion_sorted[:top_n])
    top_cosmic = set(pair for pair, _ in cosmic_sorted[:top_n])
    
    overlap_count = len(top_fusion & top_cosmic)
    enrichment_ratio = overlap_count / top_n if top_n > 0 else 0.0
    
    return {
        "top_fusion_overlap": overlap_count,
        "top_fusion_enrichment_ratio": float(enrichment_ratio),
    }


def compute_hypergeometric_enrichment(
    fusion_pairs: Dict[Tuple[str, str], float],
    cosmic_pairs: Dict[Tuple[str, str], float],
    top_n: int = 10,
) -> Dict[str, Optional[float]]:
    """
    Compute hypergeometric enrichment test for top fusion overlap.
    
    Tests whether the overlap between top N fusions in dataset and COSMIC
    is significantly greater than expected by chance.
    
    Args:
        fusion_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        cosmic_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        top_n: Number of top fusions to compare (default: 10)
    
    Returns:
        Dictionary with:
        - enrichment_p_value: p-value from hypergeometric test (or None)
        - observed_overlap: Number of overlapping fusions in top N
        - expected_overlap_random: Expected overlap by chance
    """
    if not SCIPY_AVAILABLE:
        return {
            "enrichment_p_value": None,
            "observed_overlap": None,
            "expected_overlap_random": None,
        }
    
    # Get top N fusions from each dataset
    fusion_sorted = sorted(fusion_pairs.items(), key=lambda x: x[1], reverse=True)
    cosmic_sorted = sorted(cosmic_pairs.items(), key=lambda x: x[1], reverse=True)
    
    top_fusion = set(pair for pair, _ in fusion_sorted[:top_n])
    top_cosmic = set(pair for pair, _ in cosmic_sorted[:top_n])
    
    # Observed overlap
    observed_overlap = len(top_fusion & top_cosmic)
    
    # Total population: all unique fusion pairs across both datasets
    all_pairs = set(fusion_pairs.keys()) | set(cosmic_pairs.keys())
    N = len(all_pairs)  # Population size
    
    # Hypergeometric parameters:
    # N: total population (all unique fusion pairs)
    # K: number of "success" states in population (cosmic top N)
    # n: number of draws (fusion top N)
    # k: observed successes (overlap)
    K = len(top_cosmic)
    n = len(top_fusion)
    k = observed_overlap
    
    # Expected overlap by chance: E[k] = n * K / N
    expected_overlap = (n * K / N) if N > 0 else 0.0
    
    # Compute p-value: P(X >= k) where X ~ Hypergeometric(N, K, n)
    # This is a one-tailed test for enrichment
    try:
        # sf(k-1) gives P(X >= k)
        p_value = float(hypergeom.sf(k - 1, N, K, n))
        
        return {
            "enrichment_p_value": p_value,
            "observed_overlap": k,
            "expected_overlap_random": float(expected_overlap),
        }
    except Exception as e:
        _logger.warning(f"Error computing hypergeometric enrichment: {e}")
        return {
            "enrichment_p_value": None,
            "observed_overlap": k,
            "expected_overlap_random": float(expected_overlap),
        }


def compute_distribution_metrics(
    fusion_pairs: Dict[Tuple[str, str], float],
    cosmic_pairs: Dict[Tuple[str, str], float],
) -> Dict[str, float]:
    """
    Compute distribution comparison metrics between dataset and COSMIC.
    
    Args:
        fusion_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        cosmic_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
    
    Returns:
        Dictionary with:
        - mean_recurrence_fusion: Mean recurrence count in fusion dataset
        - mean_recurrence_cosmic: Mean recurrence count in COSMIC dataset
        - variance_fusion: Variance of recurrence counts in fusion dataset
        - variance_cosmic: Variance of recurrence counts in COSMIC dataset
        - zero_inflation_rate_fusion: Fraction of fusions with recurrence_count == 0
        - zero_inflation_rate_cosmic: Fraction of fusions with recurrence_count == 0
    """
    fusion_counts = np.array(list(fusion_pairs.values()))
    cosmic_counts = np.array(list(cosmic_pairs.values()))
    
    # Mean
    mean_fusion = float(np.mean(fusion_counts)) if len(fusion_counts) > 0 else 0.0
    mean_cosmic = float(np.mean(cosmic_counts)) if len(cosmic_counts) > 0 else 0.0
    
    # Variance
    var_fusion = float(np.var(fusion_counts)) if len(fusion_counts) > 0 else 0.0
    var_cosmic = float(np.var(cosmic_counts)) if len(cosmic_counts) > 0 else 0.0
    
    # Zero inflation rate
    zero_rate_fusion = float(np.sum(fusion_counts == 0) / len(fusion_counts)) if len(fusion_counts) > 0 else 0.0
    zero_rate_cosmic = float(np.sum(cosmic_counts == 0) / len(cosmic_counts)) if len(cosmic_counts) > 0 else 0.0
    
    return {
        "mean_recurrence_fusion": mean_fusion,
        "mean_recurrence_cosmic": mean_cosmic,
        "variance_fusion": var_fusion,
        "variance_cosmic": var_cosmic,
        "zero_inflation_rate_fusion": zero_rate_fusion,
        "zero_inflation_rate_cosmic": zero_rate_cosmic,
    }


def compute_cosmic_statistical_metrics(
    fusion_pairs: Dict[Tuple[str, str], float],
    cosmic_pairs: Dict[Tuple[str, str], float],
    overlap: set,
    bootstrap_iterations: int = 1000,
    bootstrap_seed: int = 42,
) -> Dict[str, Any]:
    """
    Compute all statistical metrics for COSMIC cross-validation.
    
    Args:
        fusion_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        cosmic_pairs: Dictionary mapping (gene_1, gene_2) tuples to recurrence_count
        overlap: Set of overlapping fusion pairs
        bootstrap_iterations: Number of bootstrap iterations for CI (default: 1000)
        bootstrap_seed: Random seed for bootstrap reproducibility (default: 42)
    
    Returns:
        Dictionary containing all statistical metrics including:
        - Spearman correlation with bootstrap CI
        - Enrichment metrics
        - Distribution metrics
        - Quality gate score with component breakdown
    """
    metrics = {}
    
    # Spearman correlation with bootstrap CI
    spearman = compute_spearman_correlation(
        fusion_pairs, 
        cosmic_pairs, 
        overlap,
        compute_bootstrap_ci=True,
        bootstrap_iterations=bootstrap_iterations,
        bootstrap_seed=bootstrap_seed,
    )
    metrics.update(spearman)
    
    # Enrichment metric
    enrichment = compute_enrichment_metric(fusion_pairs, cosmic_pairs, top_n=10)
    metrics.update(enrichment)
    
    # Hypergeometric enrichment test (statistical significance)
    hypergeom_enrichment = compute_hypergeometric_enrichment(fusion_pairs, cosmic_pairs, top_n=10)
    metrics.update(hypergeom_enrichment)
    
    # Negative control test (shuffled correlation)
    negative_control = compute_negative_control_correlation(fusion_pairs, cosmic_pairs, overlap)
    metrics.update(negative_control)
    
    # Distribution metrics
    distribution = compute_distribution_metrics(fusion_pairs, cosmic_pairs)
    metrics.update(distribution)
    
    # Quality gate: compute validation score and classification with component breakdown
    quality_gate = compute_cosmic_validation_score(
        spearman_rho=metrics.get("spearman_rho"),
        spearman_p_value=metrics.get("spearman_p_value"),
        enrichment_p_value=metrics.get("enrichment_p_value"),
        negative_control_rho=metrics.get("negative_control_rho"),
        observed_overlap=metrics.get("observed_overlap"),
        expected_overlap_random=metrics.get("expected_overlap_random"),
    )
    metrics.update(quality_gate)
    
    return metrics
