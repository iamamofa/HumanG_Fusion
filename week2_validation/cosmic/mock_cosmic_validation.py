"""
Mock COSMIC Realism Validation
Week 2: Data Integrity & Statistical Validation

This module computes distribution metrics for the mock COSMIC dataset to validate
that its synthetic distribution approximates realistic recurrence patterns.

WHAT DOES THIS MODULE DO?
Computes:
- Skewness
- Kurtosis
- Tail heaviness estimate
- Zero inflation rate
- Gini coefficient (recurrence inequality)

These metrics verify that mock COSMIC preserves:
- Heavy tail recurrence distribution
- Driver gene enrichment structure
- Non-uniform recurrence behavior

DESIGN PRINCIPLE:
Mock COSMIC exists for pipeline operability, but its realism should be validated
and documented for scientific transparency.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Any

import numpy as np
import pandas as pd

_logger = logging.getLogger(__name__)

# Try to import scipy for statistical functions
try:
    from scipy.stats import skew, kurtosis
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    _logger.warning("scipy not available; skewness and kurtosis will use numpy fallback")


def compute_gini_coefficient(values: np.ndarray) -> float:
    """
    Compute the Gini coefficient for measuring inequality in recurrence counts.
    
    Gini = 0 means perfect equality (all fusions have same recurrence)
    Gini = 1 means maximum inequality (one fusion has all recurrence)
    
    Real COSMIC data typically has high Gini (0.7-0.9) due to driver fusion dominance.
    
    Args:
        values: Array of recurrence counts
    
    Returns:
        Gini coefficient (0.0 to 1.0)
    """
    if len(values) == 0:
        return 0.0
    
    # Sort values
    sorted_values = np.sort(values)
    n = len(sorted_values)
    
    # Compute Gini using the relative mean absolute difference formula
    cumulative = np.cumsum(sorted_values)
    sum_values = cumulative[-1]
    
    if sum_values == 0:
        return 0.0
    
    # Gini = (2 * sum(i * x_i) - (n + 1) * sum(x_i)) / (n * sum(x_i))
    index_sum = np.sum((np.arange(1, n + 1)) * sorted_values)
    gini = (2 * index_sum - (n + 1) * sum_values) / (n * sum_values)
    
    return float(gini)


def compute_tail_heaviness(values: np.ndarray, threshold_percentile: float = 90) -> Dict[str, float]:
    """
    Compute tail heaviness metrics for recurrence distribution.
    
    Heavy tails are characteristic of biological recurrence data where
    few fusions (driver fusions) have very high recurrence.
    
    Args:
        values: Array of recurrence counts
        threshold_percentile: Percentile threshold for defining "tail" (default: 90)
    
    Returns:
        Dictionary with:
        - tail_threshold: The value at threshold_percentile
        - tail_fraction: Fraction of total recurrence in tail
        - tail_count_fraction: Fraction of fusions in tail
        - tail_mean_ratio: Ratio of tail mean to overall mean
    """
    if len(values) == 0:
        return {
            "tail_threshold": 0.0,
            "tail_fraction": 0.0,
            "tail_count_fraction": 0.0,
            "tail_mean_ratio": 0.0,
        }
    
    threshold = float(np.percentile(values, threshold_percentile))
    tail_mask = values >= threshold
    
    total_sum = float(np.sum(values))
    tail_sum = float(np.sum(values[tail_mask]))
    
    overall_mean = float(np.mean(values))
    tail_mean = float(np.mean(values[tail_mask])) if np.any(tail_mask) else 0.0
    
    tail_fraction = tail_sum / total_sum if total_sum > 0 else 0.0
    tail_count_fraction = float(np.sum(tail_mask)) / len(values)
    tail_mean_ratio = tail_mean / overall_mean if overall_mean > 0 else 0.0
    
    return {
        "tail_threshold": threshold,
        "tail_fraction": tail_fraction,
        "tail_count_fraction": tail_count_fraction,
        "tail_mean_ratio": tail_mean_ratio,
    }


def compute_mock_cosmic_distribution_metrics(
    cosmic_df: pd.DataFrame,
    recurrence_column: str = "recurrence_count",
) -> Dict[str, Any]:
    """
    Compute comprehensive distribution metrics for mock COSMIC validation.
    
    Args:
        cosmic_df: DataFrame with COSMIC fusion data
        recurrence_column: Column name for recurrence counts
    
    Returns:
        Dictionary containing all distribution metrics
    """
    if cosmic_df is None or len(cosmic_df) == 0:
        return {
            "status": "FAILED",
            "error": "Empty or None DataFrame provided",
            "metrics": None,
        }
    
    if recurrence_column not in cosmic_df.columns:
        return {
            "status": "FAILED",
            "error": f"Column '{recurrence_column}' not found in DataFrame",
            "metrics": None,
        }
    
    values = cosmic_df[recurrence_column].values.astype(float)
    
    # Basic statistics
    n = len(values)
    mean_val = float(np.mean(values))
    median_val = float(np.median(values))
    std_val = float(np.std(values))
    min_val = float(np.min(values))
    max_val = float(np.max(values))
    
    # Skewness
    if SCIPY_AVAILABLE:
        skewness = float(skew(values, bias=False))
    else:
        # Numpy fallback for skewness
        m3 = np.mean((values - mean_val) ** 3)
        m2 = np.mean((values - mean_val) ** 2)
        skewness = float(m3 / (m2 ** 1.5)) if m2 > 0 else 0.0
    
    # Kurtosis (excess kurtosis, normal = 0)
    if SCIPY_AVAILABLE:
        kurt = float(kurtosis(values, bias=False))
    else:
        # Numpy fallback for kurtosis
        m4 = np.mean((values - mean_val) ** 4)
        m2 = np.mean((values - mean_val) ** 2)
        kurt = float(m4 / (m2 ** 2) - 3) if m2 > 0 else 0.0
    
    # Zero inflation rate
    zero_count = int(np.sum(values == 0))
    zero_inflation_rate = float(zero_count / n) if n > 0 else 0.0
    
    # Gini coefficient
    gini = compute_gini_coefficient(values)
    
    # Tail heaviness metrics
    tail_metrics = compute_tail_heaviness(values, threshold_percentile=90)
    
    # Coefficient of variation (CV = std/mean)
    cv = float(std_val / mean_val) if mean_val > 0 else 0.0
    
    # Log-scale range (orders of magnitude)
    positive_values = values[values > 0]
    if len(positive_values) > 1:
        log_range = float(np.log10(np.max(positive_values) / np.min(positive_values)))
    else:
        log_range = 0.0
    
    # Quantiles
    quantiles = {
        "q25": float(np.percentile(values, 25)),
        "q50": float(np.percentile(values, 50)),
        "q75": float(np.percentile(values, 75)),
        "q90": float(np.percentile(values, 90)),
        "q95": float(np.percentile(values, 95)),
        "q99": float(np.percentile(values, 99)),
    }
    
    # Realism assessment
    realism_assessment = assess_distribution_realism(
        skewness=skewness,
        kurtosis=kurt,
        gini=gini,
        log_range=log_range,
        tail_fraction=tail_metrics["tail_fraction"],
    )
    
    metrics = {
        "sample_size": n,
        "basic_statistics": {
            "mean": mean_val,
            "median": median_val,
            "std": std_val,
            "min": min_val,
            "max": max_val,
            "coefficient_of_variation": cv,
        },
        "distribution_shape": {
            "skewness": skewness,
            "kurtosis_excess": kurt,
            "log_range_orders_of_magnitude": log_range,
        },
        "inequality_metrics": {
            "gini_coefficient": gini,
            "zero_inflation_rate": zero_inflation_rate,
            "zero_count": zero_count,
        },
        "tail_analysis": tail_metrics,
        "quantiles": quantiles,
        "realism_assessment": realism_assessment,
    }
    
    return {
        "status": "SUCCESS",
        "error": None,
        "metrics": metrics,
    }


def assess_distribution_realism(
    skewness: float,
    kurtosis: float,
    gini: float,
    log_range: float,
    tail_fraction: float,
) -> Dict[str, Any]:
    """
    Assess whether distribution metrics are realistic for biological recurrence data.
    
    Real COSMIC-like distributions typically have:
    - High positive skewness (> 2.0)
    - High excess kurtosis (> 5.0, heavy tails)
    - High Gini coefficient (> 0.6, inequality)
    - Multiple orders of magnitude range (> 2.0)
    - Significant tail concentration (> 0.3)
    
    Args:
        skewness: Distribution skewness
        kurtosis: Excess kurtosis
        gini: Gini coefficient
        log_range: Log10 range (orders of magnitude)
        tail_fraction: Fraction of total in top 10%
    
    Returns:
        Dictionary with assessment results
    """
    checks = []
    
    # Skewness check (expect positive/right skew)
    skewness_pass = skewness > 1.0
    checks.append({
        "metric": "skewness",
        "value": skewness,
        "threshold": "> 1.0",
        "pass": skewness_pass,
        "interpretation": "Right-skewed" if skewness_pass else "Insufficient skew",
    })
    
    # Kurtosis check (expect heavy tails, excess > 3 is typical for power-law)
    kurtosis_pass = kurtosis > 3.0
    checks.append({
        "metric": "kurtosis_excess",
        "value": kurtosis,
        "threshold": "> 3.0",
        "pass": kurtosis_pass,
        "interpretation": "Heavy tails" if kurtosis_pass else "Light tails",
    })
    
    # Gini check (expect inequality in recurrence)
    gini_pass = gini > 0.5
    checks.append({
        "metric": "gini_coefficient",
        "value": gini,
        "threshold": "> 0.5",
        "pass": gini_pass,
        "interpretation": "High inequality (driver dominance)" if gini_pass else "Too uniform",
    })
    
    # Log range check (expect multiple orders of magnitude)
    range_pass = log_range > 1.5
    checks.append({
        "metric": "log_range",
        "value": log_range,
        "threshold": "> 1.5",
        "pass": range_pass,
        "interpretation": "Wide dynamic range" if range_pass else "Narrow range",
    })
    
    # Tail concentration check
    tail_pass = tail_fraction > 0.25
    checks.append({
        "metric": "tail_fraction",
        "value": tail_fraction,
        "threshold": "> 0.25",
        "pass": tail_pass,
        "interpretation": "Concentrated tail" if tail_pass else "Dispersed distribution",
    })
    
    # Overall assessment
    passes = sum(1 for c in checks if c["pass"])
    total = len(checks)
    
    if passes == total:
        overall = "REALISTIC"
        summary = "Mock COSMIC distribution matches expected biological recurrence patterns."
    elif passes >= 3:
        overall = "ACCEPTABLE"
        summary = "Mock COSMIC distribution largely matches expected patterns with minor deviations."
    else:
        overall = "DEGRADED"
        summary = "Mock COSMIC distribution deviates significantly from expected biological patterns."
    
    return {
        "overall_assessment": overall,
        "summary": summary,
        "checks_passed": passes,
        "checks_total": total,
        "individual_checks": checks,
    }


def generate_mock_cosmic_validation_report(
    cosmic_df: pd.DataFrame,
    output_path: Optional[Path] = None,
    recurrence_column: str = "recurrence_count",
    generation_seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Generate comprehensive mock COSMIC validation report.
    
    Args:
        cosmic_df: DataFrame with COSMIC fusion data
        output_path: Optional path to write JSON report
        recurrence_column: Column name for recurrence counts
        generation_seed: Random seed used for mock generation (for provenance)
    
    Returns:
        Complete validation report dictionary
    """
    result = compute_mock_cosmic_distribution_metrics(cosmic_df, recurrence_column)
    
    report = {
        "report_type": "mock_cosmic_distribution_validation",
        "generated_at": datetime.now().isoformat(),
        "validation_status": result["status"],
        "error": result["error"],
        "metrics": result.get("metrics"),
        "provenance": {
            "generation_seed": generation_seed,
            "recurrence_column": recurrence_column,
            "scipy_available": SCIPY_AVAILABLE,
        },
    }
    
    # Add compatibility statement
    if result["status"] == "SUCCESS":
        report["compatibility_statement"] = {
            "purpose": "Mock COSMIC exists for pipeline operability testing",
            "preserves": [
                "Heavy tail recurrence distribution",
                "Driver gene enrichment structure",
                "Non-uniform recurrence behavior",
                "Known fusion pair inclusion",
            ],
            "does_not_replicate": [
                "Real COSMIC cohort structure",
                "Actual variant frequencies",
                "Sample population demographics",
                "Temporal sampling patterns",
            ],
            "appropriate_use": "Pipeline testing and development only",
            "inappropriate_use": "Scientific conclusions about real fusion biology",
        }
    
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        _logger.info(f"Mock COSMIC validation report written to {output_path}")
    
    return report


def validate_mock_cosmic_from_file(
    cosmic_path: Path,
    output_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Validate mock COSMIC from a file path.
    
    Convenience function that loads the mock COSMIC CSV and runs validation.
    
    Args:
        cosmic_path: Path to mock_cosmic_census.csv
        output_dir: Optional directory for output report
    
    Returns:
        Validation report dictionary
    """
    cosmic_path = Path(cosmic_path)
    
    if not cosmic_path.exists():
        return {
            "report_type": "mock_cosmic_distribution_validation",
            "generated_at": datetime.now().isoformat(),
            "validation_status": "FAILED",
            "error": f"File not found: {cosmic_path}",
            "metrics": None,
        }
    
    try:
        cosmic_df = pd.read_csv(cosmic_path)
    except Exception as e:
        return {
            "report_type": "mock_cosmic_distribution_validation",
            "generated_at": datetime.now().isoformat(),
            "validation_status": "FAILED",
            "error": f"Failed to load CSV: {e}",
            "metrics": None,
        }
    
    output_path = None
    if output_dir:
        output_path = Path(output_dir) / "mock_cosmic_distribution_validation.json"
    
    return generate_mock_cosmic_validation_report(
        cosmic_df=cosmic_df,
        output_path=output_path,
        generation_seed=42,  # Default seed used in generate_mock_cosmic
    )


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        cosmic_path = Path(sys.argv[1])
        output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else cosmic_path.parent
    else:
        script_dir = Path(__file__).parent
        cosmic_path = script_dir / "mock_cosmic_census.csv"
        output_dir = script_dir
    
    print(f"Validating mock COSMIC at {cosmic_path}...")
    report = validate_mock_cosmic_from_file(cosmic_path, output_dir)
    
    print(f"\nValidation Status: {report['validation_status']}")
    if report.get("metrics"):
        assessment = report["metrics"]["realism_assessment"]
        print(f"Overall Assessment: {assessment['overall_assessment']}")
        print(f"Checks Passed: {assessment['checks_passed']}/{assessment['checks_total']}")
        print(f"\nSummary: {assessment['summary']}")
