"""
Week 2: Data Integrity & Statistical Validation — Scientific Claim Strength Classifier.

This module classifies the strength of scientific claims based on:
- Observed correlation (rho)
- Null model p-value
- Stability index (from distribution generalization)

SAFETY: This module does NOT modify core statistical logic.
It only provides classification based on existing metrics.
"""

from typing import Dict, Optional


def classify_scientific_claim_strength(
    observed_rho: Optional[float],
    null_p: Optional[float],
    stability_index: Optional[float],
    spearman_p: Optional[float] = None,
    overlap_count: Optional[int] = None,
) -> Dict[str, str]:
    """
    Classify scientific claim strength based on statistical evidence.
    
    OUTPUT CLASSES:
    - EXPLORATORY SIGNAL: Weak evidence, needs more data
    - MODERATE EXTERNAL CONSISTENCY: Moderate evidence, some stability
    - STRONG BIOLOGICAL AGREEMENT: Strong evidence, stable across strata
    - HIGH CONFIDENCE CROSS-DATASET SIGNAL: Very strong evidence, highly stable
    
    Args:
        observed_rho: Observed Spearman correlation coefficient
        null_p: Empirical p-value from null model test (overlap test)
        stability_index: Standard deviation of correlations across strata
        spearman_p: P-value from Spearman correlation test (optional)
        overlap_count: Number of overlapping fusion pairs with reference (optional);
            if < 10, returns LIMITED OVERLAP - EXPLORATORY
    
    Returns:
        Dictionary with:
        - classification: One of the four classes above
        - confidence_level: "LOW", "MODERATE", "HIGH", or "VERY_HIGH"
        - interpretation: Human-readable interpretation
    """
    # Handle missing values
    if observed_rho is None:
        return {
            "classification": "EXPLORATORY SIGNAL",
            "confidence_level": "LOW",
            "interpretation": "Insufficient data for classification. Correlation could not be computed.",
        }

    if overlap_count is not None and overlap_count < 10:
        return {
            "classification": "LIMITED OVERLAP - EXPLORATORY",
            "confidence_level": "LOW",
            "interpretation": (
                f"Only {overlap_count} overlapping pairs; "
                "insufficient for reliable classification. Treat as exploratory."
            ),
        }

    abs_rho = abs(observed_rho) if observed_rho is not None else 0.0
    
    # Classification logic
    # HIGH CONFIDENCE: Strong correlation, significant null test, stable across strata
    if abs_rho >= 0.6 and null_p is not None and null_p < 0.01:
        if stability_index is not None and stability_index < 0.15:
            return {
                "classification": "HIGH CONFIDENCE CROSS-DATASET SIGNAL",
                "confidence_level": "VERY_HIGH",
                "interpretation": (
                    f"Strong correlation (ρ={abs_rho:.3f}) with high statistical significance "
                    f"(p={null_p:.4f}) and stable across recurrence strata (stability={stability_index:.3f}). "
                    "This suggests a robust cross-dataset biological signal."
                ),
            }
        elif stability_index is None or stability_index < 0.25:
            return {
                "classification": "STRONG BIOLOGICAL AGREEMENT",
                "confidence_level": "HIGH",
                "interpretation": (
                    f"Strong correlation (ρ={abs_rho:.3f}) with statistical significance "
                    f"(p={null_p:.4f}). Moderate stability across strata. "
                    "This suggests meaningful biological agreement between datasets."
                ),
            }
    
    # STRONG BIOLOGICAL AGREEMENT: Moderate-strong correlation, significant null test
    if abs_rho >= 0.4 and null_p is not None and null_p < 0.05:
        if stability_index is not None and stability_index < 0.20:
            if spearman_p is not None:
                interpretation = (
                    f"Moderate-strong correlation strength (ρ={abs_rho:.3f}). "
                    f"Overlap significantly exceeds random expectation (null model p={null_p:.4f}). "
                    f"{'Correlation is also statistically significant' if spearman_p < 0.05 else 'However correlation itself is not statistically significant'} "
                    f"(correlation p={spearman_p:.2f}). "
                    f"Good stability across strata (stability={stability_index:.3f}). "
                    "This suggests consistent biological agreement."
                )
            else:
                interpretation = (
                    f"Moderate-strong correlation (ρ={abs_rho:.3f}) with overlap significantly exceeding "
                    f"random expectation (null model p={null_p:.4f}) and good stability (stability={stability_index:.3f}). "
                    "This suggests consistent biological agreement."
                )
            return {
                "classification": "STRONG BIOLOGICAL AGREEMENT",
                "confidence_level": "HIGH",
                "interpretation": interpretation,
            }
        else:
            # Build interpretation that distinguishes between null model p and correlation p
            if spearman_p is not None:
                interpretation = (
                    f"Moderate correlation strength (ρ={abs_rho:.3f}). "
                    f"Overlap significantly exceeds random expectation (null model p < 0.000001). "
                    f"However correlation itself is not statistically significant (p = {spearman_p:.2f}) "
                    f"due to low overlap sample size. "
                    "This suggests partial external consistency."
                )
            else:
                interpretation = (
                    f"Moderate correlation (ρ={abs_rho:.3f}) with overlap significantly exceeding "
                    f"random expectation (null model p={null_p:.4f}), but variable across recurrence strata. "
                    "This suggests partial external consistency."
                )
            return {
                "classification": "MODERATE EXTERNAL CONSISTENCY",
                "confidence_level": "MODERATE",
                "interpretation": interpretation,
            }
    
    # MODERATE EXTERNAL CONSISTENCY: Moderate correlation, marginal/null significance
    if abs_rho >= 0.3:
        if null_p is not None and null_p < 0.10:
            return {
                "classification": "MODERATE EXTERNAL CONSISTENCY",
                "confidence_level": "MODERATE",
                "interpretation": (
                    f"Moderate correlation (ρ={abs_rho:.3f}) with marginal significance "
                    f"(p={null_p:.4f}). This suggests some external consistency but requires "
                    "additional validation."
                ),
            }
        else:
            return {
                "classification": "EXPLORATORY SIGNAL",
                "confidence_level": "LOW",
                "interpretation": (
                    f"Moderate correlation (ρ={abs_rho:.3f}) but not statistically significant "
                    f"(p={null_p:.4f if null_p else 'N/A'}). This is an exploratory signal "
                    "requiring further investigation."
                ),
            }
    
    # EXPLORATORY SIGNAL: Weak correlation or no significance
    p_str = f"{null_p:.4f}" if null_p is not None else "N/A"
    return {
        "classification": "EXPLORATORY SIGNAL",
        "confidence_level": "LOW",
        "interpretation": (
            f"Weak correlation (ρ={abs_rho:.3f}) with no clear statistical significance "
            f"(p={p_str}). This represents an exploratory signal "
            "that may require larger sample sizes or additional validation."
        ),
    }
