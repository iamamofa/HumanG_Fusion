"""
Week 2: Data Integrity & Statistical Validation — Effect Size Interpretation.

Generates biological interpretation text for correlation effect sizes.
Uses ONLY existing pipeline metrics (no calculations).
"""

from typing import Optional


def generate_correlation_effect_interpretation(
    rho: Optional[float],
    p_value: Optional[float],
    overlap_n: Optional[int],
) -> str:
    """
    Generate biological interpretation text for correlation effect size.
    
    This function reads ONLY existing pipeline metrics and generates interpretation.
    It does NOT perform any statistical calculations.
    
    Args:
        rho: Spearman correlation coefficient (from pipeline).
        p_value: Spearman p-value (from pipeline).
        overlap_n: Number of overlapping fusion pairs (from pipeline).
    
    Returns:
        Biological interpretation paragraph text.
    """
    if rho is None:
        return (
            "Correlation effect size cannot be assessed due to insufficient overlap "
            "or computational limitations."
        )
    
    try:
        abs_rho = abs(float(rho))
    except (TypeError, ValueError):
        return (
            "Correlation effect size cannot be assessed due to invalid correlation value."
        )
    
    # Effect size interpretation
    if abs_rho < 0.3:
        effect_text = (
            "Weak biological agreement. The observed correlation suggests minimal "
            "rank-order consistency with COSMIC reference patterns, indicating "
            "potential cohort-specific effects or methodological differences."
        )
    elif abs_rho < 0.6:
        effect_text = (
            "Moderate biological pattern similarity. The observed correlation suggests "
            "partial rank-order consistency with COSMIC reference patterns, indicating "
            "some shared biological structure but with notable cohort-specific variation."
        )
    elif abs_rho < 0.8:
        effect_text = (
            "Strong shared recurrence structure. The observed correlation suggests "
            "substantial rank-order consistency with COSMIC reference patterns, "
            "indicating strong biological agreement in fusion recurrence ordering."
        )
    else:
        effect_text = (
            "Very strong rank preservation. The observed correlation suggests "
            "excellent rank-order consistency with COSMIC reference patterns, "
            "indicating highly consistent biological recurrence structure."
        )
    
    # Add uncertainty qualifiers if needed
    uncertainty_parts = []
    
    if overlap_n is not None and overlap_n < 10:
        uncertainty_parts.append(
            f"Limited overlap ({overlap_n} pairs) introduces uncertainty in "
            "correlation estimates, and results should be interpreted with caution."
        )
    
    if p_value is not None:
        try:
            p_val = float(p_value)
            if p_val > 0.05:
                uncertainty_parts.append(
                    "Statistical significance thresholds are not met within conventional "
                    "confidence intervals, suggesting the correlation signal should be "
                    "interpreted cautiously."
                )
        except (TypeError, ValueError):
            pass
    
    # Combine effect size and uncertainty
    if uncertainty_parts:
        uncertainty_text = " ".join(uncertainty_parts)
        return f"{effect_text} {uncertainty_text}"
    
    return effect_text
