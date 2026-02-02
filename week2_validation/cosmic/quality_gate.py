"""
COSMIC Quality Gate - Validation Score and Classification
Week 2: Data Integrity & Statistical Validation

This module computes a composite COSMIC validation score and classifies the
biological agreement strength based on multiple statistical metrics.

WHAT DOES THIS MODULE DO?
This module combines multiple statistical signals into a single quality score
and provides a classification of biological agreement strength:
- STRONG BIOLOGICAL AGREEMENT
- MODERATE AGREEMENT
- WEAK SIGNAL
- NO MEANINGFUL COSMIC AGREEMENT

DESIGN PRINCIPLES:
1. Quality gates provide interpretable classifications while preserving
   all underlying statistical metrics for detailed analysis.
2. Full transparency: every component score is reported individually.
3. Score component breakdown enables understanding of which aspects
   contribute to the final classification.

SCORE COMPONENT STRUCTURE:
- correlation_component: Spearman rho strength contribution
- significance_component: Spearman p-value contribution
- enrichment_component: Hypergeometric enrichment contribution
- negative_control_component: Real vs shuffled difference contribution
- overlap_component: Observed vs expected overlap contribution
"""

import logging
from typing import Dict, Optional, Any

_logger = logging.getLogger(__name__)


def compute_cosmic_validation_score(
    spearman_rho: Optional[float],
    spearman_p_value: Optional[float],
    enrichment_p_value: Optional[float],
    negative_control_rho: Optional[float],
    observed_overlap: Optional[int],
    expected_overlap_random: Optional[float],
) -> Dict[str, Any]:
    """
    Compute COSMIC validation score and classification.
    
    The validation score combines:
    1. Spearman correlation strength (rho)
    2. Statistical significance (p-values)
    3. Enrichment signal (hypergeometric p-value)
    4. Negative control validation (real vs shuffled)
    5. Overlap ratio (observed vs expected)
    
    Args:
        spearman_rho: Spearman correlation coefficient
        spearman_p_value: p-value for Spearman test
        enrichment_p_value: p-value from hypergeometric enrichment test
        negative_control_rho: Mean correlation from shuffled negative control
        observed_overlap: Observed overlap in top N fusions
        expected_overlap_random: Expected overlap by chance
    
    Returns:
        Dictionary with:
        - cosmic_validation_score: Composite score (0.0 to 1.0)
        - cosmic_validation_classification: Classification string
        - score_components: Legacy breakdown (backward compatible)
        - score_component_breakdown: Detailed transparency breakdown
    """
    score_components = {}
    component_breakdown = {}
    total_score = 0.0
    max_score = 0.0
    
    # Component 1: Spearman correlation strength (0-0.3 points)
    correlation_component = 0.0
    correlation_max = 0.3
    if spearman_rho is not None:
        abs_rho = abs(spearman_rho)
        if abs_rho >= 0.7:
            rho_score = 0.3
        elif abs_rho >= 0.5:
            rho_score = 0.2
        elif abs_rho >= 0.3:
            rho_score = 0.1
        else:
            rho_score = 0.0
        score_components["spearman_strength"] = rho_score
        correlation_component = rho_score
        total_score += rho_score
    max_score += correlation_max
    
    component_breakdown["correlation_component"] = {
        "raw_score": correlation_component,
        "max_possible": correlation_max,
        "normalized": correlation_component / correlation_max if correlation_max > 0 else 0.0,
        "input_value": spearman_rho,
        "weight_fraction": correlation_max / 1.0,
        "interpretation": _interpret_correlation_score(spearman_rho),
    }
    
    # Component 2: Spearman significance (0-0.2 points)
    significance_component = 0.0
    significance_max = 0.2
    if spearman_p_value is not None:
        if spearman_p_value < 0.001:
            p_score = 0.2
        elif spearman_p_value < 0.01:
            p_score = 0.15
        elif spearman_p_value < 0.05:
            p_score = 0.1
        else:
            p_score = 0.0
        score_components["spearman_significance"] = p_score
        significance_component = p_score
        total_score += p_score
    max_score += significance_max
    
    # Combine correlation and significance for total correlation component
    combined_correlation = correlation_component + significance_component
    
    # Component 3: Enrichment significance (0-0.2 points)
    enrichment_component = 0.0
    enrichment_max = 0.2
    if enrichment_p_value is not None:
        if enrichment_p_value < 0.001:
            enrich_score = 0.2
        elif enrichment_p_value < 0.01:
            enrich_score = 0.15
        elif enrichment_p_value < 0.05:
            enrich_score = 0.1
        else:
            enrich_score = 0.0
        score_components["enrichment_significance"] = enrich_score
        enrichment_component = enrich_score
        total_score += enrich_score
    max_score += enrichment_max
    
    component_breakdown["enrichment_component"] = {
        "raw_score": enrichment_component,
        "max_possible": enrichment_max,
        "normalized": enrichment_component / enrichment_max if enrichment_max > 0 else 0.0,
        "input_value": enrichment_p_value,
        "weight_fraction": enrichment_max / 1.0,
        "interpretation": _interpret_enrichment_score(enrichment_p_value),
    }
    
    # Component 4: Negative control validation (0-0.2 points)
    # Real correlation should be much stronger than shuffled
    negative_control_component = 0.0
    negative_control_max = 0.2
    delta_value = None
    if spearman_rho is not None and negative_control_rho is not None:
        real_abs = abs(spearman_rho)
        shuffled_abs = abs(negative_control_rho)
        if real_abs > 0 and shuffled_abs >= 0:
            delta = real_abs - shuffled_abs
            delta_value = delta
            if delta > 0.5:
                control_score = 0.2
            elif delta > 0.3:
                control_score = 0.15
            elif delta > 0.1:
                control_score = 0.1
            else:
                control_score = 0.0
            score_components["negative_control_delta"] = control_score
            negative_control_component = control_score
            total_score += control_score
    max_score += negative_control_max
    
    component_breakdown["negative_control_component"] = {
        "raw_score": negative_control_component,
        "max_possible": negative_control_max,
        "normalized": negative_control_component / negative_control_max if negative_control_max > 0 else 0.0,
        "input_value": {"real_rho": spearman_rho, "shuffled_rho": negative_control_rho, "delta": delta_value},
        "weight_fraction": negative_control_max / 1.0,
        "interpretation": _interpret_negative_control_score(delta_value),
    }
    
    # Component 5: Overlap enrichment ratio (0-0.1 points)
    overlap_component = 0.0
    overlap_max = 0.1
    enrichment_ratio = None
    if observed_overlap is not None and expected_overlap_random is not None:
        if expected_overlap_random > 0:
            enrichment_ratio = observed_overlap / expected_overlap_random
            if enrichment_ratio >= 3.0:
                overlap_score = 0.1
            elif enrichment_ratio >= 2.0:
                overlap_score = 0.075
            elif enrichment_ratio >= 1.5:
                overlap_score = 0.05
            else:
                overlap_score = 0.0
            score_components["overlap_enrichment"] = overlap_score
            overlap_component = overlap_score
            total_score += overlap_score
    max_score += overlap_max
    
    component_breakdown["overlap_component"] = {
        "raw_score": overlap_component,
        "max_possible": overlap_max,
        "normalized": overlap_component / overlap_max if overlap_max > 0 else 0.0,
        "input_value": {"observed": observed_overlap, "expected": expected_overlap_random, "ratio": enrichment_ratio},
        "weight_fraction": overlap_max / 1.0,
        "interpretation": _interpret_overlap_score(enrichment_ratio),
    }
    
    # Normalize score to 0-1 range
    normalized_score = total_score / max_score if max_score > 0 else 0.0
    
    # Classification based on score
    if normalized_score >= 0.7:
        classification = "STRONG BIOLOGICAL AGREEMENT"
    elif normalized_score >= 0.5:
        classification = "MODERATE AGREEMENT"
    elif normalized_score >= 0.3:
        classification = "WEAK SIGNAL"
    else:
        classification = "NO MEANINGFUL COSMIC AGREEMENT"
    
    # Verify score breakdown sums correctly
    breakdown_sum = (
        correlation_component + 
        significance_component + 
        enrichment_component + 
        negative_control_component + 
        overlap_component
    )
    
    return {
        "cosmic_validation_score": float(normalized_score),
        "cosmic_validation_classification": classification,
        "score_components": score_components,  # Legacy format for backward compatibility
        "score_component_breakdown": {
            "correlation_component": combined_correlation,
            "enrichment_component": enrichment_component,
            "negative_control_component": negative_control_component,
            "overlap_component": overlap_component,
            "total_raw_score": total_score,
            "max_possible_score": max_score,
            "breakdown_sum_check": abs(breakdown_sum - total_score) < 0.0001,
            "detailed_components": component_breakdown,
        },
    }


def _interpret_correlation_score(rho: Optional[float]) -> str:
    """Provide interpretation for Spearman rho score."""
    if rho is None:
        return "Unable to compute - insufficient data"
    abs_rho = abs(rho)
    if abs_rho >= 0.7:
        return "Strong rank-order agreement with COSMIC"
    elif abs_rho >= 0.5:
        return "Moderate rank-order agreement with COSMIC"
    elif abs_rho >= 0.3:
        return "Weak rank-order agreement with COSMIC"
    else:
        return "No meaningful rank-order agreement"


def _interpret_enrichment_score(p_value: Optional[float]) -> str:
    """Provide interpretation for enrichment p-value."""
    if p_value is None:
        return "Unable to compute - insufficient data"
    if p_value < 0.001:
        return "Highly significant enrichment (p < 0.001)"
    elif p_value < 0.01:
        return "Significant enrichment (p < 0.01)"
    elif p_value < 0.05:
        return "Marginally significant enrichment (p < 0.05)"
    else:
        return "No significant enrichment detected"


def _interpret_negative_control_score(delta: Optional[float]) -> str:
    """Provide interpretation for negative control delta."""
    if delta is None:
        return "Unable to compute - insufficient data"
    if delta > 0.5:
        return "Strong evidence of real biological signal (real >> shuffled)"
    elif delta > 0.3:
        return "Good evidence of biological signal (real > shuffled)"
    elif delta > 0.1:
        return "Weak evidence of biological signal"
    else:
        return "Cannot distinguish from random (real ≈ shuffled)"


def _interpret_overlap_score(ratio: Optional[float]) -> str:
    """Provide interpretation for overlap enrichment ratio."""
    if ratio is None:
        return "Unable to compute - insufficient data"
    if ratio >= 3.0:
        return "Strong overlap enrichment (>3x expected)"
    elif ratio >= 2.0:
        return "Moderate overlap enrichment (2-3x expected)"
    elif ratio >= 1.5:
        return "Weak overlap enrichment (1.5-2x expected)"
    else:
        return "No meaningful overlap enrichment"
