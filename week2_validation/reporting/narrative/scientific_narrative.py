"""
Week 2: Data Integrity & Statistical Validation — Scientific Narrative Engine.

Generates integrated scientific narrative paragraphs from pipeline metrics.
This module reads ONLY existing pipeline outputs and generates interpretation text.
No statistical calculations are performed here.
"""

from typing import Any, Dict, Optional


def generate_integrated_scientific_narrative(results_dict: Dict[str, Any]) -> str:
    """
    Generate ONE integrated scientific story paragraph using ALL metrics together.
    
    This function reads ONLY existing pipeline metrics and generates a cohesive
    scientific narrative. It does NOT perform any statistical calculations.
    
    Args:
        results_dict: Dictionary containing existing pipeline metrics:
            - skewness: float (optional)
            - mean: float (optional)
            - median: float (optional)
            - ks_p_value: float (optional)
            - ad_p_value: float (optional)
            - benford_applicable: bool (optional)
            - scale_span_orders: float (optional)
            - cosmic_overlap: int (optional)
            - spearman_rho: float (optional)
            - spearman_p_value: float (optional)
            - cosmic_quality_score: float (optional)
    
    Returns:
        Integrated scientific narrative paragraph (journal tone, evidence + interpretation + implication).
    """
    # Extract metrics (handle None values safely)
    skewness = results_dict.get("skewness")
    mean = results_dict.get("mean")
    median = results_dict.get("median")
    ks_p_value = results_dict.get("ks_p_value")
    ad_p_value = results_dict.get("ad_p_value")
    benford_applicable = results_dict.get("benford_applicable")
    scale_span = results_dict.get("scale_span_orders")
    cosmic_overlap = results_dict.get("cosmic_overlap")
    spearman_rho = results_dict.get("spearman_rho")
    spearman_p_value = results_dict.get("spearman_p_value")
    cosmic_quality_score = results_dict.get("cosmic_quality_score")
    
    narrative_parts = []
    
    # Distribution shape assessment
    if skewness is not None:
        try:
            abs_skew = abs(float(skewness))
            if abs_skew < 0.2:
                narrative_parts.append(
                    "The protein length distribution is approximately symmetric, "
                    "consistent with a well-behaved underlying process."
                )
            elif abs_skew < 0.8:
                narrative_parts.append(
                    "The protein length distribution exhibits mild asymmetry, "
                    "suggesting potential underlying biological or methodological factors."
                )
            else:
                narrative_parts.append(
                    "The protein length distribution shows strong skewness, "
                    "indicating a departure from symmetric models and suggesting "
                    "potential heavy-tailed or power-law-like behavior."
                )
        except (TypeError, ValueError):
            pass
    
    # Mean-median relationship
    if mean is not None and median is not None:
        try:
            mean_val = float(mean)
            median_val = float(median)
            if mean_val != 0:
                divergence_pct = abs(mean_val - median_val) / abs(mean_val) * 100
                if divergence_pct < 5:
                    narrative_parts.append(
                        "Mean and median values are in close agreement, "
                        "supporting the assessment of distribution symmetry."
                    )
                else:
                    narrative_parts.append(
                        "Mean-median divergence suggests the presence of outliers "
                        "or asymmetric tail behavior in the distribution."
                    )
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    
    # Log-normality assessment
    log_normal_compatible = False
    if ks_p_value is not None:
        try:
            if float(ks_p_value) > 0.05:
                log_normal_compatible = True
        except (TypeError, ValueError):
            pass
    
    if ad_p_value is not None:
        try:
            if float(ad_p_value) > 0.05:
                log_normal_compatible = True
        except (TypeError, ValueError):
            pass
    
    if log_normal_compatible:
        narrative_parts.append(
            "Distribution tests are consistent with log-normal compatibility, "
            "suggesting multiplicative processes may underlie the observed patterns."
        )
    elif ks_p_value is not None or ad_p_value is not None:
        narrative_parts.append(
            "Distribution tests indicate deviation from log-normal behavior, "
            "suggesting alternative generative models may be more appropriate."
        )
    
    # Benford assessment
    if benford_applicable is False:
        if scale_span is not None:
            try:
                span_val = float(scale_span)
                if span_val < 2.0:
                    narrative_parts.append(
                        f"The limited scale span ({span_val:.2f} orders of magnitude) "
                        "precludes reliable Benford analysis, indicating the dataset "
                        "operates within a constrained range."
                    )
            except (TypeError, ValueError):
                narrative_parts.append(
                    "Benford analysis is not applicable due to limited scale span."
                )
    
    # COSMIC cross-validation
    if spearman_rho is not None:
        try:
            rho_val = float(spearman_rho)
            abs_rho = abs(rho_val)
            
            if cosmic_overlap is not None and cosmic_overlap < 10:
                narrative_parts.append(
                    f"Cross-validation with COSMIC reference data shows "
                    f"{'positive' if rho_val > 0 else 'negative'} rank-order correlation "
                    f"(ρ={rho_val:.3f}), though limited overlap ({cosmic_overlap} pairs) "
                    "suggests caution in interpretation."
                )
            elif abs_rho > 0.7:
                narrative_parts.append(
                    f"Strong rank-order correlation with COSMIC reference data "
                    f"(ρ={rho_val:.3f}) suggests biological consistency with established "
                    "cancer fusion patterns."
                )
            elif abs_rho > 0.3:
                narrative_parts.append(
                    f"Moderate rank-order correlation with COSMIC reference data "
                    f"(ρ={rho_val:.3f}) suggests partial biological consistency, "
                    "though cohort differences may contribute to observed variation."
                )
            else:
                narrative_parts.append(
                    f"Weak correlation with COSMIC reference data (ρ={rho_val:.3f}) "
                    "suggests potential cohort-specific effects or methodological differences."
                )
        except (TypeError, ValueError):
            pass
    
    # Quality gate assessment
    if cosmic_quality_score is not None:
        try:
            score = float(cosmic_quality_score)
            if score >= 0.7:
                narrative_parts.append(
                    "COSMIC validation metrics indicate strong biological agreement, "
                    "supporting the biological plausibility of observed fusion patterns."
                )
            elif score >= 0.5:
                narrative_parts.append(
                    "COSMIC validation metrics indicate moderate biological agreement, "
                    "suggesting plausible but not definitive biological consistency."
                )
            else:
                narrative_parts.append(
                    "COSMIC validation metrics indicate weak biological agreement, "
                    "suggesting cohort-specific or methodological factors may dominate."
                )
        except (TypeError, ValueError):
            pass
    
    # Combine into integrated narrative
    if narrative_parts:
        # Join with proper flow and transitions
        if len(narrative_parts) == 1:
            return narrative_parts[0]
        elif len(narrative_parts) == 2:
            return f"{narrative_parts[0]} {narrative_parts[1]}"
        else:
            # Multiple parts: use transitions
            result = narrative_parts[0]
            for part in narrative_parts[1:-1]:
                result += f" {part}"
            result += f" {narrative_parts[-1]}"
            return result
    
    # Default if no metrics available
    return (
        "Scientific interpretation requires additional diagnostic metrics. "
        "Statistical analysis was performed, but integrated narrative generation "
        "requires distribution, Benford, and COSMIC validation results."
    )


def generate_confidence_statement(results_dict: Dict[str, Any]) -> str:
    """
    Generate confidence and uncertainty statement from pipeline metrics.
    
    This function reads ONLY existing pipeline metrics and generates confidence
    assessment text. It does NOT perform any statistical calculations.
    
    Args:
        results_dict: Dictionary containing existing pipeline metrics.
    
    Returns:
        Confidence and uncertainty statement paragraph.
    """
    confidence_parts = []
    
    # Extract metrics
    cosmic_overlap = results_dict.get("cosmic_overlap")
    spearman_p_value = results_dict.get("spearman_p_value")
    ks_p_value = results_dict.get("ks_p_value")
    ad_p_value = results_dict.get("ad_p_value")
    
    # Sample size risk
    if cosmic_overlap is not None:
        try:
            overlap_count = int(cosmic_overlap)
            if overlap_count < 10:
                confidence_parts.append(
                    f"Limited overlap with COSMIC reference ({overlap_count} pairs) "
                    "introduces instability risk in correlation estimates, and "
                    "results should be interpreted with caution."
                )
        except (TypeError, ValueError):
            pass
    
    # Correlation confidence
    if spearman_p_value is not None:
        try:
            p_val = float(spearman_p_value)
            if p_val > 0.05:
                confidence_parts.append(
                    "Correlation signal should be interpreted cautiously, as "
                    "statistical significance thresholds are not met within "
                    "conventional confidence intervals."
                )
            elif p_val > 0.01:
                confidence_parts.append(
                    "Correlation signal is statistically significant but moderate, "
                    "suggesting biological signal within statistical confidence."
                )
        except (TypeError, ValueError):
            pass
    
    # Distribution test confidence
    distribution_rejected = False
    if ks_p_value is not None:
        try:
            if float(ks_p_value) <= 0.05:
                distribution_rejected = True
        except (TypeError, ValueError):
            pass
    
    if ad_p_value is not None:
        try:
            if float(ad_p_value) <= 0.05:
                distribution_rejected = True
        except (TypeError, ValueError):
            pass
    
    if distribution_rejected:
        confidence_parts.append(
            "Distribution tests indicate deviation from theoretical models, "
            "suggesting that parametric assumptions may not hold and alternative "
            "modeling approaches should be considered."
        )
    
    # Combine confidence statements
    if confidence_parts:
        if len(confidence_parts) == 1:
            return confidence_parts[0]
        else:
            return " ".join(confidence_parts)
    
    # Default confidence statement
    return (
        "Statistical confidence is assessed through multiple diagnostic tests. "
        "Results should be interpreted within the context of sample size, "
        "methodological constraints, and biological variability."
    )


def generate_limitations_section(results_dict: Dict[str, Any]) -> str:
    """
    Generate limitations and scope statement.
    
    This function generates standard limitations text that applies to all analyses.
    It does NOT read metrics (limitations are universal).
    
    Args:
        results_dict: Dictionary (unused, kept for API consistency).
    
    Returns:
        Limitations and scope statement paragraph.
    """
    limitations = [
        "Statistical Limitations: Distribution tests are model-dependent and assess "
        "compatibility with specific theoretical distributions rather than proving "
        "biological truth. Test results indicate consistency or deviation from models, "
        "not absolute data validity.",
        
        "COSMIC Limitations: The COSMIC reference database is not population-representative "
        "and reflects inherent biases including tumor sampling bias, detection technology "
        "bias, cohort representation bias, and publication bias. Comparisons with COSMIC "
        "evaluate plausibility and consistency, not exact biological equivalence.",
        
        "Biological Limitations: Statistical validation cannot prove biological causality "
        "or establish biological truth. These analyses detect anomalies, assess distribution "
        "properties, and evaluate consistency with reference data, but cannot establish "
        "biological mechanisms or causal relationships.",
        
        "Data Integrity Scope: This validation layer detects statistical anomalies, "
        "distribution properties, and cross-validation consistency. It does not validate "
        "biological truth, establish causality, or prove data authenticity beyond "
        "statistical pattern analysis."
    ]
    
    return " ".join(limitations)


def extract_narrative_metrics(
    status_data: Dict[str, Any],
    diagnostic_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Extract metrics from status and diagnostic data for narrative generation.
    
    This function reads ONLY existing values from pipeline outputs.
    No calculations are performed.
    
    Args:
        status_data: Status JSON data.
        diagnostic_data: Diagnostic results JSON data (optional).
    
    Returns:
        Dictionary of metrics for narrative generation.
    """
    metrics = {}
    
    # Extract from diagnostic results if available
    if diagnostic_data:
        dist = diagnostic_data.get("distribution", {})
        log_norm = diagnostic_data.get("log_normality", {})
        benford = diagnostic_data.get("benford", {})
        cosmic = diagnostic_data.get("cosmic", {})
        
        # Distribution metrics
        if "mean" in dist:
            metrics["mean"] = dist["mean"]
        if "median" in dist:
            metrics["median"] = dist["median"]
        if "skewness" in dist:
            metrics["skewness"] = dist["skewness"]
        
        # Log-normality (note: p-values may not be directly available)
        # Placeholder for future enhancement if p-values are added
        
        # Benford metrics
        if "applicability" in benford:
            metrics["benford_applicable"] = benford["applicability"]
        if "scale_span_orders_of_magnitude" in benford:
            metrics["scale_span_orders"] = benford["scale_span_orders_of_magnitude"]
        
        # COSMIC metrics
        if "overlap_count" in cosmic:
            metrics["cosmic_overlap"] = cosmic["overlap_count"]
        if "spearman_rho" in cosmic:
            metrics["spearman_rho"] = cosmic["spearman_rho"]
        if "spearman_p_value" in cosmic:
            metrics["spearman_p_value"] = cosmic["spearman_p_value"]
        if "cosmic_validation_score" in cosmic:
            metrics["cosmic_quality_score"] = cosmic["cosmic_validation_score"]
    
    # Also check status data for skewness (may be stored there)
    dq = status_data.get("data_quality", {})
    if "skewness" in dq and "skewness" not in metrics:
        metrics["skewness"] = dq["skewness"]
    
    return metrics
