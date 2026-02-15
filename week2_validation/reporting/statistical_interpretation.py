"""
Week 2: Data Integrity & Statistical Validation — Statistical Interpretation Engine.

Generates rule-based interpretation of statistical metrics for PDF reports.
This module reads ONLY existing pipeline metrics and generates narrative text.
No statistical calculations are performed here.
"""

from typing import Any, Dict, Optional


def generate_distribution_interpretation(metrics_dict: Dict[str, Any]) -> str:
    """
    Generate statistical interpretation paragraph from existing pipeline metrics.
    
    This function reads ONLY existing values from the pipeline and generates
    interpretation text. It does NOT perform any statistical calculations.
    
    Args:
        metrics_dict: Dictionary containing existing pipeline metrics:
            - mean: float (optional)
            - median: float (optional)
            - std: float (optional)
            - skewness: float (optional)
            - ks_p_value: float (optional)
            - ad_p_value: float (optional)
            - benford_applicable: bool (optional)
            - scale_span_orders_of_magnitude: float (optional)
    
    Returns:
        Full paragraph text with statistical interpretation (professional scientific tone).
    """
    interpretation_parts = []
    
    # Extract metrics (handle None values safely)
    mean = metrics_dict.get("mean")
    median = metrics_dict.get("median")
    std = metrics_dict.get("std")
    skewness = metrics_dict.get("skewness")
    ks_p_value = metrics_dict.get("ks_p_value")
    ad_p_value = metrics_dict.get("ad_p_value")
    benford_applicable = metrics_dict.get("benford_applicable")
    scale_span = metrics_dict.get("scale_span_orders_of_magnitude")
    
    # Skewness Interpretation
    if skewness is not None:
        try:
            abs_skew = abs(float(skewness))
            if abs_skew < 0.2:
                interpretation_parts.append(
                    "Distribution is approximately symmetric."
                )
            elif abs_skew < 0.8:
                interpretation_parts.append(
                    "Distribution shows mild asymmetry."
                )
            else:
                interpretation_parts.append(
                    "Distribution shows strong skewness."
                )
        except (TypeError, ValueError):
            pass  # Skip if skewness cannot be interpreted
    
    # Mean vs Median Interpretation
    if mean is not None and median is not None:
        try:
            mean_val = float(mean)
            median_val = float(median)
            if mean_val != 0:
                divergence_pct = abs(mean_val - median_val) / abs(mean_val) * 100
                if divergence_pct < 5:
                    interpretation_parts.append(
                        "Mean and median agreement supports distribution symmetry."
                    )
                else:
                    interpretation_parts.append(
                        "Mean-median divergence suggests skew or outliers."
                    )
        except (TypeError, ValueError, ZeroDivisionError):
            pass  # Skip if values cannot be compared
    
    # KS Test Interpretation
    if ks_p_value is not None:
        try:
            p_val = float(ks_p_value)
            if p_val > 0.05:
                interpretation_parts.append(
                    "Kolmogorov-Smirnov test cannot reject log-normal compatibility."
                )
            else:
                interpretation_parts.append(
                    "Kolmogorov-Smirnov test indicates distribution deviates from log-normal behavior."
                )
        except (TypeError, ValueError):
            pass  # Skip if p-value cannot be interpreted
    
    # AD Test Interpretation
    if ad_p_value is not None:
        try:
            p_val = float(ad_p_value)
            if p_val > 0.05:
                interpretation_parts.append(
                    "Anderson-Darling test cannot reject log-normal compatibility."
                )
            else:
                interpretation_parts.append(
                    "Anderson-Darling test indicates distribution deviates from log-normal behavior."
                )
        except (TypeError, ValueError):
            pass  # Skip if p-value cannot be interpreted
    
    # Benford Interpretation
    if benford_applicable is False:
        if scale_span is not None:
            try:
                span_val = float(scale_span)
                if span_val < 2.0:
                    interpretation_parts.append(
                        f"Benford analysis not applicable due to limited scale span ({span_val:.2f} orders of magnitude)."
                    )
            except (TypeError, ValueError):
                interpretation_parts.append(
                    "Benford analysis not applicable due to limited scale span."
                )
        else:
            interpretation_parts.append(
                "Benford analysis not applicable due to limited scale span."
            )
    
    # Combine into paragraph
    if interpretation_parts:
        # Join with proper punctuation and flow
        if len(interpretation_parts) == 1:
            return interpretation_parts[0]
        elif len(interpretation_parts) == 2:
            return f"{interpretation_parts[0]} {interpretation_parts[1]}"
        else:
            # Multiple parts: use commas and "and" for final item
            result = ", ".join(interpretation_parts[:-1])
            result += f", and {interpretation_parts[-1]}"
            return result
    
    # Default if no metrics available
    return "Statistical interpretation requires additional diagnostic metrics."


def extract_interpretation_metrics(
    status_data: Dict[str, Any],
    diagnostic_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Extract metrics from status and diagnostic data for interpretation.
    
    This function reads ONLY existing values from pipeline outputs.
    No calculations are performed.
    
    Args:
        status_data: Status JSON data.
        diagnostic_data: Diagnostic results JSON data (optional).
    
    Returns:
        Dictionary of metrics for interpretation.
    """
    metrics = {}
    
    # Extract from diagnostic results if available
    if diagnostic_data:
        dist = diagnostic_data.get("distribution", {})
        log_norm = diagnostic_data.get("log_normality", {})
        benford = diagnostic_data.get("benford", {})
        
        # Distribution metrics
        if "mean" in dist:
            metrics["mean"] = dist["mean"]
        if "median" in dist:
            metrics["median"] = dist["median"]
        if "std" in dist:
            metrics["std"] = dist["std"]
        if "skewness" in dist:
            metrics["skewness"] = dist["skewness"]
        
        # Log-normality test p-values (if available)
        # Note: KS and AD tests don't directly provide p-values in current output
        # This is a placeholder for future enhancement
        if "ks_statistic" in log_norm:
            # KS statistic available but p-value may not be
            pass
        if "ad_statistic" in log_norm:
            # AD statistic available but p-value may not be
            pass
        
        # Benford metrics
        if "applicability" in benford:
            metrics["benford_applicable"] = benford["applicability"]
        if "scale_span_orders_of_magnitude" in benford:
            metrics["scale_span_orders_of_magnitude"] = benford["scale_span_orders_of_magnitude"]
    
    # Also check status data for skewness (may be stored there)
    dq = status_data.get("data_quality", {})
    if "skewness" in dq and "skewness" not in metrics:
        metrics["skewness"] = dq["skewness"]
    
    return metrics
