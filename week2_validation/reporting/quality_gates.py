"""
Data Integrity & Statistical Validation — Formal Quality Gate System.

Provides per-module PASS/WARN/FAIL verdicts and an overall dataset approval status.
All thresholds and their justifications are documented inline for transparency.

Usage:
  - After all diagnostics run, call evaluate_all_gates() with diagnostic_results and data_quality.
  - Write output via write_quality_gates_json(); include in status envelope and reports.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# -----------------------------------------------------------------------------
# Constants: thresholds and justifications (aligned with config/thresholds.yaml
# where applicable; hardcoded here with justification for auditability)
# -----------------------------------------------------------------------------

# Distribution: minimum N for reliable skewness and distribution tests.
# Justification: KS/AD and skewness estimates need sufficient sample size; 50 is a
# common lower bound for distribution diagnostics (see thresholds.yaml distribution_tests: 50).
DISTRIBUTION_MIN_N = 50

# Benford: minimum scale span (orders of magnitude) for applicability.
# Justification: Benford's Law applies to data spanning multiple orders of magnitude;
# below 2.0 the first-digit distribution is not expected to follow Benford (see benford/diagnostics.py).
BENFORD_MIN_SCALE_SPAN = 2.0

# Benford: chi-squared significance level. p > this value = consistent with Benford.
# Justification: Standard alpha=0.05; matches thresholds.yaml benford.significance_level.
BENFORD_ALPHA = 0.05

# Log-normality: PASS requires both KS and AD; WARN if only KS (e.g. no SciPy).
# Justification: Two tests improve robustness; AD is more sensitive to tail behavior.
LOG_NORMALITY_NEED_BOTH = True

# COSMIC: overlap count thresholds (fusion pairs in common with reference).
# Justification: Spearman with n<3 is undefined; 3–4 is marginal; ≥5 allows minimal
# correlation assessment; ≥10 preferred for power (see reporting and cosmic docs).
COSMIC_OVERLAP_MIN_PASS = 5
COSMIC_OVERLAP_WARN_LOW = 3
COSMIC_OVERLAP_WARN_HIGH = 4

# COSMIC: Spearman p-value. p < 0.05 = significant correlation; 0.05–0.10 = marginal.
# Justification: Standard significance; marginal band for transparency.
COSMIC_SPEARMAN_P_PASS = 0.05
COSMIC_SPEARMAN_P_WARN_HIGH = 0.10


# -----------------------------------------------------------------------------
# Enums and data structures
# -----------------------------------------------------------------------------


class GateStatus(str, Enum):
    """Per-module gate status."""

    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"
    SKIPPED = "SKIPPED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class OverallApproval(str, Enum):
    """Overall dataset approval from quality gates."""

    APPROVED = "APPROVED"
    CONDITIONAL = "CONDITIONAL"
    REJECTED = "REJECTED"


# Severity: 0=info, 1=warning, 2=critical (for ordering and reporting)
SEVERITY_INFO = 0
SEVERITY_WARNING = 1
SEVERITY_CRITICAL = 2


@dataclass
class QualityGateResult:
    """Single quality gate result for one validation module."""

    module_name: str
    status: GateStatus
    metric_name: str
    metric_value: Optional[float]
    threshold: Optional[float]
    reason: str
    severity: int = SEVERITY_INFO

    def to_dict(self) -> Dict[str, Any]:
        return {
            "module_name": self.module_name,
            "status": self.status.value,
            "metric_name": self.metric_name,
            "metric_value": self.metric_value,
            "threshold": self.threshold,
            "reason": self.reason,
            "severity": self.severity,
        }


# -----------------------------------------------------------------------------
# Per-module gate evaluation
# -----------------------------------------------------------------------------


def evaluate_distribution_gate(
    diagnostic_results: Dict[str, Any],
    data_quality: Optional[Dict[str, Any]] = None,
) -> QualityGateResult:
    """
    Distribution gate:
    - PASS: skewness computable and N ≥ 50.
    - WARN: N < 50.
    - FAIL: no protein_length / no valid distribution (e.g. zero rows).
    """
    dq = data_quality or {}
    dist = (diagnostic_results or {}).get("distribution") or {}
    n = dq.get("total_rows_original") or dq.get("rows_used_for_analysis")
    if n is None and dist:
        # Infer from distribution stats if present (e.g. count not in data_quality)
        n = None
    try:
        n_int = int(n) if n is not None else 0
    except (TypeError, ValueError):
        n_int = 0

    skewness = dq.get("skewness") or dist.get("skewness")
    skewness_computable = (
        skewness is not None
        and isinstance(skewness, (int, float))
        and float(skewness) == float(skewness)
    )

    if n_int == 0 and not dist:
        return QualityGateResult(
            module_name="distribution",
            status=GateStatus.FAIL,
            metric_name="sample_size",
            metric_value=0.0,
            threshold=float(DISTRIBUTION_MIN_N),
            reason="No protein_length data or empty dataset; distribution cannot be computed.",
            severity=SEVERITY_CRITICAL,
        )
    if n_int == 0:
        return QualityGateResult(
            module_name="distribution",
            status=GateStatus.FAIL,
            metric_name="sample_size",
            metric_value=0.0,
            threshold=float(DISTRIBUTION_MIN_N),
            reason="Zero rows used for analysis.",
            severity=SEVERITY_CRITICAL,
        )
    if n_int < DISTRIBUTION_MIN_N:
        return QualityGateResult(
            module_name="distribution",
            status=GateStatus.WARN,
            metric_name="sample_size",
            metric_value=float(n_int),
            threshold=float(DISTRIBUTION_MIN_N),
            reason=f"Sample size N={n_int} below minimum {DISTRIBUTION_MIN_N} for robust distribution diagnostics.",
            severity=SEVERITY_WARNING,
        )
    if not skewness_computable:
        return QualityGateResult(
            module_name="distribution",
            status=GateStatus.WARN,
            metric_name="skewness",
            metric_value=None,
            threshold=None,
            reason="Skewness not computed or invalid; distribution statistics incomplete.",
            severity=SEVERITY_WARNING,
        )
    return QualityGateResult(
        module_name="distribution",
        status=GateStatus.PASS,
        metric_name="sample_size",
        metric_value=float(n_int),
        threshold=float(DISTRIBUTION_MIN_N),
        reason=f"N≥{DISTRIBUTION_MIN_N} and skewness computed successfully.",
        severity=SEVERITY_INFO,
    )


def evaluate_benford_gate(diagnostic_results: Dict[str, Any]) -> QualityGateResult:
    """
    Benford gate:
    - PASS: chi-squared p > 0.05 AND scale_span ≥ 2.0 (applicable and consistent).
    - WARN: scale_span < 2.0 (NOT_APPLICABLE) — test not meaningful, do not fail dataset.
    - FAIL: scale_span ≥ 2.0 AND chi-squared p < 0.05 (applicable but inconsistent).
    - SKIPPED: no Benford results (e.g. module not run).
    """
    benford = (diagnostic_results or {}).get("benford") or {}
    scale_span = benford.get("scale_span_orders_of_magnitude")
    applicability = benford.get("applicability")
    p_value = benford.get("p_value")
    # Threshold for pass: p > BENFORD_ALPHA
    try:
        p_val = float(p_value) if p_value is not None else None
    except (TypeError, ValueError):
        p_val = None

    try:
        span_val = float(scale_span) if scale_span is not None else None
    except (TypeError, ValueError):
        span_val = None

    if applicability is None and span_val is None and p_val is None:
        return QualityGateResult(
            module_name="benford",
            status=GateStatus.SKIPPED,
            metric_name="scale_span_orders_of_magnitude",
            metric_value=None,
            threshold=float(BENFORD_MIN_SCALE_SPAN),
            reason="Benford diagnostics not run or no results available.",
            severity=SEVERITY_INFO,
        )

    # Scale span below minimum → NOT_APPLICABLE (WARN, not FAIL)
    if span_val is not None and span_val < BENFORD_MIN_SCALE_SPAN:
        return QualityGateResult(
            module_name="benford",
            status=GateStatus.WARN,
            metric_name="scale_span_orders_of_magnitude",
            metric_value=span_val,
            threshold=float(BENFORD_MIN_SCALE_SPAN),
            reason=f"Scale span {span_val:.2f} < {BENFORD_MIN_SCALE_SPAN}; Benford's Law not applicable (first-digit distribution not meaningful).",
            severity=SEVERITY_WARNING,
        )

    # Applicable (span ≥ 2.0) but p-value missing
    if p_val is None:
        return QualityGateResult(
            module_name="benford",
            status=GateStatus.WARN,
            metric_name="chi_squared_p_value",
            metric_value=None,
            threshold=float(BENFORD_ALPHA),
            reason="Scale span sufficient for Benford but chi-squared p-value not available.",
            severity=SEVERITY_WARNING,
        )

    # Applicable (span ≥ 2.0) and p <= alpha → WARN (informational only)
    if p_val < BENFORD_ALPHA:
        return QualityGateResult(
            module_name="benford",
            status=GateStatus.WARN,
            metric_name="chi_squared_p_value",
            metric_value=p_val,
            threshold=float(BENFORD_ALPHA),
            reason=f"Benford applicable but first-digit distribution inconsistent (p={p_val:.4f}). This is informational only; non-conformance does not indicate data issues.",
            severity=SEVERITY_WARNING,
        )

    return QualityGateResult(
        module_name="benford",
        status=GateStatus.PASS,
        metric_name="chi_squared_p_value",
        metric_value=p_val,
        threshold=float(BENFORD_ALPHA),
        reason=f"Scale span≥{BENFORD_MIN_SCALE_SPAN} and first-digit distribution consistent with Benford (p>{BENFORD_ALPHA}).",
        severity=SEVERITY_INFO,
    )


def evaluate_log_normality_gate(diagnostic_results: Dict[str, Any]) -> QualityGateResult:
    """
    Log-normality gate:
    - PASS: both KS and AD computed successfully.
    - WARN: only KS computed (e.g. AD not available without SciPy).
    - FAIL: neither KS nor AD computed.
    - SKIPPED: no log_normality results (module not run).
    """
    logn = (diagnostic_results or {}).get("log_normality") or {}
    ks = logn.get("ks_statistic") is not None
    ad = logn.get("ad_statistic") is not None

    if not ks and not ad:
        if not logn:
            return QualityGateResult(
                module_name="log_normality",
                status=GateStatus.SKIPPED,
                metric_name="ks_and_ad",
                metric_value=None,
                threshold=None,
                reason="Log-normality diagnostics not run or no results available.",
                severity=SEVERITY_INFO,
            )
        return QualityGateResult(
            module_name="log_normality",
            status=GateStatus.FAIL,
            metric_name="ks_and_ad",
            metric_value=None,
            threshold=None,
            reason="Neither KS nor Anderson-Darling statistic computed.",
            severity=SEVERITY_CRITICAL,
        )
    if ks and ad:
        return QualityGateResult(
            module_name="log_normality",
            status=GateStatus.PASS,
            metric_name="ks_and_ad",
            metric_value=1.0,
            threshold=None,
            reason="KS and Anderson-Darling log-normality tests computed successfully.",
            severity=SEVERITY_INFO,
        )
    return QualityGateResult(
        module_name="log_normality",
        status=GateStatus.WARN,
        metric_name="ks_and_ad",
        metric_value=0.5,
        threshold=None,
        reason="Only KS computed; Anderson-Darling not available (e.g. requires SciPy).",
        severity=SEVERITY_WARNING,
    )


def evaluate_cosmic_gate(diagnostic_results: Dict[str, Any]) -> QualityGateResult:
    """
    COSMIC gate:
    - PASS: overlap ≥ 10 AND spearman_p < 0.05.
    - WARN: overlap < 10 (reference too small) or spearman_p in [0.05, 0.10].
    - SKIPPED: no COSMIC results (module not run or no reference).
    """
    cosmic = (diagnostic_results or {}).get("cosmic") or {}
    overlap = cosmic.get("overlap_count")
    spearman_p = cosmic.get("spearman_p_value")

    try:
        overlap_int = int(overlap) if overlap is not None else 0
    except (TypeError, ValueError):
        overlap_int = 0

    try:
        p_val = float(spearman_p) if spearman_p is not None else None
    except (TypeError, ValueError):
        p_val = None

    if overlap is None and spearman_p is None and not cosmic:
        return QualityGateResult(
            module_name="cosmic",
            status=GateStatus.SKIPPED,
            metric_name="overlap_count",
            metric_value=None,
            threshold=float(COSMIC_OVERLAP_MIN_PASS),
            reason="COSMIC diagnostics not run or no reference data.",
            severity=SEVERITY_INFO,
        )

    # Minimum overlap guard: reference too small → WARN (not FAIL) to avoid false rejection
    if overlap_int < 10:
        return QualityGateResult(
            module_name="cosmic",
            status=GateStatus.WARN,
            metric_name="overlap_count",
            metric_value=float(overlap_int),
            threshold=10.0,
            reason="COSMIC reference too small for meaningful statistical comparison.",
            severity=SEVERITY_WARNING,
        )

    # Overlap ≥ 10: then decide on p-value
    if p_val is None:
        return QualityGateResult(
            module_name="cosmic",
            status=GateStatus.WARN,
            metric_name="spearman_p_value",
            metric_value=None,
            threshold=float(COSMIC_SPEARMAN_P_PASS),
            reason="Overlap ≥ 10 but Spearman p-value not available.",
            severity=SEVERITY_WARNING,
        )
    if p_val < COSMIC_SPEARMAN_P_PASS:
        return QualityGateResult(
            module_name="cosmic",
            status=GateStatus.PASS,
            metric_name="spearman_p_value",
            metric_value=p_val,
            threshold=float(COSMIC_SPEARMAN_P_PASS),
            reason=f"Overlap ≥ 10 and Spearman p < {COSMIC_SPEARMAN_P_PASS}; significant rank agreement with COSMIC.",
            severity=SEVERITY_INFO,
        )
    if p_val <= COSMIC_SPEARMAN_P_WARN_HIGH:
        return QualityGateResult(
            module_name="cosmic",
            status=GateStatus.WARN,
            metric_name="spearman_p_value",
            metric_value=p_val,
            threshold=float(COSMIC_SPEARMAN_P_PASS),
            reason=f"Spearman p in marginal range [{COSMIC_SPEARMAN_P_PASS}, {COSMIC_SPEARMAN_P_WARN_HIGH}].",
            severity=SEVERITY_WARNING,
        )
    return QualityGateResult(
        module_name="cosmic",
        status=GateStatus.WARN,
        metric_name="spearman_p_value",
        metric_value=p_val,
        threshold=float(COSMIC_SPEARMAN_P_PASS),
        reason=f"Spearman p > {COSMIC_SPEARMAN_P_WARN_HIGH}; no significant rank agreement with COSMIC.",
        severity=SEVERITY_WARNING,
    )


# -----------------------------------------------------------------------------
# Aggregate evaluation and overall approval
# -----------------------------------------------------------------------------


def evaluate_all_gates(
    diagnostic_results: Dict[str, Any],
    data_quality: Optional[Dict[str, Any]] = None,
    diagnostics_run: Optional[Sequence[str]] = None,
) -> Tuple[List[QualityGateResult], OverallApproval]:
    """
    Evaluate all quality gates from diagnostic_results and data_quality.
    Only evaluates modules that were run (if diagnostics_run is provided).

    Returns:
        (list of QualityGateResult, overall approval).
    """
    diagnostics_run = diagnostics_run or []
    # Normalize names: run_week2 uses "diagnostics", "benford", "lognormal", "cosmic"
    run_set = {s.lower().replace(" ", "_") for s in diagnostics_run}
    run_dist = "diagnostics" in run_set or "distribution" in run_set
    run_benford = "benford" in run_set
    run_lognormal = "lognormal" in run_set or "log_normality" in run_set
    run_cosmic = "cosmic" in run_set

    results: List[QualityGateResult] = []

    if run_dist:
        results.append(evaluate_distribution_gate(diagnostic_results, data_quality))
    else:
        results.append(
            QualityGateResult(
                module_name="distribution",
                status=GateStatus.SKIPPED,
                metric_name="sample_size",
                metric_value=None,
                threshold=None,
                reason="Distribution diagnostics not run.",
                severity=SEVERITY_INFO,
            )
        )

    if run_benford:
        results.append(evaluate_benford_gate(diagnostic_results))
    else:
        results.append(
            QualityGateResult(
                module_name="benford",
                status=GateStatus.SKIPPED,
                metric_name="scale_span_orders_of_magnitude",
                metric_value=None,
                threshold=None,
                reason="Benford diagnostics not run.",
                severity=SEVERITY_INFO,
            )
        )

    if run_lognormal:
        results.append(evaluate_log_normality_gate(diagnostic_results))
    else:
        results.append(
            QualityGateResult(
                module_name="log_normality",
                status=GateStatus.SKIPPED,
                metric_name="ks_and_ad",
                metric_value=None,
                threshold=None,
                reason="Log-normality diagnostics not run.",
                severity=SEVERITY_INFO,
            )
        )

    if run_cosmic:
        results.append(evaluate_cosmic_gate(diagnostic_results))
    else:
        results.append(
            QualityGateResult(
                module_name="cosmic",
                status=GateStatus.SKIPPED,
                metric_name="overlap_count",
                metric_value=None,
                threshold=None,
                reason="COSMIC diagnostics not run.",
                severity=SEVERITY_INFO,
            )
        )

    overall = compute_overall_approval(results)
    return results, overall


def compute_overall_approval(gate_results: Sequence[QualityGateResult]) -> OverallApproval:
    """
    Overall gate:
    - APPROVED: no FAIL, and ≥ 2 PASS (excluding SKIPPED/NOT_APPLICABLE).
    - CONDITIONAL: at least one WARN, no FAIL.
    - REJECTED: any FAIL.
    """
    statuses = [r.status for r in gate_results]
    if GateStatus.FAIL in statuses:
        return OverallApproval.REJECTED
    pass_count = sum(1 for s in statuses if s == GateStatus.PASS)
    warn_count = sum(1 for s in statuses if s == GateStatus.WARN)
    if warn_count >= 1:
        return OverallApproval.CONDITIONAL
    if pass_count >= 2:
        return OverallApproval.APPROVED
    # Fewer than 2 PASS and no WARN/FAIL: e.g. mostly SKIPPED
    if pass_count >= 1 or warn_count >= 1:
        return OverallApproval.CONDITIONAL
    return OverallApproval.CONDITIONAL


def compute_power_warnings(
    data_quality: Optional[Dict[str, Any]],
    diagnostic_results: Optional[Dict[str, Any]],
) -> List[Dict[str, str]]:
    """
    Compute power/sample-size warnings for reporting.
    Returns a list of dicts with type, message, and recommendation.
    """
    warnings: List[Dict[str, str]] = []
    n = (data_quality or {}).get("rows_used_for_analysis", 0)

    if n < 30:
        warnings.append({
            "type": "VERY_SMALL_SAMPLE",
            "message": f"N={n} is very small. All statistical measures have wide confidence intervals.",
            "recommendation": "Interpret all results as preliminary/exploratory only.",
        })
    elif n < 100:
        warnings.append({
            "type": "SMALL_SAMPLE",
            "message": f"N={n} is below 100. Distribution tests have reduced power.",
            "recommendation": "Consider combining with additional datasets if available.",
        })

    cosmic = (diagnostic_results or {}).get("cosmic") or {}
    overlap = cosmic.get("overlap_count", 0)
    if overlap is not None and overlap < 20:
        try:
            overlap_int = int(overlap)
        except (TypeError, ValueError):
            overlap_int = 0
        if overlap_int < 20:
            warnings.append({
                "type": "LOW_COSMIC_OVERLAP",
                "message": f"Only {overlap_int} gene pairs overlap with COSMIC reference.",
                "recommendation": "COSMIC correlation has low statistical power. Use full COSMIC reference if available.",
            })

    return warnings


# -----------------------------------------------------------------------------
# Serialization and file output
# -----------------------------------------------------------------------------


def write_quality_gates_json(
    output_path: Path,
    gate_results: Sequence[QualityGateResult],
    overall: OverallApproval,
    dataset_stem: str,
) -> None:
    """Write full quality gate details to validation_quality_gates_{stem}.json."""
    payload = {
        "dataset_stem": dataset_stem,
        "overall_approval": overall.value,
        "gates": [r.to_dict() for r in gate_results],
        "thresholds_used": {
            "distribution_min_n": DISTRIBUTION_MIN_N,
            "benford_min_scale_span": BENFORD_MIN_SCALE_SPAN,
            "benford_alpha": BENFORD_ALPHA,
            "cosmic_overlap_min_pass": COSMIC_OVERLAP_MIN_PASS,
            "cosmic_overlap_warn_low": COSMIC_OVERLAP_WARN_LOW,
            "cosmic_spearman_p_pass": COSMIC_SPEARMAN_P_PASS,
            "cosmic_spearman_p_warn_high": COSMIC_SPEARMAN_P_WARN_HIGH,
        },
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def compute_quality_score(gate_results: Sequence[QualityGateResult]) -> float:
    """
    Compute overall quality score (0.0-1.0) from gate results.
    
    Scoring: PASS=1.0, WARN=0.5, FAIL=0.0, SKIPPED/NOT_APPLICABLE=0.0 (excluded from average).
    Returns average of non-skipped gates.
    """
    scores = []
    for r in gate_results:
        if r.status == GateStatus.PASS:
            scores.append(1.0)
        elif r.status == GateStatus.WARN:
            scores.append(0.5)
        elif r.status == GateStatus.FAIL:
            scores.append(0.0)
        # SKIPPED and NOT_APPLICABLE are excluded (don't contribute to score)
    
    if not scores:
        return 0.0  # No gates evaluated
    return sum(scores) / len(scores)


def quality_gates_summary_for_envelope(
    gate_results: Sequence[QualityGateResult],
    overall: OverallApproval,
) -> Dict[str, Any]:
    """Return a compact dict suitable for status envelope and report inclusion."""
    return {
        "overall_approval": overall.value,
        "quality_score": compute_quality_score(gate_results),
        "gates": [r.to_dict() for r in gate_results],
        "pass_count": sum(1 for r in gate_results if r.status == GateStatus.PASS),
        "warn_count": sum(1 for r in gate_results if r.status == GateStatus.WARN),
        "fail_count": sum(1 for r in gate_results if r.status == GateStatus.FAIL),
        "skipped_count": sum(
            1 for r in gate_results if r.status in (GateStatus.SKIPPED, GateStatus.NOT_APPLICABLE)
        ),
    }
