"""
Benford's Law Diagnostic Module — Week 2 Data Integrity Pipeline

================================================================================
DIAGNOSTIC-ONLY DISCLAIMER
================================================================================

This module performs DIAGNOSTIC-ONLY analysis of first significant digit (FSD)
distributions compared to the theoretical Benford distribution.

This module:
- Does NOT determine data correctness
- Does NOT detect anomalies or irregularities
- Does NOT assess data quality
- Does NOT enforce pass/fail gates
- Does NOT block downstream pipeline execution
- Does NOT provide biological or scientific interpretation
- Does NOT assume Benford's Law is applicable to the input data

Benford's Law applicability depends on data characteristics including:
- Sufficient scale span (multiple orders of magnitude)
- Sufficient sample size
- Data arising from multiplicative processes or spanning wide ranges

A result indicating deviation from Benford's distribution is NOT indicative of
any data problem. Many naturally occurring datasets do not follow Benford's Law.

"Benford-not-applicable" is a documented and expected outcome.

All outputs are DESCRIPTIVE STATISTICS ONLY.

================================================================================
PIPELINE CONTEXT
================================================================================

This module is part of Week 2: Data Integrity & Statistical Validation.

Week 2 enforces:
- Dataset freeze gates (handled elsewhere)
- No inference before freeze
- No hypothesis testing conclusions
- No biological interpretation
- No claims of correctness or incorrectness

================================================================================
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np

try:
    from scipy import stats as scipy_stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    scipy_stats = None

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


NumericArray = Union[List[float], List[int], "np.ndarray[Any, np.dtype[np.floating[Any]]]", Any]

BENFORD_DIGITS: Tuple[int, ...] = (1, 2, 3, 4, 5, 6, 7, 8, 9)

BENFORD_EXPECTED_PROPORTIONS: Dict[int, float] = {
    d: math.log10(1 + 1 / d) for d in BENFORD_DIGITS
}

MINIMUM_SAMPLE_SIZE_THRESHOLD: int = 50

MINIMUM_SCALE_SPAN_ORDERS_OF_MAGNITUDE: float = 2.0


@dataclass(frozen=True)
class BenfordMetadata:
    """
    Metadata container for Benford diagnostic results.

    All fields are DESCRIPTIVE and carry NO inferential meaning.
    """

    diagnostic_only: bool = True
    hypothesis_tested: bool = False
    interpretation_provided: bool = False
    applicability_assessed: bool = True
    benford_applicable: Optional[bool] = None
    reason_if_not_applicable: Optional[str] = None
    scipy_available: bool = SCIPY_AVAILABLE
    sample_size: int = 0
    scale_span_orders_of_magnitude: Optional[float] = None
    module_version: str = "1.0.0"
    pipeline_stage: str = "week2_data_integrity"
    disclaimer: str = (
        "This output is DIAGNOSTIC ONLY. No inference, conclusion, or "
        "interpretation is provided or implied. Deviation from Benford's "
        "distribution does not indicate any data issue. Benford's Law "
        "applicability is context-dependent and not guaranteed."
    )


@dataclass(frozen=True)
class BenfordDiagnosticResult:
    """
    Complete result container for Benford diagnostic analysis.

    All statistics are DESCRIPTIVE ONLY.
    No pass/fail determination is made or implied.
    """

    observed_counts: Dict[int, int]
    observed_frequencies: Dict[int, float]
    expected_frequencies: Dict[int, float]
    chi_squared_statistic: Optional[float]
    degrees_of_freedom: Optional[int]
    p_value: Optional[float]
    total_extracted_digits: int
    total_input_values: int
    values_excluded_count: int
    exclusion_reasons: Dict[str, int]
    metadata: BenfordMetadata
    computation_successful: bool
    computation_notes: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class InputPreparationResult:
    """
    Result of input data preparation.

    This is a preparatory step and does NOT determine data correctness.
    """

    is_processable: bool
    positive_values: Optional[np.ndarray]
    total_input_count: int
    positive_count: int
    excluded_count: int
    exclusion_reasons: Dict[str, int]
    notes: List[str]


def prepare_input_data(data: NumericArray) -> InputPreparationResult:
    """
    Prepare and filter input data for FSD extraction.

    This function:
    - Converts input to numpy array
    - Filters to strictly positive finite values
    - Reports exclusion counts by category
    - Does NOT make quality judgments

    Parameters
    ----------
    data : NumericArray
        Input numeric data (list, numpy array, or pandas Series).

    Returns
    -------
    InputPreparationResult
        Structured preparation result with processable values and metadata.
    """
    notes: List[str] = []
    exclusion_reasons: Dict[str, int] = {
        "non_positive": 0,
        "infinite": 0,
        "nan": 0,
        "non_numeric": 0,
    }

    if PANDAS_AVAILABLE and pd is not None:
        if isinstance(data, pd.Series):
            data = data.values
            notes.append("Converted pandas Series to numpy array.")

    try:
        arr = np.asarray(data, dtype=np.float64)
    except (ValueError, TypeError) as exc:
        return InputPreparationResult(
            is_processable=False,
            positive_values=None,
            total_input_count=0,
            positive_count=0,
            excluded_count=0,
            exclusion_reasons=exclusion_reasons,
            notes=[f"Failed to convert input to numeric array: {exc}"],
        )

    arr = arr.flatten()
    total_input_count = len(arr)

    if total_input_count == 0:
        return InputPreparationResult(
            is_processable=False,
            positive_values=None,
            total_input_count=0,
            positive_count=0,
            excluded_count=0,
            exclusion_reasons=exclusion_reasons,
            notes=["Input array is empty."],
        )

    nan_mask = np.isnan(arr)
    exclusion_reasons["nan"] = int(np.sum(nan_mask))

    inf_mask = np.isinf(arr)
    exclusion_reasons["infinite"] = int(np.sum(inf_mask))

    non_positive_mask = arr <= 0
    finite_non_positive = non_positive_mask & ~nan_mask & ~inf_mask
    exclusion_reasons["non_positive"] = int(np.sum(finite_non_positive))

    processable_mask = ~nan_mask & ~inf_mask & ~non_positive_mask
    positive_values = arr[processable_mask]

    positive_count = len(positive_values)
    excluded_count = total_input_count - positive_count

    is_processable = positive_count > 0

    if not is_processable:
        notes.append("No strictly positive finite values found in input.")

    return InputPreparationResult(
        is_processable=is_processable,
        positive_values=positive_values if is_processable else None,
        total_input_count=total_input_count,
        positive_count=positive_count,
        excluded_count=excluded_count,
        exclusion_reasons=exclusion_reasons,
        notes=notes,
    )


def extract_first_significant_digits(values: np.ndarray) -> np.ndarray:
    """
    Extract the first significant digit from each value.

    The first significant digit is the leftmost non-zero digit.
    For example: 0.00456 -> 4, 123.7 -> 1, 9.999 -> 9

    Parameters
    ----------
    values : np.ndarray
        Array of strictly positive finite values.

    Returns
    -------
    np.ndarray
        Array of first significant digits (integers 1-9).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        log_values = np.log10(values)

    floor_log = np.floor(log_values)

    mantissa = values / (10.0 ** floor_log)

    first_digits = np.floor(mantissa).astype(np.int64)

    first_digits = np.clip(first_digits, 1, 9)

    return first_digits


def compute_benford_expected_frequencies() -> Dict[int, float]:
    """
    Return the theoretical Benford distribution for digits 1-9.

    This is based on the formula: P(d) = log10(1 + 1/d)

    Returns
    -------
    Dict[int, float]
        Mapping of digit to expected proportion under Benford's Law.
    """
    return dict(BENFORD_EXPECTED_PROPORTIONS)


def compute_observed_fsd_distribution(
    first_digits: np.ndarray,
) -> Tuple[Dict[int, int], Dict[int, float]]:
    """
    Compute observed counts and frequencies of first significant digits.

    Parameters
    ----------
    first_digits : np.ndarray
        Array of first significant digits (integers 1-9).

    Returns
    -------
    Tuple[Dict[int, int], Dict[int, float]]
        Tuple of (counts dict, frequencies dict) for digits 1-9.
    """
    total = len(first_digits)

    counts: Dict[int, int] = {d: 0 for d in BENFORD_DIGITS}
    frequencies: Dict[int, float] = {d: 0.0 for d in BENFORD_DIGITS}

    if total == 0:
        return counts, frequencies

    unique, unique_counts = np.unique(first_digits, return_counts=True)

    for digit, count in zip(unique, unique_counts):
        digit_int = int(digit)
        if digit_int in counts:
            counts[digit_int] = int(count)
            frequencies[digit_int] = float(count) / float(total)

    return counts, frequencies


def assess_benford_applicability(
    sample_size: int,
    positive_values: Optional[np.ndarray],
) -> Tuple[Optional[bool], Optional[str], Optional[float]]:
    """
    Assess whether Benford's Law analysis is applicable to the data.

    This is a HEURISTIC assessment based on commonly cited guidelines.
    It does NOT determine data correctness.

    Parameters
    ----------
    sample_size : int
        Number of processable values.
    positive_values : Optional[np.ndarray]
        Array of positive values (used to compute scale span).

    Returns
    -------
    Tuple[Optional[bool], Optional[str], Optional[float]]
        Tuple of (is_applicable, reason_if_not, scale_span_orders).
    """
    if sample_size < MINIMUM_SAMPLE_SIZE_THRESHOLD:
        return (
            False,
            f"Sample size ({sample_size}) below minimum threshold "
            f"({MINIMUM_SAMPLE_SIZE_THRESHOLD}) for meaningful FSD analysis.",
            None,
        )

    if positive_values is None or len(positive_values) == 0:
        return (False, "No positive values available for scale assessment.", None)

    min_val = float(np.min(positive_values))
    max_val = float(np.max(positive_values))

    if min_val <= 0:
        return (False, "Minimum value is not strictly positive.", None)

    with np.errstate(divide="ignore", invalid="ignore"):
        scale_span = math.log10(max_val) - math.log10(min_val)

    if not math.isfinite(scale_span):
        return (False, "Unable to compute finite scale span.", None)

    if scale_span < MINIMUM_SCALE_SPAN_ORDERS_OF_MAGNITUDE:
        return (
            False,
            f"Scale span ({scale_span:.2f} orders of magnitude) below minimum "
            f"threshold ({MINIMUM_SCALE_SPAN_ORDERS_OF_MAGNITUDE}) typically "
            f"associated with Benford-distributed data.",
            scale_span,
        )

    return (True, None, scale_span)


def compute_chi_squared_diagnostic(
    observed_counts: Dict[int, int],
    expected_frequencies: Dict[int, float],
    total_count: int,
) -> Tuple[Optional[float], Optional[int], Optional[float], List[str]]:
    """
    Compute chi-squared goodness-of-fit statistic for DIAGNOSTIC purposes only.

    This is a DESCRIPTIVE statistic. The p-value is provided for reference
    but does NOT constitute hypothesis testing or inference.

    Parameters
    ----------
    observed_counts : Dict[int, int]
        Observed counts for digits 1-9.
    expected_frequencies : Dict[int, float]
        Expected Benford frequencies for digits 1-9.
    total_count : int
        Total number of observations.

    Returns
    -------
    Tuple[Optional[float], Optional[int], Optional[float], List[str]]
        Tuple of (chi_squared, degrees_of_freedom, p_value, notes).
        Returns (None, None, None, notes) if computation not possible.
    """
    notes: List[str] = []

    if total_count == 0:
        notes.append("Cannot compute chi-squared: no observations.")
        return (None, None, None, notes)

    observed_array = np.array([observed_counts[d] for d in BENFORD_DIGITS], dtype=np.float64)
    expected_array = np.array(
        [expected_frequencies[d] * total_count for d in BENFORD_DIGITS],
        dtype=np.float64,
    )

    if np.any(expected_array < 5):
        notes.append(
            "Some expected frequencies are below 5. Chi-squared approximation "
            "may be less reliable. This is a diagnostic note, not a data quality issue."
        )

    chi_squared = float(np.sum((observed_array - expected_array) ** 2 / expected_array))
    degrees_of_freedom = len(BENFORD_DIGITS) - 1

    p_value: Optional[float] = None

    if SCIPY_AVAILABLE and scipy_stats is not None:
        try:
            p_value = float(scipy_stats.chi2.sf(chi_squared, degrees_of_freedom))
        except Exception as exc:
            notes.append(f"scipy p-value computation failed: {exc}")
            p_value = None
    else:
        notes.append(
            "scipy not available. p-value not computed. "
            "Chi-squared statistic is still provided for diagnostic reference."
        )

    notes.append(
        "Chi-squared statistic and p-value are DESCRIPTIVE ONLY. "
        "No hypothesis test conclusion is drawn or implied."
    )

    return (chi_squared, degrees_of_freedom, p_value, notes)


def run_benford_diagnostics(
    data: NumericArray,
    skip_applicability_check: bool = False,
) -> BenfordDiagnosticResult:
    """
    Execute Benford's Law diagnostic analysis on input data.

    ============================================================================
    DIAGNOSTIC-ONLY DISCLAIMER
    ============================================================================

    This function performs DESCRIPTIVE statistical analysis ONLY.

    It does NOT:
    - Test hypotheses
    - Draw conclusions
    - Determine data correctness
    - Detect anomalies
    - Enforce pipeline gates
    - Provide interpretation

    Deviation from Benford's distribution is NOT indicative of any issue.
    Many datasets do not follow Benford's Law.

    ============================================================================

    Parameters
    ----------
    data : NumericArray
        Input numeric data (list, numpy array, or pandas Series).
        Only strictly positive finite values will be processed.
    skip_applicability_check : bool, optional
        If True, skip heuristic applicability assessment.
        Default is False.

    Returns
    -------
    BenfordDiagnosticResult
        Structured result containing:
        - Observed and expected FSD distributions
        - Chi-squared diagnostic statistic (if computable)
        - Applicability assessment (diagnostic only)
        - Comprehensive metadata with disclaimers
    """
    computation_notes: List[str] = []

    preparation_result = prepare_input_data(data)

    if not preparation_result.is_processable:
        metadata = BenfordMetadata(
            benford_applicable=None,
            reason_if_not_applicable="Input data could not be processed.",
            sample_size=0,
            scale_span_orders_of_magnitude=None,
        )
        return BenfordDiagnosticResult(
            observed_counts={d: 0 for d in BENFORD_DIGITS},
            observed_frequencies={d: 0.0 for d in BENFORD_DIGITS},
            expected_frequencies=compute_benford_expected_frequencies(),
            chi_squared_statistic=None,
            degrees_of_freedom=None,
            p_value=None,
            total_extracted_digits=0,
            total_input_values=preparation_result.total_input_count,
            values_excluded_count=preparation_result.excluded_count,
            exclusion_reasons=preparation_result.exclusion_reasons,
            metadata=metadata,
            computation_successful=False,
            computation_notes=preparation_result.notes,
        )

    positive_values = preparation_result.positive_values
    if positive_values is None:
        raise RuntimeError(
            "Internal error: positive_values is None despite is_processable=True. "
            "This indicates a logic inconsistency in prepare_input_data() and should not occur."
        )

    sample_size = len(positive_values)

    if skip_applicability_check:
        benford_applicable: Optional[bool] = None
        reason_if_not_applicable: Optional[str] = "Applicability check skipped by request."
        scale_span: Optional[float] = None
        computation_notes.append("Applicability assessment skipped per caller request.")
    else:
        benford_applicable, reason_if_not_applicable, scale_span = assess_benford_applicability(
            sample_size=sample_size,
            positive_values=positive_values,
        )

    first_digits = extract_first_significant_digits(positive_values)

    observed_counts, observed_frequencies = compute_observed_fsd_distribution(first_digits)

    expected_frequencies = compute_benford_expected_frequencies()

    chi_squared, dof, p_value, chi_notes = compute_chi_squared_diagnostic(
        observed_counts=observed_counts,
        expected_frequencies=expected_frequencies,
        total_count=sample_size,
    )
    computation_notes.extend(chi_notes)

    metadata = BenfordMetadata(
        benford_applicable=benford_applicable,
        reason_if_not_applicable=reason_if_not_applicable,
        sample_size=sample_size,
        scale_span_orders_of_magnitude=scale_span,
    )

    return BenfordDiagnosticResult(
        observed_counts=observed_counts,
        observed_frequencies=observed_frequencies,
        expected_frequencies=expected_frequencies,
        chi_squared_statistic=chi_squared,
        degrees_of_freedom=dof,
        p_value=p_value,
        total_extracted_digits=sample_size,
        total_input_values=preparation_result.total_input_count,
        values_excluded_count=preparation_result.excluded_count,
        exclusion_reasons=preparation_result.exclusion_reasons,
        metadata=metadata,
        computation_successful=True,
        computation_notes=computation_notes,
    )


def generate_synthetic_benford_positive_control(
    size: int = 1000,
    seed: Optional[int] = None,
) -> np.ndarray:
    """
    Generate synthetic data expected to approximately follow Benford's Law.

    This is a POSITIVE CONTROL for diagnostic pipeline testing.
    It does NOT represent real data or expected behavior.

    Parameters
    ----------
    size : int
        Number of values to generate.
    seed : Optional[int]
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Array of synthetic positive values.
    """
    rng = np.random.default_rng(seed)

    log_uniform = rng.uniform(low=0.0, high=6.0, size=size)
    values = 10.0 ** log_uniform

    return values


def generate_synthetic_benford_negative_control(
    size: int = 1000,
    seed: Optional[int] = None,
) -> np.ndarray:
    """
    Generate synthetic data NOT expected to follow Benford's Law.

    This is a NEGATIVE CONTROL for diagnostic pipeline testing.
    It does NOT represent problematic data or any quality issue.

    Parameters
    ----------
    size : int
        Number of values to generate.
    seed : Optional[int]
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Array of synthetic positive values with uniform FSD distribution.
    """
    rng = np.random.default_rng(seed)

    values = rng.uniform(low=100.0, high=999.0, size=size)

    return values


__all__ = [
    "BENFORD_DIGITS",
    "BENFORD_EXPECTED_PROPORTIONS",
    "MINIMUM_SAMPLE_SIZE_THRESHOLD",
    "MINIMUM_SCALE_SPAN_ORDERS_OF_MAGNITUDE",
    "SCIPY_AVAILABLE",
    "PANDAS_AVAILABLE",
    "BenfordMetadata",
    "BenfordDiagnosticResult",
    "InputPreparationResult",
    "prepare_input_data",
    "extract_first_significant_digits",
    "compute_benford_expected_frequencies",
    "compute_observed_fsd_distribution",
    "assess_benford_applicability",
    "compute_chi_squared_diagnostic",
    "run_benford_diagnostics",
    "generate_synthetic_benford_positive_control",
    "generate_synthetic_benford_negative_control",
]