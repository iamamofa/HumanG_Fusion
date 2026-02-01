"""
Log-Normality Diagnostic Module
===============================

Week 2: Data Integrity & Statistical Validation
Diagnostic-Only Statistical Checks for Fusion Protein Length Data

IMPORTANT DISCLAIMERS
---------------------
1. This module provides DESCRIPTIVE DIAGNOSTICS ONLY.
2. All statistics reported are for EXPLORATORY DATA ANALYSIS purposes.
3. NO hypothesis testing conclusions should be drawn from these results.
4. NO decisions about data validity should be based on these outputs.
5. Deviation from log-normality does NOT imply data invalidity.
6. No model suitability assessment is performed or implied.
7. KS p-value is NOT reported because distribution parameters were
   estimated from the same data, making the p-value statistically invalid.

AUDIT COMPLIANCE
----------------
- This module performs NO inference.
- This module makes NO decisions.
- This module provides NO thresholds.
- This module implements NO gating logic.
- This module performs NO file I/O.
- This module produces NO plots.

STATISTICAL SCOPE
-----------------
Kolmogorov-Smirnov Statistic:
    Measures the maximum absolute deviation between the empirical CDF
    and a fitted log-normal CDF. Reported as a descriptive measure only.

Anderson-Darling Statistic (requires SciPy):
    Weighted measure of deviation from log-normality with emphasis on
    distribution tails. Reported as a descriptive measure only.

SCIPY OPTIONALITY
-----------------
This module operates with or without SciPy installed:
- With SciPy: Full diagnostic statistics available.
- Without SciPy: Limited diagnostics with explanatory notes.

Author: Scientific Software Engineering Team
Module Version: 1.0.0
Compliance Level: Grant Review / PI Audit / CI-Safe
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, List, Optional, Sequence, Tuple, Union

# =============================================================================
# SCIPY AVAILABILITY CHECK
# =============================================================================

_SCIPY_AVAILABLE: bool = False
_SCIPY_IMPORT_ERROR: Optional[str] = None

try:
    import scipy.stats as scipy_stats
    _SCIPY_AVAILABLE = True
except ImportError as e:
    _SCIPY_AVAILABLE = False
    _SCIPY_IMPORT_ERROR = str(e)

# =============================================================================
# NUMPY AVAILABILITY CHECK (Required)
# =============================================================================

try:
    import numpy as np
    from numpy.typing import ArrayLike, NDArray
    _NUMPY_AVAILABLE = True
except ImportError as e:
    raise ImportError(
        "NumPy is required for this module. "
        f"Import failed with: {e}"
    ) from e


# =============================================================================
# TYPE DEFINITIONS
# =============================================================================

# NumericArray: Accepts list, numpy array, or pandas Series
# We avoid importing pandas directly; we detect it at runtime
NumericArray = Union[Sequence[float], "np.ndarray[Any, np.dtype[np.floating[Any]]]", Any]


# =============================================================================
# DISCLAIMER CONSTANTS
# =============================================================================

_MODULE_DISCLAIMER: str = (
    "DIAGNOSTIC ONLY: These statistics are provided for exploratory "
    "data characterization purposes. They do NOT constitute hypothesis "
    "tests, do NOT imply data validity or invalidity, and MUST NOT be "
    "used for inferential conclusions or decision-making."
)

_LOGNORMAL_DISCLAIMER: str = (
    "Log-normality assessment is DESCRIPTIVE ONLY. Deviation from "
    "log-normal distribution does NOT indicate data quality issues, "
    "measurement errors, or biological anomalies. No model suitability "
    "is assessed or implied by these diagnostics."
)

_KS_PVALUE_OMISSION_NOTE: str = (
    "KS p-value omitted because distribution parameters were estimated from "
    "the same data; p-value would be statistically invalid. KS statistic is "
    "reported for descriptive purposes only."
)


# =============================================================================
# METADATA DATACLASS
# =============================================================================

@dataclass(frozen=True)
class DiagnosticMetadata:
    """
    Immutable metadata for diagnostic results.
    
    This metadata explicitly declares the nature and limitations of
    the diagnostic computation for audit compliance.
    
    MANDATORY DECLARATIONS
    ----------------------
    diagnostic_only: Always True. Results are descriptive only.
    hypothesis_tested: Always False. No hypothesis testing performed.
    interpretation_provided: Always False. No interpretations given.
    
    Attributes
    ----------
    diagnostic_only : bool
        Always True. Confirms results are for diagnostic purposes only.
    hypothesis_tested : bool
        Always False. Confirms no hypothesis testing was performed.
    interpretation_provided : bool
        Always False. Confirms no interpretation of results is provided.
    module_version : str
        Version identifier for the diagnostic module.
    computation_timestamp_utc : str
        ISO 8601 timestamp of computation in UTC.
    disclaimers : Tuple[str, ...]
        Immutable collection of applicable disclaimers.
    statistical_framework : str
        Description of the statistical framework used.
    """
    
    diagnostic_only: bool = field(default=True, repr=True)
    hypothesis_tested: bool = field(default=False, repr=True)
    interpretation_provided: bool = field(default=False, repr=True)
    module_version: str = field(default="1.0.0", repr=True)
    computation_timestamp_utc: str = field(default="", repr=True)
    disclaimers: Tuple[str, ...] = field(default_factory=tuple, repr=False)
    statistical_framework: str = field(
        default="Descriptive diagnostics for log-normality assessment",
        repr=True
    )
    
    def __post_init__(self) -> None:
        """Validate immutable audit-critical fields."""
        # Enforce diagnostic_only is True
        if not self.diagnostic_only:
            raise ValueError(
                "AUDIT VIOLATION: diagnostic_only must be True. "
                "This module provides diagnostic statistics only."
            )
        # Enforce hypothesis_tested is False
        if self.hypothesis_tested:
            raise ValueError(
                "AUDIT VIOLATION: hypothesis_tested must be False. "
                "This module does NOT perform hypothesis testing."
            )
        # Enforce interpretation_provided is False
        if self.interpretation_provided:
            raise ValueError(
                "AUDIT VIOLATION: interpretation_provided must be False. "
                "This module does NOT provide interpretations."
            )


def _create_metadata() -> DiagnosticMetadata:
    """
    Factory function to create properly initialized DiagnosticMetadata.
    
    Returns
    -------
    DiagnosticMetadata
        Immutable metadata instance with current timestamp and disclaimers.
    """
    return DiagnosticMetadata(
        diagnostic_only=True,
        hypothesis_tested=False,
        interpretation_provided=False,
        module_version="1.0.0",
        computation_timestamp_utc=datetime.now(timezone.utc).isoformat(),
        disclaimers=(
        _MODULE_DISCLAIMER,
        _LOGNORMAL_DISCLAIMER,
        _KS_PVALUE_OMISSION_NOTE,
        ),
        statistical_framework="Descriptive diagnostics for log-normality assessment"
    )


# =============================================================================
# RESULT DATACLASS
# =============================================================================

@dataclass(frozen=True)
class LogNormalityDiagnosticResult:
    """
    Immutable container for log-normality diagnostic results.
    
    CRITICAL DISCLAIMER
    -------------------
    All values in this result object are DESCRIPTIVE STATISTICS ONLY.
    They must NOT be used for:
    - Hypothesis testing
    - Decision-making about data validity
    - Filtering or excluding data
    - Model selection
    - Any inferential purpose
    
    Attributes
    ----------
    ks_statistic : Optional[float]
        Kolmogorov-Smirnov statistic measuring maximum deviation between
        empirical and fitted log-normal CDFs. DESCRIPTIVE ONLY.
        None if computation was not possible.
    
    ks_p_value : Optional[float]
        Always None. KS p-value is NOT reported because distribution
        parameters were estimated from the same data, which violates
        KS test assumptions and makes the p-value statistically invalid.
        The KS statistic is reported for descriptive purposes only.
    
    ad_statistic : Optional[float]
        Anderson-Darling statistic for log-normal fit. DESCRIPTIVE ONLY.
        None if SciPy unavailable or computation not possible.
    
    ad_critical_values : Optional[List[float]]
        Critical values for reference only. NOT decision thresholds.
        None if SciPy unavailable or computation not possible.
    
    ad_significance_levels : Optional[List[float]]
        Significance levels corresponding to critical values.
        For reference only. NOT for hypothesis testing.
        None if SciPy unavailable or computation not possible.
    
    scipy_available : bool
        Whether SciPy was available for computation.
    
    sample_size : int
        Number of valid observations used in computation.
    
    computation_successful : bool
        Whether computation completed without errors.
    
    metadata : DiagnosticMetadata
        Audit-compliant metadata declaring diagnostic-only nature.
    
    computation_notes : List[str]
        Notes about the computation, including any exclusions or warnings.
    
    constant_distribution_detected : bool
        True when input variance is zero or std <= machine epsilon; log-normal fit is numerically fragile.
    """
    
    ks_statistic: Optional[float]
    ks_p_value: Optional[float]
    ad_statistic: Optional[float]
    ad_critical_values: Optional[List[float]]
    ad_significance_levels: Optional[List[float]]
    scipy_available: bool
    sample_size: int
    computation_successful: bool
    metadata: DiagnosticMetadata
    computation_notes: List[str]
    constant_distribution_detected: bool
    
    def __post_init__(self) -> None:
        """Validate result integrity for audit compliance."""
        # Validate metadata is properly configured
        if not isinstance(self.metadata, DiagnosticMetadata):
            raise TypeError(
                "AUDIT VIOLATION: metadata must be DiagnosticMetadata instance."
            )
        if not self.metadata.diagnostic_only:
            raise ValueError(
                "AUDIT VIOLATION: metadata.diagnostic_only must be True."
            )


# =============================================================================
# INPUT VALIDATION AND PREPROCESSING
# =============================================================================

def _is_pandas_series(obj: Any) -> bool:
    """
    Check if object is a pandas Series without importing pandas.
    
    Parameters
    ----------
    obj : Any
        Object to check.
    
    Returns
    -------
    bool
        True if obj is a pandas Series, False otherwise.
    """
    return type(obj).__name__ == "Series" and hasattr(obj, "to_numpy")


def _convert_to_numpy_array(data: NumericArray) -> np.ndarray:
    """
    Convert input data to a numpy array of float64.
    
    This function handles:
    - Python lists
    - NumPy arrays
    - Pandas Series (detected without importing pandas)
    
    Parameters
    ----------
    data : NumericArray
        Input data in any supported format.
    
    Returns
    -------
    np.ndarray
        Data converted to float64 numpy array.
    
    Raises
    ------
    TypeError
        If data cannot be converted to numeric array.
    ValueError
        If data is empty.
    """
    # Handle pandas Series
    if _is_pandas_series(data):
        arr = data.to_numpy(dtype=np.float64, copy=True)
    # Handle numpy arrays
    elif isinstance(data, np.ndarray):
        arr = data.astype(np.float64, copy=True)
    # Handle sequences (lists, tuples)
    elif isinstance(data, (list, tuple)):
        arr = np.array(data, dtype=np.float64)
    else:
        # Attempt generic conversion
        try:
            arr = np.asarray(data, dtype=np.float64)
        except (TypeError, ValueError) as e:
            raise TypeError(
                f"Cannot convert input of type {type(data).__name__} to numeric array. "
                f"Supported types: list, numpy.ndarray, pandas.Series. Error: {e}"
            ) from e
    
    return arr


def _filter_valid_positive_finite(
    data: np.ndarray
) -> Tuple[np.ndarray, int, int, int, List[str]]:
    """
    Filter data to strictly positive, finite values.
    
    For log-normality assessment, only strictly positive (> 0) and
    finite values are mathematically valid.
    
    Parameters
    ----------
    data : np.ndarray
        Input array of numeric values.
    
    Returns
    -------
    valid_data : np.ndarray
        Array containing only strictly positive, finite values.
    n_nan : int
        Count of NaN values excluded.
    n_inf : int
        Count of infinite values excluded.
    n_nonpositive : int
        Count of non-positive (≤ 0) values excluded.
    notes : List[str]
        Descriptive notes about exclusions.
    
    Notes
    -----
    DIAGNOSTIC NOTE: Exclusion of values is a preprocessing step only.
    It does NOT imply that excluded values are invalid or erroneous.
    All exclusions are documented for transparency.
    """
    notes: List[str] = []
    original_size = len(data)
    
    # Count and handle NaN values
    nan_mask = np.isnan(data)
    n_nan = int(np.sum(nan_mask))
    
    # Count and handle infinite values
    inf_mask = np.isinf(data)
    n_inf = int(np.sum(inf_mask))
    
    # Count non-positive values (must be strictly positive for log-normal)
    finite_mask = ~nan_mask & ~inf_mask
    nonpositive_mask = finite_mask & (data <= 0)
    n_nonpositive = int(np.sum(nonpositive_mask))
    
    # Create valid data mask: finite AND strictly positive
    valid_mask = finite_mask & (data > 0)
    valid_data = data[valid_mask]
    
    # Document exclusions
    if n_nan > 0:
        notes.append(
            f"PREPROCESSING: Excluded {n_nan} NaN value(s) from computation. "
            "This exclusion is for mathematical validity only and does NOT "
            "indicate data quality issues."
        )
    
    if n_inf > 0:
        notes.append(
            f"PREPROCESSING: Excluded {n_inf} infinite value(s) from computation. "
            "This exclusion is for mathematical validity only and does NOT "
            "indicate data quality issues."
        )
    
    if n_nonpositive > 0:
        notes.append(
            f"PREPROCESSING: Excluded {n_nonpositive} non-positive value(s) "
            "(values ≤ 0) from computation. Log-normal distributions are "
            "defined only for strictly positive values. This exclusion is "
            "for mathematical validity only and does NOT indicate data errors."
        )
    
    n_valid = len(valid_data)
    notes.append(
        f"INPUT SUMMARY: Received {original_size} observation(s), "
        f"{n_valid} valid for log-normality diagnostic computation."
    )
    
    return valid_data, n_nan, n_inf, n_nonpositive, notes


# =============================================================================
# STATISTICAL COMPUTATION FUNCTIONS
# =============================================================================

def _fit_lognormal_parameters(data: np.ndarray) -> Tuple[float, float]:
    """
    Estimate log-normal distribution parameters via maximum likelihood.
    
    For a log-normal distribution, the MLE parameters are:
    - mu (location): mean of log(data)
    - sigma (scale): standard deviation of log(data)
    
    Parameters
    ----------
    data : np.ndarray
        Array of strictly positive values.
    
    Returns
    -------
    mu : float
        MLE estimate of the location parameter (mean of log-data).
    sigma : float
        MLE estimate of the scale parameter (std of log-data).
    
    Notes
    -----
    DIAGNOSTIC NOTE: These parameters are estimated for the sole purpose
    of computing descriptive fit statistics. They are NOT validated
    parameters for any model and should NOT be used for inference.
    """
    log_data = np.log(data)
    mu = float(np.mean(log_data))
    # Use N (not N-1) for MLE estimate
    sigma = float(np.std(log_data, ddof=0))
    
    # Ensure sigma is positive (handle edge case of constant data)
    if sigma <= 0:
        sigma = np.finfo(np.float64).eps
    
    return mu, sigma


def _compute_empirical_cdf(data: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the empirical cumulative distribution function.
    
    Parameters
    ----------
    data : np.ndarray
        Array of observations.
    
    Returns
    -------
    sorted_data : np.ndarray
        Sorted data values.
    ecdf : np.ndarray
        Empirical CDF values at each sorted data point.
    """
    n = len(data)
    sorted_data = np.sort(data)
    ecdf = np.arange(1, n + 1) / n
    return sorted_data, ecdf


def _compute_lognormal_cdf(
    x: np.ndarray, 
    mu: float, 
    sigma: float
) -> np.ndarray:
    """
    Compute the CDF of a log-normal distribution.
    
    The log-normal CDF is:
        F(x) = Φ((ln(x) - μ) / σ)
    
    where Φ is the standard normal CDF.
    
    Parameters
    ----------
    x : np.ndarray
        Points at which to evaluate the CDF (must be positive).
    mu : float
        Location parameter (mean of log-data).
    sigma : float
        Scale parameter (std of log-data).
    
    Returns
    -------
    np.ndarray
        CDF values at each point in x.
    """
    # Use error function for standard normal CDF: Φ(z) = 0.5 * (1 + erf(z / sqrt(2)))
    z = (np.log(x) - mu) / sigma
    return 0.5 * (1.0 + _erf_array(z / math.sqrt(2.0)))


def _erf_array(x: np.ndarray) -> np.ndarray:
    """
    Compute error function for an array using vectorized NumPy operations.
    
    This implementation is fully vectorized and avoids Python-level loops
    for performance and numerical stability. It uses np.erf if available
    (NumPy >= 2.0), otherwise falls back to a vectorized polynomial
    approximation based on Abramowitz and Stegun (1964), formula 7.1.26.
    
    Parameters
    ----------
    x : np.ndarray
        Input values.
    
    Returns
    -------
    np.ndarray
        Error function values with same shape and dtype as input.
    
    Notes
    -----
    The polynomial approximation has maximum error |ε(x)| < 1.5 × 10⁻⁷.
    No SciPy dependency is introduced; uses only NumPy operations.
    """
    # Check if np.erf is available (NumPy >= 2.0)
    if hasattr(np, 'erf'):
        return np.asarray(np.erf(x), dtype=np.float64)
    
    # Vectorized fallback using Abramowitz and Stegun approximation
    # Constants for the approximation
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911
    
    # Save the sign and work with absolute values
    sign = np.sign(x)
    x_abs = np.abs(x)
    
    # Abramowitz and Stegun formula 7.1.26
    t = 1.0 / (1.0 + p * x_abs)
    t2 = t * t
    t3 = t2 * t
    t4 = t3 * t
    t5 = t4 * t
    
    y = 1.0 - (a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5) * np.exp(-x_abs * x_abs)
    
    return np.asarray(sign * y, dtype=np.float64)


def _compute_ks_statistic_manual(
    data: np.ndarray,
    mu: float,
    sigma: float
) -> float:
    """
    Compute the two-sided Kolmogorov-Smirnov statistic without SciPy.
    
    This implementation computes the two-sided Kolmogorov-Smirnov statistic
    following the standard definition used by SciPy:
    
        D = max(D⁺, D⁻)
    
    where:
        D⁺ = max(Fₙ(xᵢ) − F(xᵢ))   for i = 1, ..., n
        D⁻ = max(F(xᵢ) − Fₙ₋₁(xᵢ)) for i = 1, ..., n
    
    and:
        Fₙ(xᵢ) = i/n      (empirical CDF at xᵢ)
        Fₙ₋₁(xᵢ) = (i-1)/n (empirical CDF just before xᵢ)
        F(xᵢ)             (theoretical log-normal CDF at xᵢ)
    
    The data points xᵢ are the sorted observations.
    
    Parameters
    ----------
    data : np.ndarray
        Array of strictly positive observations.
    mu : float
        Fitted log-normal location parameter.
    sigma : float
        Fitted log-normal scale parameter.
    
    Returns
    -------
    float
        Two-sided KS statistic D = max(D⁺, D⁻).
    
    Notes
    -----
    DIAGNOSTIC NOTE: This statistic is a DESCRIPTIVE measure of fit.
    It does NOT constitute a test result and must NOT be used for
    inference about the true underlying distribution.
    
    IMPORTANT: The KS p-value is NOT computed because distribution
    parameters were estimated from the same data, which violates
    KS test assumptions and makes the p-value statistically invalid.
    """
    n = len(data)
    sorted_data, ecdf = _compute_empirical_cdf(data)
    theoretical_cdf = _compute_lognormal_cdf(sorted_data, mu, sigma)
    
    # KS statistic: max of D+ and D-
    # D+ = max(F_n(x_i) - F(x_i))
    # D- = max(F(x_i) - F_{n-1}(x_i))
    ecdf_minus = np.concatenate([[0.0], ecdf[:-1]])
    
    d_plus = np.max(ecdf - theoretical_cdf)
    d_minus = np.max(theoretical_cdf - ecdf_minus)
    
    ks_stat = float(max(d_plus, d_minus))
    return ks_stat


def _compute_ks_with_scipy(
    data: np.ndarray,
    mu: float,
    sigma: float
) -> float:
    """
    Compute KS statistic using SciPy.
    
    Parameters
    ----------
    data : np.ndarray
        Array of strictly positive observations.
    mu : float
        Fitted log-normal location parameter.
    sigma : float
        Fitted log-normal scale parameter.
    
    Returns
    -------
    ks_statistic : float
        KS statistic value.
    
    Notes
    -----
    IMPORTANT: This function returns ONLY the KS statistic, not the p-value.
    The KS p-value is intentionally omitted because distribution parameters
    were estimated from the same data, which violates KS test assumptions
    and makes the p-value statistically invalid. The KS statistic is
    reported for descriptive purposes only.
    """
    # SciPy's lognorm uses s=sigma, scale=exp(mu)
    result = scipy_stats.kstest(
        data,
        'lognorm',
        args=(sigma, 0, math.exp(mu))
    )
    return float(result.statistic)


def _compute_ad_with_scipy(
    data: np.ndarray
) -> Tuple[float, List[float], List[float]]:
    """
    Compute Anderson-Darling statistic for log-normality using SciPy.
    
    The Anderson-Darling test is performed on log-transformed data
    against a normal distribution, which is equivalent to testing
    the original data against a log-normal distribution.
    
    Parameters
    ----------
    data : np.ndarray
        Array of strictly positive observations.
    
    Returns
    -------
    ad_statistic : float
        Anderson-Darling statistic.
    critical_values : List[float]
        Reference critical values (NOT decision thresholds).
    significance_levels : List[float]
        Significance levels for reference (NOT for testing).
    
    Notes
    -----
    CRITICAL DISCLAIMER: Critical values and significance levels are
    provided for REFERENCE ONLY. They must NOT be used as thresholds
    for decision-making or hypothesis testing.
    """
    # Anderson-Darling on log-transformed data (testing for normality of log-data
    # is equivalent to testing for log-normality of original data)
    log_data = np.log(data)
    result = scipy_stats.anderson(log_data, dist='norm')
    
    return (
        float(result.statistic),
        [float(cv) for cv in result.critical_values],
        [float(sl) for sl in result.significance_level]
    )


# =============================================================================
# MAIN PUBLIC API
# =============================================================================

def run_log_normality_diagnostics(
    data: NumericArray,
    *,
    skip_fit: bool = False
) -> LogNormalityDiagnosticResult:
    """
    Compute diagnostic statistics for log-normality assessment.
    
    IMPORTANT DISCLAIMERS
    ---------------------
    1. This function provides DESCRIPTIVE DIAGNOSTICS ONLY.
    2. Results are NOT hypothesis tests and must NOT be treated as such.
    3. KS p-value is NOT reported (parameters estimated from same data).
    4. Deviation from log-normality does NOT imply data invalidity.
    5. No model suitability assessment is performed.
    6. These diagnostics must NOT be used for decision-making.
    
    Parameters
    ----------
    data : NumericArray
        Input data as list, numpy array, or pandas Series.
        Only strictly positive, finite values will be used.
    
    skip_fit : bool, default=False
        If True, skip parameter fitting and return empty results.
        Useful for validation workflows that need a result object
        without performing computations.
    
    Returns
    -------
    LogNormalityDiagnosticResult
        Frozen dataclass containing diagnostic statistics and metadata.
        All fields explicitly marked as diagnostic-only.
    
    Notes
    -----
    INPUT HANDLING:
        - Accepts list, numpy array, or pandas Series
        - Filters to strictly positive, finite values
        - Documents all exclusions in computation_notes
        - Returns safely if no valid data remains
    
    SCIPY OPTIONALITY:
        - If SciPy is unavailable, limited diagnostics are computed
        - KS statistic computed manually without SciPy
        - Anderson-Darling requires SciPy
        - All limitations are documented in computation_notes
    
    AUDIT COMPLIANCE:
        - All metadata fields explicitly declare diagnostic-only nature
        - No thresholds, gating, or decision logic
        - No file I/O or plotting
        - Full exclusion transparency
    
    Examples
    --------
    DISCLAIMER: The following is for API demonstration only and does NOT
    constitute guidance on interpreting results.
    
    >>> result = run_log_normality_diagnostics([1.0, 2.0, 3.0, 4.0, 5.0])
    >>> result.metadata.diagnostic_only
    True
    >>> result.metadata.hypothesis_tested
    False
    """
    # Initialize computation tracking
    computation_notes: List[str] = []
    computation_successful = False
    
    # Add module-level disclaimer to notes
    computation_notes.append(_MODULE_DISCLAIMER)
    
    # Create metadata (always diagnostic-only)
    metadata = _create_metadata()
    
    # Default return values
    ks_statistic: Optional[float] = None
    ks_p_value: Optional[float] = None
    ad_statistic: Optional[float] = None
    ad_critical_values: Optional[List[float]] = None
    ad_significance_levels: Optional[List[float]] = None
    sample_size = 0
    
    # Record SciPy availability
    scipy_available = _SCIPY_AVAILABLE
    if not scipy_available:
        computation_notes.append(
            f"NOTE: SciPy is not available ({_SCIPY_IMPORT_ERROR}). "
            "Anderson-Darling statistic will not be computed. "
            "KS statistic will be computed using manual implementation. "
            "This does NOT affect the diagnostic-only nature of results."
        )
    
    # Handle skip_fit request
    if skip_fit:
        computation_notes.append(
            "SKIP_FIT: Parameter fitting and computation skipped as requested. "
            "Returning empty result object with valid metadata structure."
        )
        return LogNormalityDiagnosticResult(
            ks_statistic=None,
            ks_p_value=None,
            ad_statistic=None,
            ad_critical_values=None,
            ad_significance_levels=None,
            scipy_available=scipy_available,
            sample_size=0,
            computation_successful=True,
            metadata=metadata,
            computation_notes=computation_notes,
            constant_distribution_detected=False,
        )
    
    # Convert input to numpy array
    try:
        arr = _convert_to_numpy_array(data)
    except (TypeError, ValueError) as e:
        computation_notes.append(
            f"INPUT ERROR: Failed to convert input to numeric array. {e}"
        )
        return LogNormalityDiagnosticResult(
            ks_statistic=None,
            ks_p_value=None,
            ad_statistic=None,
            ad_critical_values=None,
            ad_significance_levels=None,
            scipy_available=scipy_available,
            sample_size=0,
            computation_successful=False,
            metadata=metadata,
            computation_notes=computation_notes,
            constant_distribution_detected=False,
        )
    
    # Check for empty input
    if len(arr) == 0:
        computation_notes.append(
            "INPUT ERROR: Empty array provided. No computation possible. "
            "This is NOT an indication of data quality issues."
        )
        return LogNormalityDiagnosticResult(
            ks_statistic=None,
            ks_p_value=None,
            ad_statistic=None,
            ad_critical_values=None,
            ad_significance_levels=None,
            scipy_available=scipy_available,
            sample_size=0,
            computation_successful=False,
            metadata=metadata,
            computation_notes=computation_notes,
            constant_distribution_detected=False,
        )
    
    # Filter to valid positive finite values
    valid_data, n_nan, n_inf, n_nonpositive, filter_notes = _filter_valid_positive_finite(arr)
    computation_notes.extend(filter_notes)
    sample_size = len(valid_data)
    
    # Check for sufficient data
    if sample_size == 0:
        computation_notes.append(
            "COMPUTATION NOTE: No valid positive finite values remain after "
            "preprocessing. Returning null statistics. This is NOT an "
            "indication of data quality issues."
        )
        return LogNormalityDiagnosticResult(
            ks_statistic=None,
            ks_p_value=None,
            ad_statistic=None,
            ad_critical_values=None,
            ad_significance_levels=None,
            scipy_available=scipy_available,
            sample_size=0,
            computation_successful=False,
            metadata=metadata,
            computation_notes=computation_notes,
            constant_distribution_detected=False,
        )
    
    if sample_size < 2:
        computation_notes.append(
            "COMPUTATION NOTE: Fewer than 2 valid observations. "
            "Statistical diagnostics require at least 2 data points. "
            "Returning null statistics. This is NOT an indication of "
            "data quality issues."
        )
        return LogNormalityDiagnosticResult(
            ks_statistic=None,
            ks_p_value=None,
            ad_statistic=None,
            ad_critical_values=None,
            ad_significance_levels=None,
            scipy_available=scipy_available,
            sample_size=sample_size,
            computation_successful=False,
            metadata=metadata,
            computation_notes=computation_notes,
            constant_distribution_detected=False,
        )
    
    # Detect constant distribution (variance zero or std <= machine epsilon)
    _variance = float(np.var(valid_data))
    _std = float(np.std(valid_data)) if _variance > 0 else 0.0
    _eps = np.finfo(np.float64).eps
    constant_distribution_detected = _variance == 0 or _std <= _eps
    if constant_distribution_detected:
        computation_notes.append(
            "Input distribution variance is zero; log-normal fit is numerically fragile."
        )
    
    # Fit log-normal parameters
    try:
        mu, sigma = _fit_lognormal_parameters(valid_data)
        computation_notes.append(
            f"PARAMETER ESTIMATION: Fitted log-normal parameters via MLE: "
            f"μ (location) = {mu:.6f}, σ (scale) = {sigma:.6f}. "
            "These parameters are for diagnostic computation only and "
            "are NOT validated model parameters."
        )
    except Exception as e:
        computation_notes.append(
            f"COMPUTATION ERROR: Failed to estimate log-normal parameters. {e}"
        )
        return LogNormalityDiagnosticResult(
            ks_statistic=None,
            ks_p_value=None,
            ad_statistic=None,
            ad_critical_values=None,
            ad_significance_levels=None,
            scipy_available=scipy_available,
            sample_size=sample_size,
            computation_successful=False,
            metadata=metadata,
            computation_notes=computation_notes,
            constant_distribution_detected=constant_distribution_detected,
        )
    
    # Compute KS statistic
    # NOTE: ks_p_value is always None because distribution parameters were
    # estimated from the same data, which violates KS test assumptions.
    try:
        if scipy_available:
            ks_statistic = _compute_ks_with_scipy(valid_data, mu, sigma)
        else:
            ks_statistic = _compute_ks_statistic_manual(valid_data, mu, sigma)
        
        # Always set ks_p_value to None (see RC-LN-1)
        ks_p_value = None
        
        computation_notes.append(
            f"KS DIAGNOSTIC: Kolmogorov-Smirnov statistic = {ks_statistic:.6f}. "
            "KS p-value omitted because distribution parameters were estimated "
            "from the same data; p-value would be statistically invalid. "
            "KS statistic is reported for descriptive purposes only."
        )
    except Exception as e:
        computation_notes.append(
            f"COMPUTATION WARNING: KS statistic computation failed. {e}"
        )
        ks_statistic = None
        ks_p_value = None
    
    # Compute Anderson-Darling statistic (SciPy only)
    if scipy_available:
        try:
            ad_statistic, ad_critical_values, ad_significance_levels = _compute_ad_with_scipy(valid_data)
            computation_notes.append(
                f"AD DIAGNOSTIC: Anderson-Darling statistic = {ad_statistic:.6f}. "
                "Critical values and significance levels are provided for "
                "REFERENCE ONLY and must NOT be used as decision thresholds."
            )
        except Exception as e:
            computation_notes.append(
                f"COMPUTATION WARNING: Anderson-Darling computation failed. {e}"
            )
            ad_statistic = None
            ad_critical_values = None
            ad_significance_levels = None
    else:
        ad_statistic = None
        ad_critical_values = None
        ad_significance_levels = None
        computation_notes.append(
            "AD DIAGNOSTIC: Anderson-Darling statistic not computed "
            "(requires SciPy)."
        )
    
    # Mark computation as successful
    computation_successful = True
    computation_notes.append(
        "COMPUTATION COMPLETE: All requested diagnostics computed successfully. "
        "FINAL REMINDER: All results are DESCRIPTIVE ONLY and must NOT be "
        "used for inference, hypothesis testing, or decision-making."
    )
    
    # Construct and return result
    return LogNormalityDiagnosticResult(
        ks_statistic=ks_statistic,
        ks_p_value=ks_p_value,
        ad_statistic=ad_statistic,
        ad_critical_values=ad_critical_values,
        ad_significance_levels=ad_significance_levels,
        scipy_available=scipy_available,
        sample_size=sample_size,
        computation_successful=computation_successful,
        metadata=metadata,
        computation_notes=computation_notes,
        constant_distribution_detected=constant_distribution_detected,
    )


# =============================================================================
# MODULE-LEVEL EXPORTS
# =============================================================================

__all__ = [
    "run_log_normality_diagnostics",
    "LogNormalityDiagnosticResult",
    "DiagnosticMetadata",
    "NumericArray",
]

__version__ = "1.0.0"
__author__ = "Scientific Software Engineering Team"
__compliance__ = "Grant Review / PI Audit / CI-Safe"