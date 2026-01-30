"""
week2_validation/distributions/visualize.py

DIAGNOSTIC-ONLY module for structural characterization of fusion protein length data.

====================================================================================
                            CRITICAL DISCLAIMERS
====================================================================================

This module is PURELY DIAGNOSTIC. It provides ONLY:
    - Descriptive statistics (no inference)
    - Visual inspection tools (no interpretation)
    - Structural characterization (no hypothesis testing)

This module does NOT:
    - Perform hypothesis testing of any kind
    - Fit statistical models or distributions
    - Conduct power-law, Zipf's Law, or Benford's Law analysis
    - Provide biological interpretation
    - Make assumptions about data correctness or significance
    - Save files to disk
    - Depend on external pipeline state or configuration

PURPOSE:
    Week 2 of the cancer genomics pipeline focuses on data integrity and
    structural validation. This module enables researchers to visually and
    numerically inspect the distributional characteristics of protein length
    data BEFORE any downstream inference or modeling.

USAGE CONTEXT:
    - Safe for local Python/Conda environments
    - Safe for CI pipelines
    - Safe for Docker containers
    - Degrades gracefully if optional dependencies are missing

====================================================================================
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple, Union

# =============================================================================
# REQUIRED DEPENDENCY: NumPy
# =============================================================================
import numpy as np

# =============================================================================
# REQUIRED DEPENDENCY: Matplotlib
# =============================================================================
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# =============================================================================
# OPTIONAL DEPENDENCY: SciPy (for bias-corrected skewness)
# =============================================================================
_SCIPY_AVAILABLE: bool = False
try:
    import scipy.stats

    _SCIPY_AVAILABLE = True
except ImportError:
    scipy = None  # type: ignore[assignment]

# =============================================================================
# OPTIONAL DEPENDENCY: pandas (for Series input support)
# =============================================================================
_PANDAS_AVAILABLE: bool = False
try:
    import pandas as pd

    _PANDAS_AVAILABLE = True
except ImportError:
    pd = None  # type: ignore[assignment]


# =============================================================================
# TYPE ALIASES
# =============================================================================
ArrayLike = Union[List[float], np.ndarray, "pd.Series"]  # type: ignore[name-defined]
LogMode = Optional[Literal["transform", "axis"]]


# =============================================================================
# DIAGNOSTIC METADATA CONTAINER
# =============================================================================
@dataclass(frozen=True)
class DiagnosticMetadata:
    """
    Immutable metadata affirming the diagnostic-only nature of results.

    This metadata is attached to all outputs to explicitly declare that:
        - No hypothesis testing was performed
        - No statistical inference was conducted
        - No biological interpretation is provided

    DIAGNOSTIC ONLY. NO INFERENCE. NO INTERPRETATION.
    """

    diagnostic_only: bool = True
    hypothesis_tested: bool = False
    interpretation_provided: bool = False
    module_version: str = "2.0.0"
    scipy_available: bool = field(default_factory=lambda: _SCIPY_AVAILABLE)


# =============================================================================
# RESULT CONTAINERS
# =============================================================================
@dataclass(frozen=True)
class DescriptiveStatistics:
    """
    Container for purely descriptive statistics.

    DIAGNOSTIC ONLY. These statistics describe structural properties of the
    data distribution. They do NOT imply correctness, significance, or
    biological meaning.

    Attributes
    ----------
    count : int
        Number of valid observations (N).
    minimum : float
        Smallest observed value.
    maximum : float
        Largest observed value.
    mean : float
        Arithmetic mean.
    median : float
        Median (50th percentile).
    std : float
        Sample standard deviation (ddof=1).
    skewness : float
        Descriptive skewness. NaN if not computable.
    skewness_method : str
        Method used to compute skewness ("scipy_bias_corrected", "numpy_fallback", "not_computed").

    Notes
    -----
    DIAGNOSTIC ONLY. NO INFERENCE. NO INTERPRETATION.
    """

    count: int
    minimum: float
    maximum: float
    mean: float
    median: float
    std: float
    skewness: float
    skewness_method: str
    metadata: DiagnosticMetadata = field(default_factory=DiagnosticMetadata)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary containing all statistics and metadata.

        Notes
        -----
        DIAGNOSTIC ONLY. NO INFERENCE. NO INTERPRETATION.
        """
        return {
            "count": self.count,
            "minimum": self.minimum,
            "maximum": self.maximum,
            "mean": self.mean,
            "median": self.median,
            "std": self.std,
            "skewness": self.skewness,
            "skewness_method": self.skewness_method,
            "metadata": {
                "diagnostic_only": self.metadata.diagnostic_only,
                "hypothesis_tested": self.metadata.hypothesis_tested,
                "interpretation_provided": self.metadata.interpretation_provided,
                "module_version": self.metadata.module_version,
                "scipy_available": self.metadata.scipy_available,
            },
        }


@dataclass
class DiagnosticResult:
    """
    Complete result container for diagnostic analysis.

    DIAGNOSTIC ONLY. NO INFERENCE. NO INTERPRETATION.

    Attributes
    ----------
    statistics : DescriptiveStatistics
        Computed descriptive statistics.
    figure_linear : Figure
        Matplotlib Figure with linear-scale histogram.
    figure_log : Optional[Figure]
        Matplotlib Figure with log-scale histogram, or None if skipped.
    log_mode_used : Optional[str]
        The log mode that was applied ("transform", "axis", or None).
    metadata : DiagnosticMetadata
        Explicit metadata affirming diagnostic-only status.
    """

    statistics: DescriptiveStatistics
    figure_linear: Figure
    figure_log: Optional[Figure]
    log_mode_used: Optional[str]
    metadata: DiagnosticMetadata = field(default_factory=DiagnosticMetadata)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary representation (excluding Figure objects).

        Returns
        -------
        Dict[str, Any]
            Dictionary containing statistics, log mode, and metadata.
            Figures are not serializable and are excluded.

        Notes
        -----
        DIAGNOSTIC ONLY. NO INFERENCE. NO INTERPRETATION.
        """
        return {
            "statistics": self.statistics.to_dict(),
            "log_mode_used": self.log_mode_used,
            "figures_generated": {
                "linear": self.figure_linear is not None,
                "log": self.figure_log is not None,
            },
            "metadata": {
                "diagnostic_only": self.metadata.diagnostic_only,
                "hypothesis_tested": self.metadata.hypothesis_tested,
                "interpretation_provided": self.metadata.interpretation_provided,
                "module_version": self.metadata.module_version,
                "scipy_available": self.metadata.scipy_available,
            },
        }


# =============================================================================
# INPUT VALIDATION
# =============================================================================
class ValidationError(ValueError):
    """
    Raised when input data fails validation requirements.

    This is a fail-fast mechanism to ensure data integrity before
    any diagnostic computation.
    """

    pass


def validate_input(data: ArrayLike) -> np.ndarray:
    """
    Validate and convert input data to a clean NumPy array.

    This function enforces strict requirements on input data:
        - Must be array-like (list, numpy array, or pandas Series)
        - Must contain numeric values
        - Must contain only finite values (no NaN, no Inf)
        - Must contain only strictly positive values (> 0)
        - Must not be empty

    Parameters
    ----------
    data : ArrayLike
        Input data as Python list, NumPy array, or pandas Series.

    Returns
    -------
    np.ndarray
        Validated, 1-dimensional NumPy array of float64 values.

    Raises
    ------
    ValidationError
        If any validation requirement is violated. Error messages are
        explicit and actionable.
    TypeError
        If input type is not supported.

    Notes
    -----
    DIAGNOSTIC ONLY. This validation ensures data integrity for
    descriptive analysis. It does NOT validate biological correctness
    or scientific significance.

    Examples
    --------
    >>> arr = validate_input([100, 200, 300])
    >>> arr.dtype
    dtype('float64')
    """
    # -------------------------------------------------------------------------
    # Step 1: Type checking and conversion
    # -------------------------------------------------------------------------
    if isinstance(data, np.ndarray):
        arr = data.copy()
    elif isinstance(data, list):
        try:
            arr = np.array(data)
        except (ValueError, TypeError) as exc:
            raise ValidationError(
                f"Failed to convert list to array. Ensure all elements are numeric. "
                f"Original error: {exc}"
            ) from exc
    elif _PANDAS_AVAILABLE and isinstance(data, pd.Series):
        arr = data.to_numpy().copy()
    else:
        supported = "list, numpy.ndarray"
        if _PANDAS_AVAILABLE:
            supported += ", pandas.Series"
        raise TypeError(
            f"Unsupported input type: {type(data).__name__}. "
            f"Supported types: {supported}."
        )

    # -------------------------------------------------------------------------
    # Step 2: Ensure 1-dimensional
    # -------------------------------------------------------------------------
    if arr.ndim == 0:
        raise ValidationError(
            "Input is a 0-dimensional scalar. Expected a 1-dimensional array."
        )
    if arr.ndim > 1:
        raise ValidationError(
            f"Input has {arr.ndim} dimensions. Expected a 1-dimensional array. "
            f"Shape received: {arr.shape}."
        )

    # -------------------------------------------------------------------------
    # Step 3: Ensure non-empty
    # -------------------------------------------------------------------------
    if arr.size == 0:
        raise ValidationError("Input array is empty. At least one value is required.")

    # -------------------------------------------------------------------------
    # Step 4: Convert to float64 for consistent numeric handling
    # -------------------------------------------------------------------------
    try:
        arr = arr.astype(np.float64)
    except (ValueError, TypeError) as exc:
        raise ValidationError(
            f"Failed to convert input to numeric (float64). "
            f"Ensure all values are numeric. Original error: {exc}"
        ) from exc

    # -------------------------------------------------------------------------
    # Step 5: Check for NaN values
    # -------------------------------------------------------------------------
    nan_count = np.isnan(arr).sum()
    if nan_count > 0:
        raise ValidationError(
            f"Input contains {nan_count} NaN value(s). "
            f"All values must be finite. Remove or impute NaN values before analysis."
        )

    # -------------------------------------------------------------------------
    # Step 6: Check for infinite values
    # -------------------------------------------------------------------------
    inf_count = np.isinf(arr).sum()
    if inf_count > 0:
        raise ValidationError(
            f"Input contains {inf_count} infinite value(s). "
            f"All values must be finite."
        )

    # -------------------------------------------------------------------------
    # Step 7: Check for strictly positive values
    # -------------------------------------------------------------------------
    non_positive_count = (arr <= 0).sum()
    if non_positive_count > 0:
        min_val = arr.min()
        raise ValidationError(
            f"Input contains {non_positive_count} non-positive value(s). "
            f"All values must be strictly positive (> 0). "
            f"Minimum value found: {min_val}."
        )

    return arr


# =============================================================================
# SKEWNESS COMPUTATION (WITH FALLBACK)
# =============================================================================
def _compute_skewness_scipy(arr: np.ndarray) -> Tuple[float, str]:
    """
    Compute bias-corrected skewness using scipy.stats.skew.

    Parameters
    ----------
    arr : np.ndarray
        Validated input array.

    Returns
    -------
    Tuple[float, str]
        (skewness_value, method_name)

    Notes
    -----
    DIAGNOSTIC ONLY. Skewness is a descriptive measure of asymmetry.
    NO INFERENCE is made about the underlying distribution.
    """
    # scipy.stats.skew with bias=False applies the bias correction
    # (adjusts for sample size)
    skew_val = scipy.stats.skew(arr, bias=False)
    return float(skew_val), "scipy_bias_corrected"


def _compute_skewness_numpy(arr: np.ndarray) -> Tuple[float, str]:
    """
    Compute skewness using NumPy-only implementation (fallback).

    Uses the adjusted Fisher-Pearson standardized moment coefficient
    (bias-corrected formula) to match scipy's bias=False behavior.

    Parameters
    ----------
    arr : np.ndarray
        Validated input array.

    Returns
    -------
    Tuple[float, str]
        (skewness_value, method_name)

    Notes
    -----
    DIAGNOSTIC ONLY. This is a fallback implementation for environments
    where scipy is not available. The formula used is:

        G1 = (sqrt(n*(n-1)) / (n-2)) * m3 / m2^(3/2)

    where m3 is the third central moment and m2 is the variance.
    This is the adjusted Fisher-Pearson coefficient.

    Returns NaN if n < 3 (insufficient data for bias correction).
    """
    n = len(arr)

    # Bias-corrected skewness requires at least 3 observations
    if n < 3:
        return float("nan"), "not_computed"

    mean = np.mean(arr)
    deviations = arr - mean

    # Second and third central moments
    m2 = np.mean(deviations**2)
    m3 = np.mean(deviations**3)

    # Guard against zero variance
    if m2 == 0:
        return float("nan"), "not_computed"

    # Unadjusted (biased) skewness
    g1 = m3 / (m2 ** (3 / 2))

    # Bias correction factor (Fisher-Pearson adjustment)
    correction = math.sqrt(n * (n - 1)) / (n - 2)
    skew_val = correction * g1

    return float(skew_val), "numpy_fallback"


def compute_skewness(arr: np.ndarray) -> Tuple[float, str]:
    """
    Compute descriptive skewness with automatic method selection.

    Attempts to use scipy.stats.skew (bias-corrected) if available.
    Falls back to NumPy-only implementation if scipy is unavailable.
    Returns NaN if skewness cannot be computed.

    Parameters
    ----------
    arr : np.ndarray
        Validated input array.

    Returns
    -------
    Tuple[float, str]
        Tuple of (skewness_value, method_used).
        method_used is one of: "scipy_bias_corrected", "numpy_fallback", "not_computed"

    Notes
    -----
    DIAGNOSTIC ONLY. Skewness is reported as a descriptive measure of
    distributional asymmetry. NO INFERENCE is made. NO INTERPRETATION
    is provided. The value should be used for visual inspection and
    structural characterization only.
    """
    n = len(arr)

    # Cannot compute skewness for very small samples
    if n < 3:
        return float("nan"), "not_computed"

    # Check for zero variance (all identical values)
    if np.var(arr) == 0:
        return float("nan"), "not_computed"

    if _SCIPY_AVAILABLE:
        return _compute_skewness_scipy(arr)
    else:
        return _compute_skewness_numpy(arr)


# =============================================================================
# DESCRIPTIVE STATISTICS COMPUTATION
# =============================================================================
def compute_descriptive_statistics(arr: np.ndarray) -> DescriptiveStatistics:
    """
    Compute purely descriptive statistics for validated data.

    This function computes structural characteristics of the data
    distribution WITHOUT any inference, hypothesis testing, or
    biological interpretation.

    Parameters
    ----------
    arr : np.ndarray
        Validated input array (must pass validate_input first).

    Returns
    -------
    DescriptiveStatistics
        Frozen dataclass containing all computed statistics and metadata.

    Notes
    -----
    DIAGNOSTIC ONLY. These statistics describe what IS in the data.
    They do NOT imply what SHOULD be, what is CORRECT, or what is
    SIGNIFICANT.

    Statistics computed:
        - count: Number of observations
        - minimum: Smallest value
        - maximum: Largest value
        - mean: Arithmetic mean
        - median: 50th percentile
        - std: Sample standard deviation (ddof=1)
        - skewness: Descriptive skewness (bias-corrected if scipy available)

    Examples
    --------
    >>> arr = validate_input([100, 150, 200, 250, 300])
    >>> stats = compute_descriptive_statistics(arr)
    >>> stats.mean
    200.0
    """
    n = len(arr)

    # Compute skewness with fallback
    skewness_val, skewness_method = compute_skewness(arr)

    return DescriptiveStatistics(
        count=n,
        minimum=float(np.min(arr)),
        maximum=float(np.max(arr)),
        mean=float(np.mean(arr)),
        median=float(np.median(arr)),
        std=float(np.std(arr, ddof=1)) if n > 1 else float("nan"),
        skewness=skewness_val,
        skewness_method=skewness_method,
        metadata=DiagnosticMetadata(),
    )


# =============================================================================
# HISTOGRAM GENERATION
# =============================================================================
def _determine_bin_count(n: int) -> int:
    """
    Determine appropriate bin count using Sturges' rule as baseline.

    Parameters
    ----------
    n : int
        Number of observations.

    Returns
    -------
    int
        Number of bins to use (minimum 5, maximum 100).

    Notes
    -----
    This is a DIAGNOSTIC choice for visualization only.
    It does NOT imply an optimal or correct binning strategy.
    """
    # Sturges' rule: k = ceil(log2(n) + 1)
    if n <= 1:
        return 5

    sturges = int(math.ceil(math.log2(n) + 1))

    # Clamp to reasonable range
    return max(5, min(sturges, 100))


def generate_histogram(
    arr: np.ndarray,
    *,
    log_mode: LogMode = None,
    bins: Optional[int] = None,
    title_suffix: str = "",
    figsize: Tuple[float, float] = (8, 5),
) -> Figure:
    """
    Generate a diagnostic histogram for visual inspection.

    This function creates a single histogram figure with configurable
    scaling for diagnostic purposes. It supports both linear and
    logarithmic visualization strategies.

    Parameters
    ----------
    arr : np.ndarray
        Validated input array.
    log_mode : LogMode, optional
        Logarithmic scaling strategy:
            - None: Linear scale (no log transformation)
            - "transform": Apply log10 transformation to data values
            - "axis": Plot raw data with logarithmic x-axis scale
        Default is None (linear scale).
    bins : int, optional
        Number of histogram bins. If None, automatically determined
        using Sturges' rule. Default is None.
    title_suffix : str, optional
        Additional text to append to the figure title. Default is "".
    figsize : Tuple[float, float], optional
        Figure size in inches (width, height). Default is (8, 5).

    Returns
    -------
    Figure
        Matplotlib Figure object. The figure is NOT saved to disk.
        The caller is responsible for display or serialization.

    Raises
    ------
    ValueError
        If log_mode is not one of: None, "transform", "axis".
    RuntimeError
        If log transformation fails (e.g., due to non-positive values
        after floating-point issues).

    Notes
    -----
    DIAGNOSTIC ONLY. This histogram is for VISUAL INSPECTION of
    distributional structure. It does NOT imply any statistical
    inference or biological interpretation.

    The histogram shows:
        - Data distribution shape
        - Approximate range and concentration
        - Potential structural features (for inspection only)

    NO assumptions are made about the underlying distribution.
    NO fit lines or theoretical distributions are overlaid.

    Examples
    --------
    >>> arr = validate_input([100, 150, 200, 250, 300, 350, 400])
    >>> fig = generate_histogram(arr, log_mode=None)
    >>> type(fig)
    <class 'matplotlib.figure.Figure'>
    """
    # -------------------------------------------------------------------------
    # Validate log_mode parameter
    # -------------------------------------------------------------------------
    valid_log_modes = (None, "transform", "axis")
    if log_mode not in valid_log_modes:
        raise ValueError(
            f"Invalid log_mode: {log_mode!r}. "
            f"Valid options are: {valid_log_modes}."
        )

    # -------------------------------------------------------------------------
    # Prepare data based on log_mode
    # -------------------------------------------------------------------------
    if log_mode == "transform":
        # Apply log10 transformation to data
        # Data was validated as strictly positive, but guard against edge cases
        if np.any(arr <= 0):
            raise RuntimeError(
                "Log transformation failed: data contains non-positive values. "
                "This should not occur with validated input."
            )
        plot_data = np.log10(arr)
        x_label = "log₁₀(Protein Length)"
        scale_note = "Log₁₀-Transformed"
    elif log_mode == "axis":
        # Keep raw data, will set axis to log scale
        plot_data = arr
        x_label = "Protein Length"
        scale_note = "Log-Scaled Axis"
    else:
        # Linear scale
        plot_data = arr
        x_label = "Protein Length"
        scale_note = "Linear Scale"

    # -------------------------------------------------------------------------
    # Determine bin count
    # -------------------------------------------------------------------------
    if bins is None:
        bins = _determine_bin_count(len(arr))

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)

    # -------------------------------------------------------------------------
    # Generate histogram
    # -------------------------------------------------------------------------
    if log_mode == "axis":
        # For log-scaled axis, use logarithmically-spaced bins
        # to avoid visual artifacts
        bin_edges = np.logspace(
            np.log10(plot_data.min()),
            np.log10(plot_data.max()),
            bins + 1,
        )
        ax.hist(
            plot_data,
            bins=bin_edges,
            edgecolor="black",
            linewidth=0.5,
            alpha=0.7,
            color="#4878CF",
        )
        ax.set_xscale("log")
    else:
        ax.hist(
            plot_data,
            bins=bins,
            edgecolor="black",
            linewidth=0.5,
            alpha=0.7,
            color="#4878CF",
        )

    # -------------------------------------------------------------------------
    # Labels and title
    # -------------------------------------------------------------------------
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel("Frequency", fontsize=10)

    title = f"Diagnostic Histogram ({scale_note})"
    if title_suffix:
        title = f"{title} - {title_suffix}"
    ax.set_title(title, fontsize=11, fontweight="medium")

    # -------------------------------------------------------------------------
    # Add diagnostic disclaimer annotation
    # -------------------------------------------------------------------------
    disclaimer = "DIAGNOSTIC ONLY • No inference • No interpretation"
    fig.text(
        0.99,
        0.01,
        disclaimer,
        ha="right",
        va="bottom",
        fontsize=7,
        style="italic",
        color="gray",
        transform=fig.transFigure,
    )

    # -------------------------------------------------------------------------
    # Finalize layout
    # -------------------------------------------------------------------------
    fig.tight_layout()

    return fig


# =============================================================================
# TOP-LEVEL ORCHESTRATION FUNCTION
# =============================================================================
def run_diagnostics(
    data: ArrayLike,
    *,
    log_mode: LogMode = "transform",
    bins: Optional[int] = None,
    title_suffix: str = "",
    figsize: Tuple[float, float] = (8, 5),
    include_log_histogram: bool = True,
) -> DiagnosticResult:
    """
    Run complete diagnostic analysis on protein length data.

    This is the primary entry point for Week 2 diagnostic validation.
    It performs input validation, computes descriptive statistics,
    and generates diagnostic histograms in a single call.

    Parameters
    ----------
    data : ArrayLike
        Input protein length data as Python list, NumPy array, or
        pandas Series. Values must be numeric, finite, and strictly
        positive.
    log_mode : LogMode, optional
        Logarithmic scaling strategy for the log histogram:
            - "transform": Apply log10 transformation to data (default)
            - "axis": Plot raw data with logarithmic x-axis
            - None: Skip log histogram entirely
        Default is "transform" for explicit, safe behavior.
    bins : int, optional
        Number of histogram bins. If None, automatically determined.
        Default is None.
    title_suffix : str, optional
        Additional text for figure titles. Default is "".
    figsize : Tuple[float, float], optional
        Figure size in inches. Default is (8, 5).
    include_log_histogram : bool, optional
        Whether to generate a log-scale histogram. If False, only the
        linear histogram is generated regardless of log_mode. Default is True.

    Returns
    -------
    DiagnosticResult
        Complete diagnostic result containing:
            - statistics: DescriptiveStatistics object
            - figure_linear: Matplotlib Figure (linear scale)
            - figure_log: Matplotlib Figure (log scale) or None
            - log_mode_used: The log mode that was applied
            - metadata: DiagnosticMetadata affirming diagnostic-only status

    Raises
    ------
    ValidationError
        If input data fails validation requirements.
    TypeError
        If input type is not supported.
    ValueError
        If log_mode is invalid.

    Notes
    -----
    DIAGNOSTIC ONLY. This function provides structural characterization
    of protein length data for Week 2 validation. It does NOT:
        - Perform hypothesis testing
        - Fit statistical models
        - Provide biological interpretation
        - Make assumptions about correctness or significance

    All outputs include explicit metadata affirming their diagnostic-only
    nature. This metadata should be preserved in downstream pipeline
    stages to maintain audit trail integrity.

    The function is safe to run in:
        - Local Python/Conda environments
        - CI pipelines
        - Docker containers

    It does NOT:
        - Save files to disk
        - Access environment variables
        - Depend on global state
        - Log sensitive data

    Examples
    --------
    >>> lengths = [150, 200, 180, 220, 190, 210, 175, 195, 205, 185]
    >>> result = run_diagnostics(lengths)
    >>> result.statistics.count
    10
    >>> result.statistics.metadata.diagnostic_only
    True
    >>> type(result.figure_linear)
    <class 'matplotlib.figure.Figure'>
    """
    # =========================================================================
    # STEP 1: INPUT VALIDATION (fail fast)
    # =========================================================================
    arr = validate_input(data)

    # =========================================================================
    # STEP 2: COMPUTE DESCRIPTIVE STATISTICS
    # =========================================================================
    statistics = compute_descriptive_statistics(arr)

    # =========================================================================
    # STEP 3: GENERATE LINEAR HISTOGRAM
    # =========================================================================
    figure_linear = generate_histogram(
        arr,
        log_mode=None,
        bins=bins,
        title_suffix=title_suffix,
        figsize=figsize,
    )

    # =========================================================================
    # STEP 4: GENERATE LOG HISTOGRAM (if requested)
    # =========================================================================
    figure_log: Optional[Figure] = None
    effective_log_mode: Optional[str] = None

    if include_log_histogram and log_mode is not None:
        figure_log = generate_histogram(
            arr,
            log_mode=log_mode,
            bins=bins,
            title_suffix=title_suffix,
            figsize=figsize,
        )
        effective_log_mode = log_mode

    # =========================================================================
    # STEP 5: ASSEMBLE AND RETURN RESULT
    # =========================================================================
    return DiagnosticResult(
        statistics=statistics,
        figure_linear=figure_linear,
        figure_log=figure_log,
        log_mode_used=effective_log_mode,
        metadata=DiagnosticMetadata(),
    )


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================
def get_statistics_dict(data: ArrayLike) -> Dict[str, Any]:
    """
    Compute descriptive statistics and return as dictionary.

    Convenience function for cases where only statistics are needed
    without histogram generation.

    Parameters
    ----------
    data : ArrayLike
        Input protein length data.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing all descriptive statistics and metadata.

    Notes
    -----
    DIAGNOSTIC ONLY. NO INFERENCE. NO INTERPRETATION.
    """
    arr = validate_input(data)
    stats = compute_descriptive_statistics(arr)
    return stats.to_dict()


def check_scipy_availability() -> Dict[str, Any]:
    """
    Check availability of optional scipy dependency.

    This function is provided for pipeline introspection and debugging.
    It allows callers to verify which skewness computation method will
    be used.

    Returns
    -------
    Dict[str, Any]
        Dictionary with:
            - scipy_available: bool
            - skewness_method: str describing which method will be used
            - note: str with diagnostic information

    Notes
    -----
    This is a DIAGNOSTIC utility function. It does not perform any
    statistical analysis.
    """
    return {
        "scipy_available": _SCIPY_AVAILABLE,
        "skewness_method": (
            "scipy_bias_corrected" if _SCIPY_AVAILABLE else "numpy_fallback"
        ),
        "note": (
            "scipy.stats.skew will be used with bias=False"
            if _SCIPY_AVAILABLE
            else "NumPy-only bias-corrected implementation will be used"
        ),
    }


def check_pandas_availability() -> Dict[str, Any]:
    """
    Check availability of optional pandas dependency.

    This function is provided for pipeline introspection and debugging.
    It allows callers to verify whether pandas Series input is supported.

    Returns
    -------
    Dict[str, Any]
        Dictionary with:
            - pandas_available: bool
            - series_input_supported: bool
            - note: str with diagnostic information

    Notes
    -----
    This is a DIAGNOSTIC utility function. It does not perform any
    statistical analysis.
    """
    return {
        "pandas_available": _PANDAS_AVAILABLE,
        "series_input_supported": _PANDAS_AVAILABLE,
        "note": (
            "pandas.Series input is supported"
            if _PANDAS_AVAILABLE
            else "pandas is not available; use list or numpy.ndarray input"
        ),
    }


# =============================================================================
# MODULE-LEVEL EXPORTS
# =============================================================================
__all__ = [
    # Primary entry point
    "run_diagnostics",
    # Result containers
    "DiagnosticResult",
    "DescriptiveStatistics",
    "DiagnosticMetadata",
    # Core functions
    "validate_input",
    "compute_descriptive_statistics",
    "compute_skewness",
    "generate_histogram",
    # Convenience functions
    "get_statistics_dict",
    "check_scipy_availability",
    "check_pandas_availability",
    # Exceptions
    "ValidationError",
    # Type aliases
    "ArrayLike",
    "LogMode",
]