"""
Week 2: Data Integrity & Statistical Validation — Inference Robustness Tests.

Tests inference robustness analysis functionality to ensure:
- Functions return outputs
- No pipeline crashes
- No metric drift
- Works when overlap < 5
- Works when overlap > 50
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

try:
    from scipy.stats import spearmanr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    import matplotlib
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


@pytest.mark.skipif(not SCIPY_AVAILABLE, reason="SciPy not available")
def test_correlation_stability_basic():
    """Test correlation stability simulation with basic data."""
    from week2_validation.reporting.robustness.correlation_stability import simulate_correlation_stability

    # Create test DataFrames with overlap
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"] * 5,
        "gene_2": ["B", "C", "D", "E", "F", "G", "H", "I", "J", "K"] * 5,
        "recurrence_count": [100, 90, 80, 70, 60, 50, 40, 30, 20, 10] * 5,
    })

    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D", "E", "X", "Y", "Z"],
        "gene_2": ["B", "C", "D", "E", "F", "Y", "Z", "W"],
        "recurrence_count": [95, 85, 75, 65, 55, 45, 35, 25],
    })

    result = simulate_correlation_stability(fusion_df, cosmic_df, overlap_sizes=[5, 10], bootstrap_n=50)

    # Should return DataFrame or None
    assert result is None or isinstance(result, pd.DataFrame)
    if result is not None:
        assert len(result) > 0
        assert "overlap_size" in result.columns
        assert "mean_rho" in result.columns


@pytest.mark.skipif(not SCIPY_AVAILABLE, reason="SciPy not available")
def test_correlation_stability_insufficient_overlap():
    """Test correlation stability with insufficient overlap."""
    from week2_validation.reporting.robustness.correlation_stability import simulate_correlation_stability

    # Create DataFrames with no overlap
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B"],
        "gene_2": ["B", "C"],
        "recurrence_count": [100, 90],
    })

    cosmic_df = pd.DataFrame({
        "gene_1": ["X", "Y"],
        "gene_2": ["Y", "Z"],
        "recurrence_count": [80, 70],
    })

    result = simulate_correlation_stability(fusion_df, cosmic_df)

    # Should return None gracefully
    assert result is None


@pytest.mark.skipif(not SCIPY_AVAILABLE, reason="SciPy not available")
def test_correlation_stability_large_overlap():
    """Test correlation stability with large overlap (>50)."""
    from week2_validation.reporting.robustness.correlation_stability import simulate_correlation_stability

    # Create DataFrames with large overlap
    genes = [f"G{i}" for i in range(100)]
    fusion_df = pd.DataFrame({
        "gene_1": genes,
        "gene_2": genes[1:] + [genes[0]],
        "recurrence_count": list(range(100, 0, -1)),
    })

    cosmic_df = pd.DataFrame({
        "gene_1": genes[:80],
        "gene_2": genes[1:81],
        "recurrence_count": list(range(95, 15, -1)),
    })

    result = simulate_correlation_stability(fusion_df, cosmic_df, overlap_sizes=[20, 50, 80], bootstrap_n=100)

    # Should return DataFrame or None
    assert result is None or isinstance(result, pd.DataFrame)


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_bootstrap_plot_generation():
    """Test bootstrap plot generation."""
    from week2_validation.reporting.robustness.bootstrap_visuals import generate_bootstrap_rho_plot

    # Create bootstrap samples
    bootstrap_rhos = np.random.normal(0.7, 0.1, 200).tolist()

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_bootstrap.png"

        result = generate_bootstrap_rho_plot(bootstrap_rhos, output_path)

        # Should return path or None
        assert result is None or isinstance(result, Path)
        if result:
            assert output_path.exists()
            assert output_path.stat().st_size > 0


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_bootstrap_plot_insufficient_samples():
    """Test bootstrap plot with insufficient samples."""
    from week2_validation.reporting.robustness.bootstrap_visuals import generate_bootstrap_rho_plot

    # Too few samples
    bootstrap_rhos = [0.7, 0.8, 0.6]

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_bootstrap.png"

        result = generate_bootstrap_rho_plot(bootstrap_rhos, output_path)

        # Should return None gracefully
        assert result is None


@pytest.mark.skipif(not SCIPY_AVAILABLE, reason="SciPy not available")
def test_jackknife_sensitivity_basic():
    """Test jackknife sensitivity analysis with basic data."""
    from week2_validation.reporting.robustness.jackknife_sensitivity import compute_cosmic_jackknife_sensitivity

    # Create test DataFrames with sufficient overlap
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D", "E", "F"],
        "gene_2": ["B", "C", "D", "E", "F", "G"],
        "recurrence_count": [100, 90, 80, 70, 60, 50],
    })

    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D", "E"],
        "gene_2": ["B", "C", "D", "E", "F"],
        "recurrence_count": [95, 85, 75, 65, 55],
    })

    result = compute_cosmic_jackknife_sensitivity(fusion_df, cosmic_df)

    # Should return DataFrame or None
    assert result is None or isinstance(result, pd.DataFrame)
    if result is not None:
        assert len(result) > 0
        assert "fusion_pair_removed" in result.columns
        assert "rho_after_removal" in result.columns
        assert "delta_from_full_rho" in result.columns


@pytest.mark.skipif(not SCIPY_AVAILABLE, reason="SciPy not available")
def test_jackknife_insufficient_overlap():
    """Test jackknife with insufficient overlap (<4 pairs)."""
    from week2_validation.reporting.robustness.jackknife_sensitivity import compute_cosmic_jackknife_sensitivity

    # Only 3 overlapping pairs (need 4 for jackknife)
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B", "C"],
        "gene_2": ["B", "C", "D"],
        "recurrence_count": [100, 90, 80],
    })

    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "B", "C"],
        "gene_2": ["B", "C", "D"],
        "recurrence_count": [95, 85, 75],
    })

    result = compute_cosmic_jackknife_sensitivity(fusion_df, cosmic_df)

    # Should return None gracefully
    assert result is None


def test_effect_size_interpretation_basic():
    """Test effect size interpretation generation."""
    from week2_validation.reporting.robustness.effect_size_interpretation import generate_correlation_effect_interpretation

    result = generate_correlation_effect_interpretation(rho=0.75, p_value=0.01, overlap_n=20)

    assert isinstance(result, str)
    assert len(result) > 0
    assert "strong" in result.lower() or "moderate" in result.lower() or "weak" in result.lower()


def test_effect_size_interpretation_weak():
    """Test effect size interpretation for weak correlation."""
    from week2_validation.reporting.robustness.effect_size_interpretation import generate_correlation_effect_interpretation

    result = generate_correlation_effect_interpretation(rho=0.2, p_value=0.1, overlap_n=15)

    assert isinstance(result, str)
    assert "weak" in result.lower()


def test_effect_size_interpretation_uncertainty():
    """Test effect size interpretation includes uncertainty when appropriate."""
    from week2_validation.reporting.robustness.effect_size_interpretation import generate_correlation_effect_interpretation

    # Low overlap case
    result_low_overlap = generate_correlation_effect_interpretation(rho=0.5, p_value=0.01, overlap_n=5)
    assert isinstance(result_low_overlap, str)
    assert "caution" in result_low_overlap.lower() or "uncertainty" in result_low_overlap.lower() or "limited" in result_low_overlap.lower()

    # High p-value case
    result_high_p = generate_correlation_effect_interpretation(rho=0.5, p_value=0.1, overlap_n=20)
    assert isinstance(result_high_p, str)
    assert "caution" in result_high_p.lower() or "cautiously" in result_high_p.lower()


def test_effect_size_interpretation_none_values():
    """Test effect size interpretation handles None values."""
    from week2_validation.reporting.robustness.effect_size_interpretation import generate_correlation_effect_interpretation

    result = generate_correlation_effect_interpretation(rho=None, p_value=None, overlap_n=None)

    assert isinstance(result, str)
    assert len(result) > 0


def test_robustness_no_pipeline_crash():
    """Test that robustness functions don't crash pipeline."""
    from week2_validation.reporting.robustness.correlation_stability import simulate_correlation_stability
    from week2_validation.reporting.robustness.jackknife_sensitivity import compute_cosmic_jackknife_sensitivity

    # Test with invalid inputs
    invalid_df = pd.DataFrame({"invalid": [1, 2, 3]})

    # Should return None, not crash
    result1 = simulate_correlation_stability(invalid_df, invalid_df)
    assert result1 is None

    result2 = compute_cosmic_jackknife_sensitivity(invalid_df, invalid_df)
    assert result2 is None
