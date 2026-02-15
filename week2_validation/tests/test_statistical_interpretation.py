"""
Week 2: Data Integrity & Statistical Validation — Statistical Interpretation Tests.

Tests statistical interpretation functionality to ensure:
- Interpretation text is generated correctly
- Handles missing metrics gracefully
- Skewness visualization file is created
- File size > 0
"""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

try:
    import matplotlib
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def test_interpretation_basic():
    """Test that interpretation text is generated from basic metrics."""
    from week2_validation.reporting.statistical_interpretation import generate_distribution_interpretation

    metrics = {
        "mean": 500.0,
        "median": 510.0,
        "skewness": 0.3,
    }

    result = generate_distribution_interpretation(metrics)

    # Should return a string
    assert isinstance(result, str)
    assert len(result) > 0


def test_interpretation_symmetric_distribution():
    """Test interpretation for symmetric distribution."""
    from week2_validation.reporting.statistical_interpretation import generate_distribution_interpretation

    metrics = {
        "mean": 500.0,
        "median": 500.0,
        "skewness": 0.1,
    }

    result = generate_distribution_interpretation(metrics)

    assert isinstance(result, str)
    assert "symmetric" in result.lower() or "symmetry" in result.lower()


def test_interpretation_skewed_distribution():
    """Test interpretation for skewed distribution."""
    from week2_validation.reporting.statistical_interpretation import generate_distribution_interpretation

    metrics = {
        "mean": 500.0,
        "median": 400.0,
        "skewness": 1.2,
    }

    result = generate_distribution_interpretation(metrics)

    assert isinstance(result, str)
    assert "skew" in result.lower() or "asymmetry" in result.lower()


def test_interpretation_missing_metrics():
    """Test that interpretation handles missing metrics gracefully."""
    from week2_validation.reporting.statistical_interpretation import generate_distribution_interpretation

    # Empty metrics
    result = generate_distribution_interpretation({})

    assert isinstance(result, str)
    assert len(result) > 0


def test_interpretation_partial_metrics():
    """Test interpretation with partial metrics."""
    from week2_validation.reporting.statistical_interpretation import generate_distribution_interpretation

    # Only skewness available
    metrics = {
        "skewness": 0.5,
    }

    result = generate_distribution_interpretation(metrics)

    assert isinstance(result, str)
    assert len(result) > 0


def test_interpretation_benford_not_applicable():
    """Test interpretation when Benford is not applicable."""
    from week2_validation.reporting.statistical_interpretation import generate_distribution_interpretation

    metrics = {
        "benford_applicable": False,
        "scale_span_orders_of_magnitude": 1.5,
    }

    result = generate_distribution_interpretation(metrics)

    assert isinstance(result, str)
    assert "benford" in result.lower() or "scale span" in result.lower()


def test_extract_interpretation_metrics():
    """Test metric extraction from status and diagnostic data."""
    from week2_validation.reporting.statistical_interpretation import extract_interpretation_metrics

    status_data = {
        "data_quality": {
            "skewness": 0.3,
        },
    }

    diagnostic_data = {
        "distribution": {
            "mean": 500.0,
            "median": 510.0,
            "std": 100.0,
            "skewness": 0.3,
        },
        "benford": {
            "applicability": True,
            "scale_span_orders_of_magnitude": 2.5,
        },
    }

    metrics = extract_interpretation_metrics(status_data, diagnostic_data)

    assert isinstance(metrics, dict)
    assert "mean" in metrics
    assert "median" in metrics
    assert "skewness" in metrics


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_skewness_visualization_generation():
    """Test that skewness visualization is generated successfully."""
    from week2_validation.reporting.visuals.skewness_plots import generate_skewness_visualization

    # Create test DataFrame
    df = pd.DataFrame({
        "protein_length": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] * 10,
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_skewness.png"

        result = generate_skewness_visualization(df, output_path)

        # Verify visualization was created
        assert result is not None
        assert output_path.exists()
        assert output_path.stat().st_size > 0


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_skewness_visualization_missing_column():
    """Test that skewness visualization handles missing column safely."""
    from week2_validation.reporting.visuals.skewness_plots import generate_skewness_visualization

    # Create DataFrame without protein_length column
    df = pd.DataFrame({
        "other_column": [1, 2, 3, 4, 5],
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_skewness.png"

        result = generate_skewness_visualization(df, output_path)

        # Should return None gracefully
        assert result is None


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_skewness_visualization_empty_dataframe():
    """Test that skewness visualization handles empty DataFrame safely."""
    from week2_validation.reporting.visuals.skewness_plots import generate_skewness_visualization

    # Create empty DataFrame
    df = pd.DataFrame({
        "protein_length": [],
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_skewness.png"

        result = generate_skewness_visualization(df, output_path)

        # Should return None gracefully
        assert result is None


def test_skewness_visualization_without_matplotlib():
    """Test that skewness visualization gracefully handles missing matplotlib."""
    # This test runs even if matplotlib is not available
    try:
        from week2_validation.reporting.visuals.skewness_plots import generate_skewness_visualization

        df = pd.DataFrame({
            "protein_length": [100, 200, 300],
        })

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_path = tmp_path / "test_skewness.png"

            result = generate_skewness_visualization(df, output_path)

            # If matplotlib is not available, should return None gracefully
            if not MATPLOTLIB_AVAILABLE:
                assert result is None
            else:
                # If matplotlib is available, should generate visualization
                assert result is not None
                assert output_path.exists()
    except ImportError:
        # If module can't be imported due to missing matplotlib, that's OK
        pass
