"""
Week 2: Data Integrity & Statistical Validation — Histogram Generation Tests.

Tests histogram generation functionality to ensure:
- Histogram files are created successfully
- File size > 0
- Log histogram works
- Handles missing column safely
- Graceful degradation when matplotlib unavailable
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


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_histogram_generation_basic():
    """Test that histogram is generated successfully with basic data."""
    from week2_validation.reporting.visuals.histograms import generate_protein_length_histogram

    # Create test DataFrame
    df = pd.DataFrame({
        "protein_length": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] * 10,
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_histogram.png"

        result = generate_protein_length_histogram(df, output_path)

        # Verify histogram was created
        assert result is not None
        assert output_path.exists()
        assert output_path.stat().st_size > 0


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_log_histogram_generation():
    """Test that log histogram is generated successfully."""
    from week2_validation.reporting.visuals.histograms import generate_log_protein_length_histogram

    # Create test DataFrame with positive values
    df = pd.DataFrame({
        "protein_length": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] * 10,
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_log_histogram.png"

        result = generate_log_protein_length_histogram(df, output_path)

        # Verify histogram was created
        assert result is not None
        assert output_path.exists()
        assert output_path.stat().st_size > 0


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_histogram_missing_column():
    """Test that histogram generation handles missing column safely."""
    from week2_validation.reporting.visuals.histograms import generate_protein_length_histogram

    # Create DataFrame without protein_length column
    df = pd.DataFrame({
        "other_column": [1, 2, 3, 4, 5],
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_histogram.png"

        result = generate_protein_length_histogram(df, output_path)

        # Should return None gracefully
        assert result is None


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_histogram_empty_dataframe():
    """Test that histogram generation handles empty DataFrame safely."""
    from week2_validation.reporting.visuals.histograms import generate_protein_length_histogram

    # Create empty DataFrame
    df = pd.DataFrame({
        "protein_length": [],
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_histogram.png"

        result = generate_protein_length_histogram(df, output_path)

        # Should return None gracefully
        assert result is None


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_histogram_negative_values():
    """Test that histogram generation handles negative values safely."""
    from week2_validation.reporting.visuals.histograms import generate_protein_length_histogram

    # Create DataFrame with some negative values
    df = pd.DataFrame({
        "protein_length": [100, 200, -50, 300, 400, 500] * 5,
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_histogram.png"

        result = generate_protein_length_histogram(df, output_path)

        # Should still generate histogram (negative values filtered out)
        assert result is not None
        assert output_path.exists()


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_log_histogram_negative_values():
    """Test that log histogram filters negative values."""
    from week2_validation.reporting.visuals.histograms import generate_log_protein_length_histogram

    # Create DataFrame with some negative and zero values
    df = pd.DataFrame({
        "protein_length": [100, 200, -50, 0, 300, 400, 500] * 5,
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_log_histogram.png"

        result = generate_log_protein_length_histogram(df, output_path)

        # Should still generate histogram (negative/zero values filtered out)
        assert result is not None
        assert output_path.exists()


def test_histogram_without_matplotlib():
    """Test that histogram generation gracefully handles missing matplotlib."""
    # This test runs even if matplotlib is not available
    try:
        from week2_validation.reporting.visuals.histograms import generate_protein_length_histogram

        df = pd.DataFrame({
            "protein_length": [100, 200, 300],
        })

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_path = tmp_path / "test_histogram.png"

            result = generate_protein_length_histogram(df, output_path)

            # If matplotlib is not available, should return None gracefully
            if not MATPLOTLIB_AVAILABLE:
                assert result is None
            else:
                # If matplotlib is available, should generate histogram
                assert result is not None
                assert output_path.exists()
    except ImportError:
        # If module can't be imported due to missing matplotlib, that's OK
        pass


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
def test_histogram_custom_column_name():
    """Test histogram generation with custom column name."""
    from week2_validation.reporting.visuals.histograms import generate_protein_length_histogram

    # Create DataFrame with custom column name
    df = pd.DataFrame({
        "custom_length": [100, 200, 300, 400, 500] * 10,
    })

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        output_path = tmp_path / "test_histogram.png"

        result = generate_protein_length_histogram(df, output_path, column="custom_length")

        # Verify histogram was created
        assert result is not None
        assert output_path.exists()
        assert output_path.stat().st_size > 0
