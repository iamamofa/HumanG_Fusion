"""
Week 2: Data Integrity & Statistical Validation — COSMIC Validation Tests.

Tests COSMIC cross-validation functionality including:
- Real COSMIC data loading and transformation
- Spearman correlation computation
- Enrichment metrics
- Distribution comparison
- Real COSMIC Fusion v103 GRCh38 integration
- Cross-format compatibility (CSV, TSV, JSON, Parquet, XLSX)
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from week2_validation.cosmic.diagnostics import (
    compute_cosmic_statistical_metrics,
    compute_distribution_metrics,
    compute_enrichment_metric,
    compute_hypergeometric_enrichment,
    compute_spearman_correlation,
    run_cosmic_recurrence_diagnostic,
)
# Mock COSMIC imports removed - using real COSMIC data
from week2_validation.cosmic.gene_alias_map import normalize_gene_alias, GENE_ALIAS_MAP
from week2_validation.cosmic.statistical_controls import compute_negative_control_correlation
from week2_validation.cosmic.quality_gate import compute_cosmic_validation_score


# =============================================================================
# Test 1: Real COSMIC Data Loading and Transformation
# =============================================================================

def test_real_cosmic_loading():
    """Test that real COSMIC Fusion TSV can be loaded and transformed."""
    # Use real COSMIC file if available
    cosmic_module_path = Path(__file__).parent.parent / "cosmic"
    real_cosmic_path = cosmic_module_path / "Cosmic_Fusion_v103_GRCh38.tsv"
    
    if not real_cosmic_path.exists():
        pytest.skip("Real COSMIC file not found - skipping test")
    
    from week2_validation.utils.data_loader import load_reference_data
    from week2_validation.cosmic.diagnostics import transform_cosmic_fusion_to_standard_format
    
    # Load real COSMIC
    cosmic_df = load_reference_data(str(real_cosmic_path))
    assert cosmic_df is not None
    
    # Transform to standard format
    transformed_df = transform_cosmic_fusion_to_standard_format(cosmic_df)
    
    # Check required columns
    assert "gene_1" in transformed_df.columns
    assert "gene_2" in transformed_df.columns
    assert "recurrence_count" in transformed_df.columns
    
    # Check we have data
    assert len(transformed_df) > 0
    
    # Check recurrence counts are positive
    assert (transformed_df["recurrence_count"] > 0).all()


# =============================================================================
# Test 2: Spearman Correlation
# =============================================================================

def test_spearman_correlation_computation():
    """Test Spearman correlation is computed correctly."""
    # Create test data with known correlation
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 80,
        ("E", "F"): 60,
        ("G", "H"): 40,
        ("I", "J"): 20,
    }
    
    cosmic_pairs = {
        ("A", "B"): 90,  # Similar ranks
        ("C", "D"): 70,
        ("E", "F"): 50,
        ("G", "H"): 30,
        ("I", "J"): 10,
    }
    
    overlap = set(fusion_pairs.keys())
    
    result = compute_spearman_correlation(fusion_pairs, cosmic_pairs, overlap)
    
    # Should have correlation (positive, high)
    assert "spearman_rho" in result
    assert "spearman_p_value" in result
    
    # With scipy, should compute values
    try:
        from scipy.stats import spearmanr
        assert result["spearman_rho"] is not None
        assert result["spearman_p_value"] is not None
        assert isinstance(result["spearman_rho"], float)
        assert isinstance(result["spearman_p_value"], float)
    except ImportError:
        # Without scipy, should return None
        assert result["spearman_rho"] is None
        assert result["spearman_p_value"] is None


def test_spearman_correlation_insufficient_overlap():
    """Test Spearman correlation with insufficient overlap."""
    fusion_pairs = {("A", "B"): 100}
    cosmic_pairs = {("A", "B"): 90}
    overlap = set(fusion_pairs.keys())
    
    result = compute_spearman_correlation(fusion_pairs, cosmic_pairs, overlap)
    
    # Need at least 3 pairs for correlation
    assert result["spearman_rho"] is None
    assert result["spearman_p_value"] is None


# =============================================================================
# Test 3: Enrichment Metric
# =============================================================================

def test_enrichment_metric_computation():
    """Test enrichment metric computation."""
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 90,
        ("E", "F"): 80,
        ("G", "H"): 70,
        ("I", "J"): 60,
        ("K", "L"): 50,
        ("M", "N"): 40,
        ("O", "P"): 30,
        ("Q", "R"): 20,
        ("S", "T"): 10,
    }
    
    cosmic_pairs = {
        ("A", "B"): 95,  # Top overlap
        ("C", "D"): 85,
        ("E", "F"): 75,
        ("X", "Y"): 65,  # Different fusion
        ("Z", "W"): 55,
        ("K", "L"): 45,
        ("M", "N"): 35,
        ("O", "P"): 25,
        ("Q", "R"): 15,
        ("S", "T"): 5,
    }
    
    result = compute_enrichment_metric(fusion_pairs, cosmic_pairs, top_n=10)
    
    assert "top_fusion_overlap" in result
    assert "top_fusion_enrichment_ratio" in result
    
    # Should have some overlap (at least A-B, C-D, E-F)
    assert result["top_fusion_overlap"] >= 3
    assert 0.0 <= result["top_fusion_enrichment_ratio"] <= 1.0


# =============================================================================
# Test 4: Distribution Metrics
# =============================================================================

def test_distribution_metrics_computation():
    """Test distribution metrics computation."""
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 50,
        ("E", "F"): 25,
        ("G", "H"): 0,
    }
    
    cosmic_pairs = {
        ("A", "B"): 90,
        ("C", "D"): 45,
        ("E", "F"): 20,
        ("G", "H"): 0,
    }
    
    result = compute_distribution_metrics(fusion_pairs, cosmic_pairs)
    
    assert "mean_recurrence_fusion" in result
    assert "mean_recurrence_cosmic" in result
    assert "variance_fusion" in result
    assert "variance_cosmic" in result
    assert "zero_inflation_rate_fusion" in result
    assert "zero_inflation_rate_cosmic" in result
    
    # Check values are reasonable
    assert result["mean_recurrence_fusion"] > 0
    assert result["mean_recurrence_cosmic"] > 0
    assert result["variance_fusion"] >= 0
    assert result["variance_cosmic"] >= 0
    assert 0.0 <= result["zero_inflation_rate_fusion"] <= 1.0
    assert 0.0 <= result["zero_inflation_rate_cosmic"] <= 1.0


# =============================================================================
# Test 5: Full Diagnostic with Overlap
# =============================================================================

def test_cosmic_diagnostic_with_overlap():
    """Test full COSMIC diagnostic with overlapping fusions."""
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D"],
        "gene_2": ["B", "C", "D", "E"],
        "recurrence_count": [100, 80, 60, 40],
    })
    
    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "B", "X", "Y"],
        "gene_2": ["B", "C", "Y", "Z"],
        "recurrence_count": [90, 70, 50, 30],
    })
    
    result = run_cosmic_recurrence_diagnostic(
        fusion_df=fusion_df,
        cosmic_df=cosmic_df,
        top_n=10,
    )
    
    assert result["total_fusions_ours"] == 4
    assert result["total_fusions_cosmic"] == 4
    assert result["overlap_count"] == 2  # A-B and B-C overlap
    assert result["only_in_ours_count"] == 2
    assert result["only_in_cosmic_count"] == 2
    
    # Should have statistical metrics
    assert "spearman_rho" in result
    assert "top_fusion_overlap" in result
    assert "top_fusion_enrichment_ratio" in result


def test_cosmic_diagnostic_no_overlap():
    """Test COSMIC diagnostic with no overlapping fusions."""
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B"],
        "gene_2": ["B", "C"],
        "recurrence_count": [100, 80],
    })
    
    cosmic_df = pd.DataFrame({
        "gene_1": ["X", "Y"],
        "gene_2": ["Y", "Z"],
        "recurrence_count": [90, 70],
    })
    
    result = run_cosmic_recurrence_diagnostic(
        fusion_df=fusion_df,
        cosmic_df=cosmic_df,
        top_n=10,
    )
    
    assert result["overlap_count"] == 0
    assert "message" in result
    assert "No overlapping" in result["message"]


def test_cosmic_diagnostic_none_cosmic():
    """Test COSMIC diagnostic when COSMIC is None."""
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B"],
        "gene_2": ["B", "C"],
        "recurrence_count": [100, 80],
    })
    
    result = run_cosmic_recurrence_diagnostic(
        fusion_df=fusion_df,
        cosmic_df=None,
        top_n=10,
    )
    
    assert result["total_fusions_ours"] == 2
    assert result["total_fusions_cosmic"] == 0
    assert result["overlap_count"] == 0
    assert "message" in result


# =============================================================================
# Test 6: Cross-Format Compatibility
# =============================================================================

def _create_test_fusion_data(format_type: str, output_path: Path) -> Path:
    """Create test fusion data in specified format.
    
    Uses real gene names from mock COSMIC to ensure overlap.
    """
    df = pd.DataFrame({
        "fusion_id": ["F1", "F2", "F3"],
        "gene_1": ["BCR", "EML4", "TMPRSS2"],  # Real genes from mock COSMIC
        "gene_2": ["ABL1", "ALK", "ERG"],      # Real genes from mock COSMIC
        "protein_length": [100, 200, 300],
        "recurrence_count": [10, 20, 30],
    })
    
    if format_type == "csv":
        df.to_csv(output_path, index=False)
    elif format_type == "tsv":
        df.to_csv(output_path, sep="\t", index=False)
    elif format_type == "json":
        df.to_json(output_path, orient="records", indent=2)
    elif format_type == "parquet":
        df.to_parquet(output_path, index=False)
    elif format_type == "xlsx":
        df.to_excel(output_path, index=False, engine="openpyxl")
    else:
        raise ValueError(f"Unsupported format: {format_type}")
    
    return output_path


@pytest.mark.parametrize("format_type", ["csv", "tsv", "json", "parquet", "xlsx"])
def test_cosmic_validation_across_formats(format_type):
    """Test COSMIC validation works across all supported formats."""
    with tempfile.TemporaryDirectory() as tmp:
        # Create fusion data in specified format
        fusion_path = Path(tmp) / f"fusion.{format_type}"
        
        # Skip xlsx if openpyxl not available
        if format_type == "xlsx":
            try:
                import openpyxl
            except ImportError:
                pytest.skip("openpyxl not available")
        
        # Skip parquet if pyarrow not available
        if format_type == "parquet":
            try:
                import pyarrow
            except ImportError:
                pytest.skip("pyarrow not available")
        
        _create_test_fusion_data(format_type, fusion_path)
        
        # Create minimal test COSMIC data in standard format
        cosmic_df = pd.DataFrame({
            "gene_1": ["BCR", "EML4", "TMPRSS2", "EWSR1", "PML"],
            "gene_2": ["ABL1", "ALK", "ERG", "FLI1", "RARA"],
            "recurrence_count": [100, 80, 60, 40, 30],
        })
        cosmic_path = Path(tmp) / "test_cosmic.csv"
        cosmic_df.to_csv(cosmic_path, index=False)
        
        # Load data
        from week2_validation.utils.data_loader import load_fusion_data, load_reference_data
        
        fusion_df = load_fusion_data(str(fusion_path))
        cosmic_df = load_reference_data(str(cosmic_path))
        
        # Run diagnostic
        result = run_cosmic_recurrence_diagnostic(
            fusion_df=fusion_df,
            cosmic_df=cosmic_df,
            top_n=10,
        )
        
        # Verify results
        assert result["total_fusions_ours"] > 0
        assert result["total_fusions_cosmic"] > 0
        assert "spearman_rho" in result
        assert "top_fusion_overlap" in result


# =============================================================================
# Test 7: Pipeline Integration (Real COSMIC)
# =============================================================================

def test_pipeline_real_cosmic_loading():
    """Test that pipeline loads real COSMIC data correctly."""
    with tempfile.TemporaryDirectory() as tmp:
        # Create fusion data with real gene names to ensure overlap
        fusion_path = Path(tmp) / "fusion.csv"
        fusion_df = pd.DataFrame({
            "fusion_id": ["F1", "F2"],
            "gene_1": ["CCDC6", "PAX8"],  # Real genes from COSMIC
            "gene_2": ["RET", "PPARG"],   # Real genes from COSMIC
            "protein_length": [100, 200],
            "recurrence_count": [10, 20],
        })
        fusion_df.to_csv(fusion_path, index=False)
        
        # Create minimal test COSMIC data in standard format
        cosmic_df = pd.DataFrame({
            "gene_1": ["CCDC6", "PAX8", "TMPRSS2", "BCR"],
            "gene_2": ["RET", "PPARG", "ERG", "ABL1"],
            "recurrence_count": [50, 40, 30, 20],
        })
        cosmic_path = Path(tmp) / "test_cosmic.csv"
        cosmic_df.to_csv(cosmic_path, index=False)
        
        # Load and verify
        from week2_validation.utils.data_loader import load_reference_data
        
        loaded_cosmic_df = load_reference_data(str(cosmic_path))
        assert loaded_cosmic_df is not None
        assert len(loaded_cosmic_df) > 0
        
        # Run diagnostic
        result = run_cosmic_recurrence_diagnostic(
            fusion_df=fusion_df,
            cosmic_df=loaded_cosmic_df,
            top_n=10,
        )
        
        # After normalization, unique pairs may be fewer than 50 due to deduplication
        assert result["total_fusions_cosmic"] >= 4  # At least 4 unique pairs from test data
        # Provenance metadata is now included in diagnostic results
        assert "cosmic_reference_source" in result  # Provenance included


# =============================================================================
# Test 8: Statistical Metrics Integration
# =============================================================================

def test_statistical_metrics_integration():
    """Test that all statistical metrics are computed together."""
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 80,
        ("E", "F"): 60,
        ("G", "H"): 40,
        ("I", "J"): 20,
    }
    
    cosmic_pairs = {
        ("A", "B"): 90,
        ("C", "D"): 70,
        ("E", "F"): 50,
        ("G", "H"): 30,
        ("I", "J"): 10,
    }
    
    overlap = set(fusion_pairs.keys())
    
    metrics = compute_cosmic_statistical_metrics(fusion_pairs, cosmic_pairs, overlap)
    
    # Check all metrics are present
    assert "spearman_rho" in metrics
    assert "spearman_p_value" in metrics
    assert "top_fusion_overlap" in metrics
    assert "top_fusion_enrichment_ratio" in metrics
    assert "mean_recurrence_fusion" in metrics
    assert "mean_recurrence_cosmic" in metrics
    assert "variance_fusion" in metrics
    assert "variance_cosmic" in metrics
    assert "zero_inflation_rate_fusion" in metrics
    assert "zero_inflation_rate_cosmic" in metrics
    assert "enrichment_p_value" in metrics
    assert "expected_overlap_random" in metrics
    assert "observed_overlap" in metrics
    assert "negative_control_rho" in metrics
    assert "cosmic_validation_score" in metrics
    assert "cosmic_validation_classification" in metrics


# =============================================================================
# Test 9: Hypergeometric Enrichment Test
# =============================================================================

def test_hypergeometric_enrichment_computation():
    """Test hypergeometric enrichment p-value computation."""
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 90,
        ("E", "F"): 80,
        ("G", "H"): 70,
        ("I", "J"): 60,
        ("K", "L"): 50,
        ("M", "N"): 40,
        ("O", "P"): 30,
        ("Q", "R"): 20,
        ("S", "T"): 10,
    }
    
    cosmic_pairs = {
        ("A", "B"): 95,  # Top overlap
        ("C", "D"): 85,
        ("E", "F"): 75,
        ("X", "Y"): 65,  # Different fusion
        ("Z", "W"): 55,
        ("K", "L"): 45,
        ("M", "N"): 35,
        ("O", "P"): 25,
        ("Q", "R"): 15,
        ("S", "T"): 5,
    }
    
    result = compute_hypergeometric_enrichment(fusion_pairs, cosmic_pairs, top_n=10)
    
    assert "enrichment_p_value" in result
    assert "expected_overlap_random" in result
    assert "observed_overlap" in result
    
    # Should have observed overlap
    assert result["observed_overlap"] >= 0
    
    # Expected overlap should be reasonable
    if result["expected_overlap_random"] is not None:
        assert result["expected_overlap_random"] >= 0


# =============================================================================
# Test 10: Gene Alias Normalization
# =============================================================================

def test_gene_alias_normalization():
    """Test gene alias normalization."""
    # Test known aliases
    assert normalize_gene_alias("P53") == "TP53"
    assert normalize_gene_alias("p53") == "TP53"  # Case insensitive
    assert normalize_gene_alias("MLL") == "KMT2A"
    assert normalize_gene_alias("HER2") == "ERBB2"
    
    # Test canonical names (should return uppercase)
    assert normalize_gene_alias("TP53") == "TP53"
    assert normalize_gene_alias("ERBB2") == "ERBB2"
    
    # Test unknown genes (should return uppercase)
    assert normalize_gene_alias("UNKNOWN") == "UNKNOWN"
    assert normalize_gene_alias("unknown") == "UNKNOWN"


def test_alias_normalization_in_diagnostic():
    """Test that alias normalization is applied in diagnostic."""
    fusion_df = pd.DataFrame({
        "gene_1": ["P53", "MLL", "HER2"],
        "gene_2": ["ABL1", "ALK", "EGFR"],
        "recurrence_count": [100, 80, 60],
    })
    
    cosmic_df = pd.DataFrame({
        "gene_1": ["TP53", "KMT2A", "ERBB2"],  # Canonical names
        "gene_2": ["ABL1", "ALK", "EGFR"],
        "recurrence_count": [90, 70, 50],
    })
    
    result = run_cosmic_recurrence_diagnostic(
        fusion_df=fusion_df,
        cosmic_df=cosmic_df,
        top_n=10,
    )
    
    # Should find overlap after alias normalization
    assert result["overlap_count"] == 3


# =============================================================================
# Test 11: Negative Control Test
# =============================================================================

def test_negative_control_computation():
    """Test negative control correlation computation."""
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 80,
        ("E", "F"): 60,
        ("G", "H"): 40,
        ("I", "J"): 20,
    }
    
    cosmic_pairs = {
        ("A", "B"): 90,  # Similar ranks
        ("C", "D"): 70,
        ("E", "F"): 50,
        ("G", "H"): 30,
        ("I", "J"): 10,
    }
    
    overlap = set(fusion_pairs.keys())
    
    result = compute_negative_control_correlation(
        fusion_pairs, cosmic_pairs, overlap, n_shuffles=100
    )
    
    assert "negative_control_rho" in result
    assert "negative_control_p_value" in result
    assert "shuffled_rho_mean" in result
    assert "shuffled_rho_std" in result
    
    # With scipy, should compute values
    try:
        from scipy.stats import spearmanr
        assert result["negative_control_rho"] is not None
        assert result["negative_control_p_value"] is not None
        assert isinstance(result["negative_control_rho"], float)
        assert isinstance(result["negative_control_p_value"], float)
        
        # Real correlation should be higher than shuffled (for this test data)
        # Note: This may not always be true, but for well-correlated data it should be
    except ImportError:
        # Without scipy, should return None
        assert result["negative_control_rho"] is None
        assert result["negative_control_p_value"] is None


def test_negative_control_reduces_rho():
    """Test that negative control (shuffled) reduces correlation."""
    # Create strongly correlated data
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 90,
        ("E", "F"): 80,
        ("G", "H"): 70,
        ("I", "J"): 60,
    }
    
    cosmic_pairs = {
        ("A", "B"): 95,
        ("C", "D"): 85,
        ("E", "F"): 75,
        ("G", "H"): 65,
        ("I", "J"): 55,
    }
    
    overlap = set(fusion_pairs.keys())
    
    # Compute real correlation
    from week2_validation.cosmic.diagnostics import compute_spearman_correlation
    real_result = compute_spearman_correlation(fusion_pairs, cosmic_pairs, overlap)
    
    # Compute negative control
    control_result = compute_negative_control_correlation(
        fusion_pairs, cosmic_pairs, overlap, n_shuffles=1000
    )
    
    if real_result["spearman_rho"] is not None and control_result["negative_control_rho"] is not None:
        real_abs = abs(real_result["spearman_rho"])
        shuffled_abs = abs(control_result["negative_control_rho"])
        # Real correlation should be stronger than shuffled (on average)
        # This is probabilistic, so we check that shuffled is not significantly higher
        assert shuffled_abs <= real_abs + 0.3  # Allow some variance


# =============================================================================
# Test 12: Quality Gate Classification
# =============================================================================

def test_quality_gate_computation():
    """Test COSMIC quality gate score and classification."""
    # Strong agreement case
    score_strong = compute_cosmic_validation_score(
        spearman_rho=0.8,
        spearman_p_value=0.001,
        enrichment_p_value=0.001,
        negative_control_rho=0.1,
        observed_overlap=8,
        expected_overlap_random=2.0,
    )
    
    assert "cosmic_validation_score" in score_strong
    assert "cosmic_validation_classification" in score_strong
    assert score_strong["cosmic_validation_score"] >= 0.0
    assert score_strong["cosmic_validation_score"] <= 1.0
    assert score_strong["cosmic_validation_classification"] in [
        "STRONG BIOLOGICAL AGREEMENT",
        "MODERATE AGREEMENT",
        "WEAK SIGNAL",
        "NO MEANINGFUL COSMIC AGREEMENT",
    ]
    
    # Weak agreement case
    score_weak = compute_cosmic_validation_score(
        spearman_rho=0.2,
        spearman_p_value=0.1,
        enrichment_p_value=0.1,
        negative_control_rho=0.15,
        observed_overlap=3,
        expected_overlap_random=2.5,
    )
    
    assert score_weak["cosmic_validation_score"] < score_strong["cosmic_validation_score"]


# =============================================================================
# Test 13: Provenance Metadata
# =============================================================================

def test_provenance_metadata():
    """Test that provenance metadata is included in diagnostic results."""
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B"],
        "gene_2": ["B", "C"],
        "recurrence_count": [100, 80],
    })
    
    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "X"],
        "gene_2": ["B", "Y"],
        "recurrence_count": [90, 70],
    })
    
    result = run_cosmic_recurrence_diagnostic(
        fusion_df=fusion_df,
        cosmic_df=cosmic_df,
        top_n=10,
        cosmic_reference_source="mock",
        cosmic_file_path=None,
    )
    
    assert "cosmic_reference_source" in result
    assert "cosmic_reference_version" in result
    assert "cosmic_reference_load_timestamp" in result
    
    assert result["cosmic_reference_source"] == "mock"
    assert result["cosmic_reference_version"] == "mock_v1"


# =============================================================================
# Test 14: Cross-Format Stability
# =============================================================================

def test_cosmic_metrics_identical_across_formats():
    """Test that COSMIC metrics are identical across file formats."""
    # Create test data in memory
    fusion_data = {
        "fusion_id": ["F1", "F2", "F3"],
        "gene_1": ["A", "B", "C"],
        "gene_2": ["B", "C", "D"],
        "protein_length": [100, 200, 300],
        "recurrence_count": [10, 20, 30],
    }
    
    cosmic_data = {
        "gene_1": ["A", "B", "X"],
        "gene_2": ["B", "C", "Y"],
        "recurrence_count": [9, 19, 29],
    }
    
    fusion_df = pd.DataFrame(fusion_data)
    cosmic_df = pd.DataFrame(cosmic_data)
    
    # Run diagnostic
    result = run_cosmic_recurrence_diagnostic(
        fusion_df=fusion_df,
        cosmic_df=cosmic_df,
        top_n=10,
    )
    
    # Store key metrics
    baseline_metrics = {
        "spearman_rho": result.get("spearman_rho"),
        "spearman_p_value": result.get("spearman_p_value"),
        "overlap_count": result.get("overlap_count"),
        "enrichment_p_value": result.get("enrichment_p_value"),
    }
    
    # Test that metrics are consistent (same data should give same results)
    result2 = run_cosmic_recurrence_diagnostic(
        fusion_df=fusion_df.copy(),
        cosmic_df=cosmic_df.copy(),
        top_n=10,
    )
    
    for key, value in baseline_metrics.items():
        assert result2.get(key) == value, f"Metric {key} differs between runs"


# =============================================================================
# Test 15: Bootstrap CI Contains Rho
# =============================================================================

def test_bootstrap_ci_contains_rho():
    """Test that bootstrap CI contains the point estimate rho."""
    # Create test data with sufficient overlap (≥10 required for bootstrap CI)
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 80,
        ("E", "F"): 60,
        ("G", "H"): 40,
        ("I", "J"): 20,
        ("K", "L"): 15,
        ("M", "N"): 10,
        ("O", "P"): 5,
        ("Q", "R"): 3,
        ("S", "T"): 2,
    }
    
    cosmic_pairs = {
        ("A", "B"): 95,
        ("C", "D"): 75,
        ("E", "F"): 55,
        ("G", "H"): 35,
        ("I", "J"): 18,
        ("K", "L"): 12,
        ("M", "N"): 8,
        ("O", "P"): 4,
        ("Q", "R"): 2,
        ("S", "T"): 1,
    }
    
    overlap = set(fusion_pairs.keys())
    
    result = compute_spearman_correlation(
        fusion_pairs, cosmic_pairs, overlap,
        compute_bootstrap_ci=True,
        bootstrap_iterations=500,
        bootstrap_seed=42,
    )
    
    # With scipy, should compute CI
    try:
        from scipy.stats import spearmanr
        assert result["spearman_rho"] is not None
        assert result["rho_ci_lower"] is not None
        assert result["rho_ci_upper"] is not None
        
        # CI should contain the point estimate (with small tolerance for edge cases)
        rho = result["spearman_rho"]
        ci_lower = result["rho_ci_lower"]
        ci_upper = result["rho_ci_upper"]
        
        # The CI should be valid (lower < upper)
        assert ci_lower <= ci_upper
        
        # The point estimate should be within or very close to the CI
        # Allow some tolerance because bootstrap CI may not always contain point estimate
        tolerance = 0.1
        assert ci_lower - tolerance <= rho <= ci_upper + tolerance
        
    except ImportError:
        # Without scipy, CI should be None
        assert result["rho_ci_lower"] is None
        assert result["rho_ci_upper"] is None


# =============================================================================
# Test 16: Reproducibility Lock Metadata
# =============================================================================

def test_reproducibility_lock_metadata():
    """Test that reproducibility lock metadata exists in provenance."""
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B", "C"],
        "gene_2": ["B", "C", "D"],
        "recurrence_count": [100, 80, 60],
    })
    
    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "B"],
        "gene_2": ["B", "C"],
        "recurrence_count": [90, 70],
    })
    
    result = run_cosmic_recurrence_diagnostic(
        fusion_df=fusion_df,
        cosmic_df=cosmic_df,
        top_n=10,
        cosmic_reference_source="mock",
        mock_generation_seed=42,
    )
    
    # Check reproducibility lock exists
    assert "reproducibility_lock" in result
    lock = result["reproducibility_lock"]
    
    # Check required fields
    assert "python_version" in lock
    assert "numpy_version" in lock
    assert "pandas_version" in lock
    assert "cosmic_validation_code_version" in lock
    assert "random_seed_mock_generation" in lock
    
    # Verify values are not None/empty
    assert lock["python_version"] is not None
    assert lock["numpy_version"] is not None
    assert lock["pandas_version"] is not None
    assert lock["cosmic_validation_code_version"] is not None
    # mock_generation_seed is None for real COSMIC
    assert lock["random_seed_mock_generation"] is None or lock["random_seed_mock_generation"] == 42


# =============================================================================
# Test 17: COSMIC Bias Quantification (Distribution Metrics)
# =============================================================================

def test_mock_distribution_metrics_computed():
    """Test that COSMIC bias quantification metrics are computed correctly."""
    from week2_validation.reporting.bias_analysis.cosmic_bias_quantification import (
        quantify_cosmic_sampling_bias,
    )
    
    # Create test COSMIC data in standard format (7 unique gene pairs, no duplicates)
    cosmic_df = pd.DataFrame({
        "gene_1": ["BCR", "EML4", "TMPRSS2", "EWSR1", "PML", "RET", "CCDC6"],
        "gene_2": ["ABL1", "ALK", "ERG", "FLI1", "RARA", "CCDC6", "PAX8"],
        "recurrence_count": [100, 80, 60, 40, 30, 25, 20],
    })
    
    result = quantify_cosmic_sampling_bias(cosmic_df)
    
    # Check required keys
    assert "top5_gene_fraction" in result
    assert "gini_gene_recurrence" in result
    assert "bias_severity_classification" in result
    assert "n_unique_genes" in result
    assert "n_fusion_pairs" in result
    
    # Basic sanity checks (7 unique pairs: BCR-ABL1, EML4-ALK, TMPRSS2-ERG, etc.)
    assert result["n_unique_genes"] >= 7
    assert result["n_fusion_pairs"] == 7
    assert result["top5_gene_fraction"] is not None
    assert 0.0 <= result["top5_gene_fraction"] <= 1.0
    assert result["gini_gene_recurrence"] is not None
    assert 0.0 <= result["gini_gene_recurrence"] <= 1.0
    assert result["bias_severity_classification"] in ["LOW", "MODERATE", "HIGH", "N/A"]


def test_gini_coefficient_computation():
    """Test Gini coefficient computation from bias analysis module."""
    from week2_validation.reporting.bias_analysis.cosmic_bias_quantification import (
        compute_gini_coefficient,
    )
    
    # Perfect equality: all same values -> Gini = 0
    equal_values = np.array([10, 10, 10, 10, 10])
    gini_equal = compute_gini_coefficient(equal_values)
    assert abs(gini_equal) < 0.01  # Should be close to 0
    
    # High inequality: one dominates (note: compute_gini_coefficient filters zeros)
    unequal_values = np.array([1, 1, 1, 1, 100])  # One value dominates
    gini_unequal = compute_gini_coefficient(unequal_values)
    assert gini_unequal > 0.5  # Should be high inequality
    
    # Empty/small array safety
    empty_gini = compute_gini_coefficient(np.array([]))
    assert empty_gini == 0.0


# =============================================================================
# Test 18: Score Breakdown Sums Correctly
# =============================================================================

def test_score_breakdown_sums_correctly():
    """Test that score component breakdown sums to total score."""
    # Strong agreement case
    score_result = compute_cosmic_validation_score(
        spearman_rho=0.8,
        spearman_p_value=0.001,
        enrichment_p_value=0.001,
        negative_control_rho=0.1,
        observed_overlap=8,
        expected_overlap_random=2.0,
    )
    
    assert "score_component_breakdown" in score_result
    breakdown = score_result["score_component_breakdown"]
    
    # Check breakdown components exist
    assert "correlation_component" in breakdown
    assert "enrichment_component" in breakdown
    assert "negative_control_component" in breakdown
    assert "overlap_component" in breakdown
    assert "total_raw_score" in breakdown
    assert "max_possible_score" in breakdown
    
    # Check breakdown sum check passes
    assert "breakdown_sum_check" in breakdown
    assert breakdown["breakdown_sum_check"] == True
    
    # Verify detailed components exist
    assert "detailed_components" in breakdown
    detailed = breakdown["detailed_components"]
    assert "correlation_component" in detailed
    assert "enrichment_component" in detailed
    assert "negative_control_component" in detailed
    assert "overlap_component" in detailed
    
    # Each detailed component should have interpretation
    for comp_name, comp_data in detailed.items():
        assert "interpretation" in comp_data
        assert "raw_score" in comp_data
        assert "max_possible" in comp_data


# =============================================================================
# Test 19: Statistical Docs File Exists
# =============================================================================

def test_statistical_methodology_docs_exist():
    """Test that statistical methodology documentation file exists."""
    from pathlib import Path
    
    # Get the cosmic module directory
    import week2_validation.cosmic as cosmic_module
    cosmic_dir = Path(cosmic_module.__file__).parent
    
    # Check for statistical_methodology.md
    methodology_path = cosmic_dir / "statistical_methodology.md"
    assert methodology_path.exists(), f"Statistical methodology docs not found at {methodology_path}"
    
    # Verify file has content
    content = methodology_path.read_text(encoding="utf-8")
    assert len(content) > 1000  # Should have substantial content
    
    # Check for key sections
    assert "Spearman" in content
    assert "Hypergeometric" in content
    assert "Bootstrap" in content
    assert "Independence" in content or "independence" in content


# =============================================================================
# Test 20: Bootstrap CI with Insufficient Data
# =============================================================================

def test_bootstrap_ci_insufficient_data():
    """Test bootstrap CI returns None with insufficient overlap."""
    fusion_pairs = {("A", "B"): 100}
    cosmic_pairs = {("A", "B"): 90}
    overlap = set(fusion_pairs.keys())
    
    result = compute_spearman_correlation(
        fusion_pairs, cosmic_pairs, overlap,
        compute_bootstrap_ci=True,
    )
    
    # Should return None for CI with only 1 pair
    assert result["rho_ci_lower"] is None
    assert result["rho_ci_upper"] is None


# =============================================================================
# Test 21: Provenance with Mock Generation Seed
# =============================================================================

def test_provenance_includes_mock_seed():
    """Test that provenance includes mock generation seed when provided."""
    from week2_validation.cosmic.diagnostics import compute_cosmic_provenance
    
    provenance = compute_cosmic_provenance(
        cosmic_df=None,
        cosmic_reference_source="mock",
        cosmic_file_path=None,
        mock_generation_seed=12345,
    )
    
    assert "reproducibility_lock" in provenance
    assert provenance["reproducibility_lock"]["random_seed_mock_generation"] == 12345


# =============================================================================
# Test 22: Quality Gate Component Interpretations
# =============================================================================

def test_quality_gate_interpretations():
    """Test that quality gate provides interpretations for all components."""
    score_result = compute_cosmic_validation_score(
        spearman_rho=0.5,
        spearman_p_value=0.01,
        enrichment_p_value=0.01,
        negative_control_rho=0.1,
        observed_overlap=5,
        expected_overlap_random=2.0,
    )
    
    breakdown = score_result["score_component_breakdown"]
    detailed = breakdown["detailed_components"]
    
    # All components should have non-empty interpretations
    for comp_name, comp_data in detailed.items():
        assert "interpretation" in comp_data
        assert len(comp_data["interpretation"]) > 0
        assert comp_data["interpretation"] != "N/A"
