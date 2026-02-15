"""
Week 2: Data Integrity & Statistical Validation — Tier 2 Statistical Robustness Tests.

Tests for:
- Null model simulation
- External validity stability
- COSMIC bias quantification
- Scientific claim strength classifier
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from week2_validation.reporting.bias_analysis.cosmic_bias_quantification import (
    quantify_cosmic_sampling_bias,
)
from week2_validation.reporting.external_validity.distribution_generalization import (
    assess_distribution_generalization,
)
from week2_validation.reporting.robustness.null_model_simulation import (
    run_cosmic_null_model_test,
)
from week2_validation.reporting.robustness.scientific_claim_classifier import (
    classify_scientific_claim_strength,
)


# =============================================================================
# Test 1: Null Model Simulation
# =============================================================================

def test_null_model_returns_valid_p_values():
    """Test that null model returns valid p-values."""
    # Create test data with known correlation
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"],
        "gene_2": ["B", "C", "D", "E", "F", "G", "H", "I", "J", "K"],
        "recurrence_count": [100, 90, 80, 70, 60, 50, 40, 30, 20, 10],
    })
    
    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"],
        "gene_2": ["B", "C", "D", "E", "F", "G", "H", "I", "J", "K"],
        "recurrence_count": [95, 85, 75, 65, 55, 45, 35, 25, 15, 5],
    })
    
    result = run_cosmic_null_model_test(fusion_df, cosmic_df, n_simulations=100)
    
    assert "observed_rho" in result
    assert "null_mean_rho" in result
    assert "null_std_rho" in result
    assert "empirical_p_value" in result
    assert "z_score_vs_null" in result
    assert "n_simulations" in result
    
    # Check values are valid
    if result["observed_rho"] is not None:
        assert isinstance(result["observed_rho"], float)
        assert -1.0 <= result["observed_rho"] <= 1.0
    
    if result["empirical_p_value"] is not None:
        assert 0.0 <= result["empirical_p_value"] <= 1.0
    
    if result["null_mean_rho"] is not None:
        assert isinstance(result["null_mean_rho"], float)
    
    if result["null_std_rho"] is not None:
        assert result["null_std_rho"] >= 0.0


def test_null_model_insufficient_overlap():
    """Test null model with insufficient overlap."""
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
    
    result = run_cosmic_null_model_test(fusion_df, cosmic_df, n_simulations=100)
    
    # Should return None values for insufficient overlap
    assert result["observed_rho"] is None or result["n_simulations"] == 0


# =============================================================================
# Test 2: External Validity Stability
# =============================================================================

def test_external_validity_returns_stability_index():
    """Test that external validity returns stability index."""
    # Create test data with sufficient fusions for stratification
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"],
        "gene_2": ["B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"],
        "recurrence_count": [100, 90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 5],
    })
    
    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"],
        "gene_2": ["B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"],
        "recurrence_count": [95, 85, 75, 65, 55, 45, 35, 25, 18, 12, 8, 3],
    })
    
    result = assess_distribution_generalization(fusion_df, cosmic_df)
    
    assert "high_recurrence_rho" in result
    assert "mid_recurrence_rho" in result
    assert "low_recurrence_rho" in result
    assert "stability_index" in result
    assert "n_high" in result
    assert "n_mid" in result
    assert "n_low" in result
    
    # Check stability index if computed
    if result["stability_index"] is not None:
        assert result["stability_index"] >= 0.0


def test_external_validity_insufficient_data():
    """Test external validity with insufficient data."""
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
    
    result = assess_distribution_generalization(fusion_df, cosmic_df)
    
    # Should return None or zero values
    assert result["stability_index"] is None or result["n_high"] == 0


# =============================================================================
# Test 3: COSMIC Bias Quantification
# =============================================================================

def test_bias_quantification_returns_metrics():
    """Test that bias quantification returns metrics."""
    # Create test COSMIC data with known bias (top-heavy)
    cosmic_df = pd.DataFrame({
        "gene_1": ["A"] * 100 + ["B"] * 50 + ["C"] * 30 + ["D"] * 20 + ["E"] * 10 + ["F"] * 5,
        "gene_2": ["X"] * 100 + ["Y"] * 50 + ["Z"] * 30 + ["W"] * 20 + ["V"] * 10 + ["U"] * 5,
        "recurrence_count": [100] * 100 + [50] * 50 + [30] * 30 + [20] * 20 + [10] * 10 + [5] * 5,
    })
    
    result = quantify_cosmic_sampling_bias(cosmic_df)
    
    assert "top5_gene_fraction" in result
    assert "gini_gene_recurrence" in result
    assert "tail_coverage_ratio" in result
    assert "bias_severity_classification" in result
    assert "n_unique_genes" in result
    assert "n_fusion_pairs" in result
    
    # Check values are valid
    if result["top5_gene_fraction"] is not None:
        assert 0.0 <= result["top5_gene_fraction"] <= 1.0
    
    if result["gini_gene_recurrence"] is not None:
        assert 0.0 <= result["gini_gene_recurrence"] <= 1.0
    
    assert result["bias_severity_classification"] in ["LOW", "MODERATE", "HIGH", "N/A"]


def test_bias_quantification_empty_data():
    """Test bias quantification with empty data."""
    cosmic_df = pd.DataFrame({
        "gene_1": [],
        "gene_2": [],
        "recurrence_count": [],
    })
    
    result = quantify_cosmic_sampling_bias(cosmic_df)
    
    assert result["top5_gene_fraction"] is None
    assert result["n_unique_genes"] == 0


# =============================================================================
# Test 4: Scientific Claim Strength Classifier
# =============================================================================

def test_claim_classifier_returns_classification():
    """Test that claim classifier returns classification."""
    # Strong evidence case
    result = classify_scientific_claim_strength(
        observed_rho=0.75,
        null_p=0.001,
        stability_index=0.10,
    )
    
    assert "classification" in result
    assert "confidence_level" in result
    assert "interpretation" in result
    
    assert result["classification"] in [
        "EXPLORATORY SIGNAL",
        "MODERATE EXTERNAL CONSISTENCY",
        "STRONG BIOLOGICAL AGREEMENT",
        "HIGH CONFIDENCE CROSS-DATASET SIGNAL",
    ]
    
    assert result["confidence_level"] in ["LOW", "MODERATE", "HIGH", "VERY_HIGH"]
    assert len(result["interpretation"]) > 0


def test_claim_classifier_weak_evidence():
    """Test claim classifier with weak evidence."""
    result = classify_scientific_claim_strength(
        observed_rho=0.2,
        null_p=0.5,
        stability_index=0.3,
    )
    
    assert result["classification"] == "EXPLORATORY SIGNAL"
    assert result["confidence_level"] == "LOW"


def test_claim_classifier_missing_values():
    """Test claim classifier with missing values."""
    result = classify_scientific_claim_strength(
        observed_rho=None,
        null_p=None,
        stability_index=None,
    )
    
    assert result["classification"] == "EXPLORATORY SIGNAL"
    assert result["confidence_level"] == "LOW"


# =============================================================================
# Test 5: Integration - No Crashes with Small Overlap
# =============================================================================

def test_no_crashes_small_overlap():
    """Test that all functions handle small overlap gracefully."""
    fusion_df = pd.DataFrame({
        "gene_1": ["A", "B", "C"],
        "gene_2": ["B", "C", "D"],
        "recurrence_count": [100, 80, 60],
    })
    
    cosmic_df = pd.DataFrame({
        "gene_1": ["A", "B", "X"],
        "gene_2": ["B", "C", "Y"],
        "recurrence_count": [90, 70, 50],
    })
    
    # All functions should complete without crashing
    null_result = run_cosmic_null_model_test(fusion_df, cosmic_df, n_simulations=10)
    assert null_result is not None
    
    validity_result = assess_distribution_generalization(fusion_df, cosmic_df)
    assert validity_result is not None
    
    bias_result = quantify_cosmic_sampling_bias(cosmic_df)
    assert bias_result is not None
    
    claim_result = classify_scientific_claim_strength(
        observed_rho=0.5,
        null_p=0.05,
        stability_index=0.2,
    )
    assert claim_result is not None


# =============================================================================
# Test 6: All Existing Tests Still Pass (Integration Check)
# =============================================================================

def test_existing_correlation_function_unchanged():
    """Test that existing correlation function is still accessible and unchanged."""
    from week2_validation.cosmic.diagnostics import compute_spearman_correlation
    
    fusion_pairs = {
        ("A", "B"): 100,
        ("C", "D"): 80,
        ("E", "F"): 60,
    }
    
    cosmic_pairs = {
        ("A", "B"): 90,
        ("C", "D"): 70,
        ("E", "F"): 50,
    }
    
    overlap = set(fusion_pairs.keys())
    
    result = compute_spearman_correlation(fusion_pairs, cosmic_pairs, overlap)
    
    assert "spearman_rho" in result
    assert "spearman_p_value" in result
