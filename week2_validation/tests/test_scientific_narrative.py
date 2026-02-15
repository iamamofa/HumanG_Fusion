"""
Week 2: Data Integrity & Statistical Validation — Scientific Narrative Tests.

Tests scientific narrative generation functionality to ensure:
- Returns string
- Handles missing fields gracefully
- No exceptions raised
- No numeric mutation
"""

import json
import tempfile
from pathlib import Path

import pytest


def test_narrative_generation_basic():
    """Test that narrative is generated successfully with basic metrics."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_integrated_scientific_narrative

    metrics = {
        "skewness": 0.3,
        "mean": 500.0,
        "median": 510.0,
        "spearman_rho": 0.75,
        "cosmic_overlap": 20,
    }

    result = generate_integrated_scientific_narrative(metrics)

    # Should return a string
    assert isinstance(result, str)
    assert len(result) > 0


def test_narrative_generation_missing_fields():
    """Test that narrative handles missing fields gracefully."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_integrated_scientific_narrative

    # Empty metrics
    result = generate_integrated_scientific_narrative({})

    assert isinstance(result, str)
    assert len(result) > 0
    # Should return default message
    assert "requires additional diagnostic metrics" in result.lower() or "interpretation" in result.lower()


def test_narrative_generation_partial_metrics():
    """Test narrative generation with partial metrics."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_integrated_scientific_narrative

    # Only skewness available
    metrics = {
        "skewness": 0.5,
    }

    result = generate_integrated_scientific_narrative(metrics)

    assert isinstance(result, str)
    assert len(result) > 0


def test_narrative_no_exceptions():
    """Test that narrative generation never raises exceptions."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_integrated_scientific_narrative

    # Test with various invalid inputs
    test_cases = [
        {},  # Empty
        {"skewness": None},  # None values
        {"skewness": "invalid"},  # Invalid types
        {"mean": float("inf")},  # Infinity
        {"mean": float("nan")},  # NaN
    ]

    for metrics in test_cases:
        result = generate_integrated_scientific_narrative(metrics)
        assert isinstance(result, str)
        assert len(result) > 0


def test_narrative_no_numeric_mutation():
    """Test that narrative generation does not mutate numeric values."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_integrated_scientific_narrative

    original_metrics = {
        "skewness": 0.3,
        "mean": 500.0,
        "median": 510.0,
        "spearman_rho": 0.75,
    }

    metrics_copy = original_metrics.copy()
    result = generate_integrated_scientific_narrative(metrics_copy)

    # Verify original values unchanged
    assert metrics_copy["skewness"] == 0.3
    assert metrics_copy["mean"] == 500.0
    assert metrics_copy["median"] == 510.0
    assert metrics_copy["spearman_rho"] == 0.75


def test_confidence_statement_generation():
    """Test confidence statement generation."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_confidence_statement

    metrics = {
        "cosmic_overlap": 5,
        "spearman_p_value": 0.1,
    }

    result = generate_confidence_statement(metrics)

    assert isinstance(result, str)
    assert len(result) > 0


def test_confidence_statement_missing_fields():
    """Test confidence statement with missing fields."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_confidence_statement

    result = generate_confidence_statement({})

    assert isinstance(result, str)
    assert len(result) > 0


def test_limitations_section_generation():
    """Test limitations section generation."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_limitations_section

    result = generate_limitations_section({})

    assert isinstance(result, str)
    assert len(result) > 0
    # Should contain key limitation topics
    assert "statistical" in result.lower() or "limitations" in result.lower()
    assert "cosmic" in result.lower()
    assert "biological" in result.lower()


def test_extract_narrative_metrics():
    """Test metric extraction from status and diagnostic data."""
    from week2_validation.reporting.narrative.scientific_narrative import extract_narrative_metrics

    status_data = {
        "data_quality": {
            "skewness": 0.3,
        },
    }

    diagnostic_data = {
        "distribution": {
            "mean": 500.0,
            "median": 510.0,
            "skewness": 0.3,
        },
        "cosmic": {
            "overlap_count": 20,
            "spearman_rho": 0.75,
            "spearman_p_value": 0.01,
            "cosmic_validation_score": 0.8,
        },
    }

    metrics = extract_narrative_metrics(status_data, diagnostic_data)

    assert isinstance(metrics, dict)
    assert "skewness" in metrics
    assert "mean" in metrics
    assert "cosmic_overlap" in metrics
    assert "spearman_rho" in metrics


def test_narrative_scientific_tone():
    """Test that narrative uses appropriate scientific tone."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_integrated_scientific_narrative

    metrics = {
        "skewness": 0.3,
        "spearman_rho": 0.75,
    }

    result = generate_integrated_scientific_narrative(metrics)

    # Should use scientific language (not marketing)
    assert "consistent with" in result.lower() or "suggests" in result.lower() or "indicates" in result.lower()
    
    # Should NOT use certainty language
    assert "proves" not in result.lower()
    assert "guarantees" not in result.lower()
    assert "confirms biological truth" not in result.lower()


def test_confidence_statement_caution():
    """Test that confidence statement includes caution when appropriate."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_confidence_statement

    # Low overlap case
    metrics_low_overlap = {
        "cosmic_overlap": 5,
    }

    result = generate_confidence_statement(metrics_low_overlap)

    assert isinstance(result, str)
    assert "caution" in result.lower() or "instability" in result.lower() or "limited" in result.lower()


def test_narrative_integrated_story():
    """Test that narrative generates integrated story, not metric-by-metric."""
    from week2_validation.reporting.narrative.scientific_narrative import generate_integrated_scientific_narrative

    metrics = {
        "skewness": 0.3,
        "mean": 500.0,
        "median": 510.0,
        "spearman_rho": 0.75,
        "cosmic_overlap": 20,
        "cosmic_quality_score": 0.8,
    }

    result = generate_integrated_scientific_narrative(metrics)

    # Should be integrated paragraph, not bullet points
    assert isinstance(result, str)
    assert len(result) > 50  # Should be substantial
    # Should flow as paragraph (multiple sentences)
    assert result.count(".") >= 2
