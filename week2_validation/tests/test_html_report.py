"""
Data Integrity & Statistical Validation — HTML Report Generation Tests.

Tests HTML report generation to ensure:
- HTML file is created successfully
- HTML contains expected content
- Same data as PDF report
"""

import json
import tempfile
from pathlib import Path

import pytest


def test_html_generation_basic():
    """Test that HTML is generated successfully with basic data."""
    from week2_validation.reporting.html_report_generator import generate_week2_html_report

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        status_json = {
            "week2_version": "1.0.0",
            "run_timestamp_utc": "2025-02-03T12:00:00Z",
            "dataset_hash": "test_hash",
            "exit_code": 0,
            "status": "SUCCESS",
            "data_quality": {
                "total_rows_original": 100,
                "rows_used_for_analysis": 95,
                "excluded_fraction": 0.05,
                "warning_flag": False,
                "high_risk_flag": False,
            },
        }

        status_path = tmp_path / "week2_status_test.json"
        with open(status_path, "w", encoding="utf-8") as f:
            json.dump(status_json, f)

        html_path = tmp_path / "test_report.html"
        result = generate_week2_html_report(
            results_json_path=status_path,
            output_html_path=html_path,
            run_metadata={},
        )

        assert result is not None
        assert html_path.exists()
        assert html_path.stat().st_size > 0

        content = html_path.read_text(encoding="utf-8")
        assert "Data Integrity" in content
        assert "General Statistics" in content
        assert "mqc-" in content or "Data Integrity" in content


def test_html_generation_with_diagnostics():
    """Test HTML generation with diagnostic results (same structure as PDF)."""
    from week2_validation.reporting.html_report_generator import generate_week2_html_report

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        status_json = {
            "run_timestamp_utc": "2025-02-03T12:00:00Z",
            "dataset_hash": "test_hash",
            "exit_code": 0,
            "status": "SUCCESS",
            "data_quality": {
                "total_rows_original": 100,
                "rows_used_for_analysis": 100,
                "excluded_fraction": 0,
                "warning_flag": False,
                "high_risk_flag": False,
            },
        }

        diagnostic_json = {
            "distribution": {"mean": 500, "median": 510, "std": 100, "skewness": 0.2},
            "cosmic": {
                "spearman_rho": 0.75,
                "overlap_count": 10,
                "cosmic_validation_classification": "MODERATE AGREEMENT",
                "cosmic_reference_version": "v103",
            },
        }

        status_path = tmp_path / "week2_status_test.json"
        diagnostic_path = tmp_path / "week2_diagnostic_results_test.json"
        with open(status_path, "w", encoding="utf-8") as f:
            json.dump(status_json, f)
        with open(diagnostic_path, "w", encoding="utf-8") as f:
            json.dump(diagnostic_json, f)

        html_path = tmp_path / "report.html"
        result = generate_week2_html_report(
            results_json_path=status_path,
            output_html_path=html_path,
            run_metadata={},
        )

        assert result is not None
        content = html_path.read_text(encoding="utf-8")
        assert "COSMIC" in content
        assert "0.750" in content or "0.75" in content
        assert "Distribution" in content
