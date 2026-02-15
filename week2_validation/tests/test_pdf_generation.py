"""
Week 2: Data Integrity & Statistical Validation — PDF Generation Tests.

Tests PDF report generation functionality to ensure:
- PDF file is created successfully
- PDF file size > 0
- Pipeline does not crash when PDF disabled
- PDF generation works with minimal results JSON
"""

import json
import tempfile
from pathlib import Path

import pytest

try:
    from reportlab.lib.pagesizes import letter
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False


@pytest.mark.skipif(not REPORTLAB_AVAILABLE, reason="ReportLab not available")
def test_pdf_generation_basic():
    """Test that PDF is generated successfully with basic data."""
    from week2_validation.reporting.pdf_report_generator import generate_week2_pdf_report

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # Create minimal status JSON
        status_json = {
            "week2_version": "1.0.0",
            "run_timestamp_utc": "2025-02-03T12:00:00Z",
            "dataset_hash": "test_hash",
            "diagnostics_run": [],
            "exit_code": 0,
            "status": "SUCCESS",
            "notes": [],
            "data_quality": {
                "total_rows_original": 100,
                "rows_used_for_analysis": 95,
                "rows_excluded": 5,
                "excluded_fraction": 0.05,
                "warning_flag": False,
                "high_risk_flag": False,
            },
        }

        status_path = tmp_path / "week2_status_test.json"
        with open(status_path, "w", encoding="utf-8") as f:
            json.dump(status_json, f)

        pdf_path = tmp_path / "test_report.pdf"
        result = generate_week2_pdf_report(
            results_json_path=status_path,
            output_pdf_path=pdf_path,
            run_metadata={},
        )

        # Verify PDF was created
        assert result is not None
        assert pdf_path.exists()
        assert pdf_path.stat().st_size > 0


@pytest.mark.skipif(not REPORTLAB_AVAILABLE, reason="ReportLab not available")
def test_pdf_generation_with_diagnostics():
    """Test PDF generation with diagnostic results."""
    from week2_validation.reporting.pdf_report_generator import generate_week2_pdf_report

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # Create status JSON
        status_json = {
            "week2_version": "1.0.0",
            "run_timestamp_utc": "2025-02-03T12:00:00Z",
            "dataset_hash": "test_hash",
            "diagnostics_run": ["diagnostics", "benford", "cosmic"],
            "exit_code": 0,
            "status": "SUCCESS",
            "notes": [],
            "data_quality": {
                "total_rows_original": 1000,
                "rows_used_for_analysis": 950,
                "rows_excluded": 50,
                "excluded_fraction": 0.05,
                "warning_flag": False,
                "high_risk_flag": False,
            },
        }

        status_path = tmp_path / "week2_status_test.json"
        with open(status_path, "w", encoding="utf-8") as f:
            json.dump(status_json, f)

        # Create diagnostic results JSON
        diagnostic_json = {
            "distribution": {
                "minimum": 10.0,
                "maximum": 1000.0,
                "mean": 250.5,
                "median": 200.0,
                "std": 150.2,
                "skewness": 1.5,
            },
            "log_normality": {
                "ks_statistic": 0.05,
                "ad_statistic": 0.8,
                "scipy_available": True,
            },
            "benford": {
                "observed_frequencies": {1: 0.3, 2: 0.2, 3: 0.15},
                "expected_frequencies": {1: 0.301, 2: 0.176, 3: 0.125},
                "applicability": True,
                "scale_span_orders_of_magnitude": 2.5,
            },
            "cosmic": {
                "total_fusions_ours": 50,
                "total_fusions_cosmic": 100,
                "overlap_count": 20,
                "only_in_ours_count": 30,
                "only_in_cosmic_count": 80,
                "spearman_rho": 0.75,
                "spearman_p_value": 0.001,
                "enrichment_p_value": 0.002,
                "cosmic_validation_score": 0.85,
                "cosmic_validation_classification": "STRONG BIOLOGICAL AGREEMENT",
                "cosmic_reference_source": "cosmic_v103_grch38_real",
                "cosmic_reference_version": "v103",
            },
        }

        diagnostic_path = tmp_path / "week2_diagnostic_results_test_hash.json"
        with open(diagnostic_path, "w", encoding="utf-8") as f:
            json.dump(diagnostic_json, f)

        pdf_path = tmp_path / "test_report.pdf"
        result = generate_week2_pdf_report(
            results_json_path=status_path,
            output_pdf_path=pdf_path,
            run_metadata={},
        )

        # Verify PDF was created
        assert result is not None
        assert pdf_path.exists()
        assert pdf_path.stat().st_size > 0


@pytest.mark.skipif(not REPORTLAB_AVAILABLE, reason="ReportLab not available")
def test_pdf_generation_missing_diagnostics():
    """Test PDF generation when diagnostic results are missing."""
    from week2_validation.reporting.pdf_report_generator import generate_week2_pdf_report

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # Create minimal status JSON without diagnostic results
        status_json = {
            "week2_version": "1.0.0",
            "run_timestamp_utc": "2025-02-03T12:00:00Z",
            "dataset_hash": "test_hash",
            "diagnostics_run": [],
            "exit_code": 0,
            "status": "SUCCESS",
            "notes": [],
            "data_quality": {
                "total_rows_original": 100,
                "rows_used_for_analysis": 100,
                "rows_excluded": 0,
                "excluded_fraction": 0.0,
                "warning_flag": False,
                "high_risk_flag": False,
            },
        }

        status_path = tmp_path / "week2_status_test.json"
        with open(status_path, "w", encoding="utf-8") as f:
            json.dump(status_json, f)

        pdf_path = tmp_path / "test_report.pdf"
        result = generate_week2_pdf_report(
            results_json_path=status_path,
            output_pdf_path=pdf_path,
            run_metadata={},
        )

        # PDF should still be generated even without diagnostics
        assert result is not None
        assert pdf_path.exists()
        assert pdf_path.stat().st_size > 0


@pytest.mark.skipif(not REPORTLAB_AVAILABLE, reason="ReportLab not available")
def test_pdf_generation_invalid_json():
    """Test PDF generation handles invalid JSON gracefully."""
    from week2_validation.reporting.pdf_report_generator import generate_week2_pdf_report

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # Create invalid JSON file
        status_path = tmp_path / "week2_status_test.json"
        status_path.write_text("invalid json content", encoding="utf-8")

        pdf_path = tmp_path / "test_report.pdf"
        result = generate_week2_pdf_report(
            results_json_path=status_path,
            output_pdf_path=pdf_path,
            run_metadata={},
        )

        # Should return None on error, not crash
        assert result is None


@pytest.mark.skipif(not REPORTLAB_AVAILABLE, reason="ReportLab not available")
def test_pdf_generation_missing_file():
    """Test PDF generation handles missing status file gracefully."""
    from week2_validation.reporting.pdf_report_generator import generate_week2_pdf_report

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # Use non-existent file
        status_path = tmp_path / "nonexistent.json"
        pdf_path = tmp_path / "test_report.pdf"

        result = generate_week2_pdf_report(
            results_json_path=status_path,
            output_pdf_path=pdf_path,
            run_metadata={},
        )

        # Should return None on error, not crash
        assert result is None


def test_pdf_generation_without_reportlab():
    """Test that PDF generation gracefully handles missing ReportLab."""
    # This test runs even if ReportLab is not available
    try:
        from week2_validation.reporting.pdf_report_generator import generate_week2_pdf_report

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)

            status_json = {
                "week2_version": "1.0.0",
                "run_timestamp_utc": "2025-02-03T12:00:00Z",
                "dataset_hash": "test_hash",
                "diagnostics_run": [],
                "exit_code": 0,
                "status": "SUCCESS",
                "notes": [],
            }

            status_path = tmp_path / "week2_status_test.json"
            with open(status_path, "w", encoding="utf-8") as f:
                json.dump(status_json, f)

            pdf_path = tmp_path / "test_report.pdf"
            result = generate_week2_pdf_report(
                results_json_path=status_path,
                output_pdf_path=pdf_path,
                run_metadata={},
            )

            # If ReportLab is not available, should return None gracefully
            if not REPORTLAB_AVAILABLE:
                assert result is None
            else:
                # If ReportLab is available, should generate PDF
                assert result is not None
                assert pdf_path.exists()
    except ImportError:
        # If module can't be imported due to missing ReportLab, that's OK
        pass


@pytest.mark.skipif(not REPORTLAB_AVAILABLE, reason="ReportLab not available")
def test_pdf_generation_with_certification():
    """Test PDF generation with certification data."""
    from week2_validation.reporting.pdf_report_generator import generate_week2_pdf_report

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        status_json = {
            "week2_version": "1.0.0",
            "run_timestamp_utc": "2025-02-03T12:00:00Z",
            "dataset_hash": "test_hash",
            "diagnostics_run": ["diagnostics"],
            "exit_code": 0,
            "status": "SUCCESS",
            "notes": [],
            "data_quality": {
                "total_rows_original": 100,
                "rows_used_for_analysis": 100,
                "rows_excluded": 0,
                "excluded_fraction": 0.0,
                "warning_flag": False,
                "high_risk_flag": False,
            },
        }

        status_path = tmp_path / "week2_status_test.json"
        with open(status_path, "w", encoding="utf-8") as f:
            json.dump(status_json, f)

        # Create certification JSON
        cert_json = {
            "approved_for_powerlaw_modeling": True,
            "approval_timestamp": "2025-02-03T12:00:00Z",
        }

        cert_path = tmp_path / "week2_dataset_certification_test_hash.json"
        with open(cert_path, "w", encoding="utf-8") as f:
            json.dump(cert_json, f)

        pdf_path = tmp_path / "test_report.pdf"
        result = generate_week2_pdf_report(
            results_json_path=status_path,
            output_pdf_path=pdf_path,
            run_metadata={},
        )

        assert result is not None
        assert pdf_path.exists()
        assert pdf_path.stat().st_size > 0
