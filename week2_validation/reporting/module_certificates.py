"""
Data Integrity & Statistical Validation — Per-Module Validation Certificates.

Generates individual validation certificates for each diagnostic module under
results/{output}/certificates/, building up to the final aggregate certificate
(cleaned, validated dataset approved for power-law modeling).

Certificate order:
  01_distribution_validation.json
  02_benford_validation.json
  03_lognormality_validation.json
  04_cosmic_validation.json
  05_aggregate_validation.json
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

try:
    from week2_validation import __version__ as MODULE_VERSION
except ImportError:
    MODULE_VERSION = "unknown"


CERTIFICATES_SUBDIR = "certificates"

# Filenames (order preserved for report build-up)
CERT_DISTRIBUTION = "01_distribution_validation.json"
CERT_BENFORD = "02_benford_validation.json"
CERT_LOGNORMALITY = "03_lognormality_validation.json"
CERT_COSMIC = "04_cosmic_validation.json"
CERT_AGGREGATE = "05_aggregate_validation.json"


def _certificates_dir(output_dir: Path) -> Path:
    """Return path to certificates subdir; create if needed."""
    cert_dir = Path(output_dir) / CERTIFICATES_SUBDIR
    cert_dir.mkdir(parents=True, exist_ok=True)
    return cert_dir


def _timestamp_iso8601() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _provenance(dataset_hash: str, input_file_hash: Optional[str] = None) -> Dict[str, Any]:
    """Build provenance block. Frozen file hash = dataset_hash (freeze dir name)."""
    return {
        "frozen_file_hash": dataset_hash,
        "input_file_hash": input_file_hash if input_file_hash is not None else dataset_hash,
    }


def _gate_status_for_distribution(
    diagnostic_results: Dict[str, Any], data_quality: Optional[Dict[str, Any]]
) -> str:
    from week2_validation.reporting.quality_gates import (
        GateStatus,
        evaluate_distribution_gate,
    )
    result = evaluate_distribution_gate(diagnostic_results, data_quality)
    return result.status.value


def _gate_status_for_benford(diagnostic_results: Dict[str, Any]) -> str:
    from week2_validation.reporting.quality_gates import (
        GateStatus,
        evaluate_benford_gate,
    )
    result = evaluate_benford_gate(diagnostic_results)
    return result.status.value


def _gate_status_for_lognormality(diagnostic_results: Dict[str, Any]) -> str:
    from week2_validation.reporting.quality_gates import (
        GateStatus,
        evaluate_log_normality_gate,
    )
    result = evaluate_log_normality_gate(diagnostic_results)
    return result.status.value


def _gate_status_for_cosmic(diagnostic_results: Dict[str, Any]) -> str:
    from week2_validation.reporting.quality_gates import (
        GateStatus,
        evaluate_cosmic_gate,
    )
    result = evaluate_cosmic_gate(diagnostic_results)
    return result.status.value


def _serialize(value: Any) -> Any:
    """Convert numpy/types for JSON."""
    if value is None:
        return None
    try:
        import numpy as np
        if isinstance(value, (np.integer, np.int64, np.int32)):
            return int(value)
        if isinstance(value, (np.floating, np.float64, np.float32)):
            return float(value)
        if isinstance(value, np.bool_):
            return bool(value)
    except ImportError:
        pass
    if isinstance(value, (int, float, str, bool)):
        return value
    if isinstance(value, dict):
        return {k: _serialize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_serialize(v) for v in value]
    return value


# -----------------------------------------------------------------------------
# Per-module certificate writers
# -----------------------------------------------------------------------------


def write_distribution_certificate(
    output_dir: Path,
    run_metadata: Dict[str, Any],
    histogram_path_relative: str = "protein_distribution.png",
) -> Optional[Path]:
    """
    Write 01_distribution_validation.json after distribution diagnostics run.
    Contains: N, min, max, mean, median, std, skewness, histogram_path, gate_status.
    """
    cert_dir = _certificates_dir(output_dir)
    dataset_hash = run_metadata.get("dataset_hash", "")
    diag = run_metadata.get("diagnostic_results") or {}
    dq = run_metadata.get("data_quality") or {}
    dist = diag.get("distribution") or {}

    n = dq.get("total_rows_original") or dq.get("rows_used_for_analysis") or dist.get("n")
    if n is None and dist:
        n = None  # Keep None if not present

    metrics = {
        "N": _serialize(n),
        "min": dist.get("minimum"),
        "max": dist.get("maximum"),
        "mean": dist.get("mean"),
        "median": dist.get("median"),
        "std": dist.get("std"),
        "skewness": dq.get("skewness") or dist.get("skewness"),
        "histogram_path": histogram_path_relative,
    }
    metrics = _serialize(metrics)

    try:
        gate_status = _gate_status_for_distribution(diag, dq)
    except Exception:
        gate_status = "SKIPPED"

    payload = {
        "module_name": "distribution",
        "module_version": MODULE_VERSION,
        "timestamp": _timestamp_iso8601(),
        "dataset_hash": dataset_hash,
        "gate_status": gate_status,
        "metrics": metrics,
        "provenance": _provenance(dataset_hash, run_metadata.get("input_file_hash")),
    }
    path = cert_dir / CERT_DISTRIBUTION
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return path
    except OSError:
        return None


def write_benford_certificate(
    output_dir: Path,
    run_metadata: Dict[str, Any],
    control_results: Optional[Dict[str, Any]] = None,
) -> Optional[Path]:
    """
    Write 02_benford_validation.json after Benford diagnostics run.
    Contains: observed_freq, expected_freq, chi_squared, p_value, scale_span,
    applicability, gate_status, control_results (positive + negative).
    """
    cert_dir = _certificates_dir(output_dir)
    dataset_hash = run_metadata.get("dataset_hash", "")
    diag = run_metadata.get("diagnostic_results") or {}
    benford = diag.get("benford") or {}

    metrics = {
        "observed_freq": benford.get("observed_frequencies"),
        "expected_freq": benford.get("expected_frequencies"),
        "chi_squared": benford.get("chi_squared_statistic"),
        "p_value": benford.get("p_value"),
        "scale_span": benford.get("scale_span_orders_of_magnitude"),
        "applicability": benford.get("applicability"),
    }
    metrics = _serialize(metrics)
    if control_results is not None:
        metrics["control_results"] = _serialize(control_results)
    else:
        metrics["control_results"] = None

    try:
        gate_status = _gate_status_for_benford(diag)
    except Exception:
        gate_status = "SKIPPED"

    payload = {
        "module_name": "benford",
        "module_version": MODULE_VERSION,
        "timestamp": _timestamp_iso8601(),
        "dataset_hash": dataset_hash,
        "gate_status": gate_status,
        "metrics": metrics,
        "provenance": _provenance(dataset_hash, run_metadata.get("input_file_hash")),
    }
    path = cert_dir / CERT_BENFORD
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return path
    except OSError:
        return None


def write_lognormality_certificate(
    output_dir: Path,
    run_metadata: Dict[str, Any],
) -> Optional[Path]:
    """
    Write 03_lognormality_validation.json after log-normality diagnostics run.
    Contains: ks_statistic, ad_statistic, scipy_available, gate_status.
    """
    cert_dir = _certificates_dir(output_dir)
    dataset_hash = run_metadata.get("dataset_hash", "")
    diag = run_metadata.get("diagnostic_results") or {}
    logn = diag.get("log_normality") or {}

    metrics = {
        "ks_statistic": logn.get("ks_statistic"),
        "ad_statistic": logn.get("ad_statistic"),
        "scipy_available": logn.get("scipy_available"),
    }
    metrics = _serialize(metrics)

    try:
        gate_status = _gate_status_for_lognormality(diag)
    except Exception:
        gate_status = "SKIPPED"

    payload = {
        "module_name": "log_normality",
        "module_version": MODULE_VERSION,
        "timestamp": _timestamp_iso8601(),
        "dataset_hash": dataset_hash,
        "gate_status": gate_status,
        "metrics": metrics,
        "provenance": _provenance(dataset_hash, run_metadata.get("input_file_hash")),
    }
    path = cert_dir / CERT_LOGNORMALITY
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return path
    except OSError:
        return None


def write_cosmic_certificate(
    output_dir: Path,
    run_metadata: Dict[str, Any],
) -> Optional[Path]:
    """
    Write 04_cosmic_validation.json after COSMIC diagnostics run.
    Contains: spearman_rho, spearman_p, overlap_count, enrichment_p,
    validation_score, classification, gate_status.
    """
    cert_dir = _certificates_dir(output_dir)
    dataset_hash = run_metadata.get("dataset_hash", "")
    diag = run_metadata.get("diagnostic_results") or {}
    cosmic = diag.get("cosmic") or {}

    metrics = {
        "spearman_rho": cosmic.get("spearman_rho"),
        "spearman_p": cosmic.get("spearman_p_value"),
        "overlap_count": cosmic.get("overlap_count"),
        "enrichment_p": cosmic.get("enrichment_p_value"),
        "validation_score": cosmic.get("cosmic_validation_score"),
        "classification": cosmic.get("cosmic_validation_classification"),
    }
    metrics = _serialize(metrics)

    try:
        gate_status = _gate_status_for_cosmic(diag)
    except Exception:
        gate_status = "SKIPPED"

    payload = {
        "module_name": "cosmic",
        "module_version": MODULE_VERSION,
        "timestamp": _timestamp_iso8601(),
        "dataset_hash": dataset_hash,
        "gate_status": gate_status,
        "metrics": metrics,
        "provenance": _provenance(dataset_hash, run_metadata.get("input_file_hash")),
    }
    path = cert_dir / CERT_COSMIC
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return path
    except OSError:
        return None


def write_aggregate_certificate(
    output_dir: Path,
    run_metadata: Dict[str, Any],
) -> Optional[Path]:
    """
    Write 05_aggregate_validation.json after all modules complete.
    Contains: all individual gate results, overall approval status, dataset_hash,
    timestamp. This is the "cleaned, validated dataset approved for power-law
    modeling" deliverable summary.
    """
    cert_dir = _certificates_dir(output_dir)
    dataset_hash = run_metadata.get("dataset_hash", "")
    quality_gates = run_metadata.get("quality_gates") or {}
    gates = quality_gates.get("gates", [])
    overall_approval = quality_gates.get("overall_approval", "N/A")

    # Load individual certificates if present (for full build-up record)
    individual = {}
    for name, filename in [
        ("distribution", CERT_DISTRIBUTION),
        ("benford", CERT_BENFORD),
        ("log_normality", CERT_LOGNORMALITY),
        ("cosmic", CERT_COSMIC),
    ]:
        p = cert_dir / filename
        if p.is_file():
            try:
                with open(p, "r", encoding="utf-8") as f:
                    individual[name] = json.load(f)
            except (OSError, json.JSONDecodeError):
                individual[name] = None
        else:
            individual[name] = None

    payload = {
        "module_name": "aggregate",
        "module_version": MODULE_VERSION,
        "timestamp": _timestamp_iso8601(),
        "dataset_hash": dataset_hash,
        "gate_status": overall_approval,
        "overall_approval": overall_approval,
        "gate_results": [g if isinstance(g, dict) else g for g in gates],
        "individual_certificates": individual,
        "provenance": _provenance(dataset_hash, run_metadata.get("input_file_hash")),
    }
    path = cert_dir / CERT_AGGREGATE
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return path
    except OSError:
        return None
