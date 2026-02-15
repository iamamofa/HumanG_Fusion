"""
Data Integrity & Statistical Validation — Machine Output Layer (Structured JSON Summary).

This module writes a single structured JSON file (validation_integrity_summary.json)
when called. It is passive: it does not pull data from the pipeline, import
diagnostics, or modify state. It only writes JSON when explicitly invoked.
"""

import json
from pathlib import Path
from typing import Any, Dict

from week2_validation.output_names import VALIDATION_VERSION_KEY

# Output filename (contract)
SUMMARY_FILENAME = "validation_integrity_summary.json"


def write_week2_summary_json(
    output_dir: Path,
    run_metadata: dict,
    diagnostic_results: dict,
) -> Path:
    """
    Write Data Integrity & Statistical Validation integrity summary as structured JSON.

    Creates output_dir/validation_integrity_summary.json with required keys:
    validation_version, run_timestamp_utc, dataset_hash, diagnostics_run,
    results, notes.

    This function is passive: it does not fetch data, import diagnostics,
    or modify pipeline state. It only writes the provided dicts to JSON.

    Args:
        output_dir: Directory where the JSON file will be written.
        run_metadata: Dict with validation_version, run_timestamp_utc, dataset_hash.
        diagnostic_results: Dict with diagnostics_run, results; notes optional.

    Returns:
        Path to the written JSON file.

    Raises:
        OSError: If the file cannot be written.
    """
    payload: Dict[str, Any] = {
        VALIDATION_VERSION_KEY: run_metadata.get(VALIDATION_VERSION_KEY, run_metadata.get("week2_version", "")),
        "run_timestamp_utc": run_metadata.get("run_timestamp_utc", ""),
        "dataset_hash": run_metadata.get("dataset_hash", ""),
        "diagnostics_run": diagnostic_results.get("diagnostics_run", []),
        "results": diagnostic_results.get("results", {}),
        "notes": diagnostic_results.get("notes", run_metadata.get("notes", [])),
    }
    out_path = output_dir / SUMMARY_FILENAME
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    return out_path
