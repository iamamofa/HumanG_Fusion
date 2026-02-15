"""
Data Integrity & Statistical Validation — Approved Dataset Writer (Deliverable Completion Layer).

When validation succeeds, optionally writes:
  1. validation_cleaned_dataset_{stem}.csv — copy of frozen validated dataset (no modifications).
  2. validation_dataset_certification_{stem}.json — certification record for power-law modeling.

Additive only. Does not affect pipeline success/failure, freeze, or diagnostics.
"""

import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Required columns per fusion_schema.yaml (contract only; no import from schemas)
_REQUIRED_COLUMNS = frozenset({
    "fusion_id",
    "gene_1",
    "gene_2",
    "protein_length",
    "recurrence_count",
})

from week2_validation.output_names import (
    CLEANED_DATASET_FILENAME_PATTERN,
    CERTIFICATION_FILENAME_PATTERN,
)
CERTIFICATION_VERSION = "1.0"


def _resolve_output_path(output_dir: Path, filename: str) -> Path:
    """
    Resolve output path under output_dir. Ensures no path traversal.
    """
    output_dir = output_dir.resolve()
    candidate = (output_dir / filename).resolve()
    try:
        candidate.relative_to(output_dir)
    except ValueError:
        return None  # path escaped output_dir
    return candidate


def _read_frozen_dataset(path: Path):
    """
    Read frozen dataset in any supported format (CSV, TSV, JSON, Parquet, Excel).
    Uses data_loader.load_data for format-agnostic loading.
    Does not modify values, drop rows, or change dtypes (preserve as read).
    """
    from week2_validation.utils.data_loader import load_data
    return load_data(str(path))


def write_approved_dataset(
    frozen_dataset_path: Path,
    output_dir: Path,
    dataset_hash: str,
    validation_passed: bool,
    diagnostics_run: List[str],
    runtime_metadata: Dict[str, Any],
    dataset_stem: Optional[str] = None,
) -> Tuple[Path | None, Path | None]:
    """
    When validation_passed is True, write cleaned dataset snapshot and certification record.

    Behavior:
      - If validation_passed is False: return (None, None); no writes.
      - If True: write validation_cleaned_dataset_{stem}.csv (copy of frozen data, no modifications)
                and validation_dataset_certification_{stem}.json (exact structure as specified).

    Processing rules for CSV:
      - Do NOT modify values, drop rows, recalculate, reorder columns, or change dtypes.
      - ONLY ensure required columns exist; write CSV safely; use atomic write (temp + replace).

    Security: atomic writes, JSON safe dump, path resolve + root check, no path traversal,
    no user input injection, no stack traces in JSON.

    Args:
        frozen_dataset_path: Path to frozen validated dataset file.
        output_dir: Directory for output artifacts.
        dataset_hash: Hash identifying the frozen dataset (e.g. frozen dir name).
        validation_passed: If True, write artifacts; if False, return (None, None).
        diagnostics_run: List of diagnostic names that were run.
        runtime_metadata: Dict with scipy_available, psutil_available (bool).

    Returns:
        (path_to_csv, path_to_json) or (None, None) on skip/failure.
    """
    if not validation_passed:
        return (None, None)

    output_dir = Path(output_dir).resolve()
    if not output_dir.is_dir():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            return (None, None)

    frozen_dataset_path = Path(frozen_dataset_path).resolve()
    if not frozen_dataset_path.is_file():
        return (None, None)

    stem = dataset_stem if dataset_stem else "default"
    cleaned_name = CLEANED_DATASET_FILENAME_PATTERN.format(stem=stem)
    cert_name = CERTIFICATION_FILENAME_PATTERN.format(stem=stem)
    csv_path = _resolve_output_path(output_dir, cleaned_name)
    cert_path = _resolve_output_path(output_dir, cert_name)
    if csv_path is None or cert_path is None:
        return (None, None)

    try:
        df = _read_frozen_dataset(frozen_dataset_path)
    except Exception:
        return (None, None)

    required = _REQUIRED_COLUMNS
    if not required.issubset(set(df.columns)):
        return (None, None)

    # Write CSV atomically: temp in output_dir then replace
    try:
        fd, tmp_path = tempfile.mkstemp(
            suffix=".csv",
            prefix=".validation_cleaned.",
            dir=str(output_dir),
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="") as f:
                df.to_csv(f, sep=",", index=False, date_format=None)
            os.replace(tmp_path, str(csv_path))
        except OSError:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
            return (None, None)
    except OSError:
        return (None, None)

    # Build certification payload (exact structure; no stack traces, no path traversal)
    frozen_source_display = frozen_dataset_path.name

    payload: Dict[str, Any] = {
        "certification_version": CERTIFICATION_VERSION,
        "approved_for_powerlaw_modeling": True,
        "dataset_hash": str(dataset_hash),
        "frozen_dataset_source": frozen_source_display,
        "validation_layer": "week2_validation",
        "validation_passed": True,
        "diagnostics_run": list(diagnostics_run),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_environment": {
            "scipy_available": bool(runtime_metadata.get("scipy_available", False)),
            "psutil_available": bool(runtime_metadata.get("psutil_available", False)),
        },
    }

    # Write certification JSON atomically
    try:
        fd, tmp_path = tempfile.mkstemp(
            suffix=".json",
            prefix=".validation_cert.",
            dir=str(output_dir),
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
            os.replace(tmp_path, str(cert_path))
        except (OSError, TypeError, ValueError):
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
            return (csv_path, None)
    except OSError:
        return (csv_path, None)

    return (csv_path, cert_path)
