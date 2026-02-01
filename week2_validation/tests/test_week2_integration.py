"""
Week 2 Integration Tests - Stability Infrastructure.

Tests freeze determinism, schema validation contract, SciPy optional mode,
and JSON summary structure. Non-destructive, diagnostic-only.
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

# Schema path relative to week2_validation package
_SCHEMA_PATH = Path(__file__).resolve().parent.parent / "schemas" / "fusion_schema.yaml"


def _load_fusion_schema():
    """Load fusion_schema.yaml. Used only within tests."""
    import yaml
    with open(_SCHEMA_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _validate_dataframe_against_schema(df: pd.DataFrame, schema: dict) -> None:
    """
    Validate DataFrame against fusion_schema contract. Validation-only; no auto-fix.
    Raises ValueError if validation fails.
    """
    required = schema.get("required_columns", [])
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    for col in required:
        if df[col].isnull().any():
            raise ValueError(f"Nulls not allowed in required column: {col}")
    constraints = schema.get("constraints", {})
    if "protein_length" in df.columns:
        pl = constraints.get("protein_length", {})
        if pl.get("min_exclusive") is not None and (df["protein_length"] <= 0).any():
            raise ValueError("protein_length must be > 0")
    if "recurrence_count" in df.columns:
        rc = constraints.get("recurrence_count", {})
        if rc.get("min_inclusive") is not None and (df["recurrence_count"] < 0).any():
            raise ValueError("recurrence_count must be >= 0")


# -----------------------------------------------------------------------------
# Test 1 — Freeze Determinism
# -----------------------------------------------------------------------------
def test_freeze_determinism_same_input_twice_same_frozen_directory():
    """Same input twice must produce the same frozen directory path."""
    from week2_validation.utils.freeze import ensure_frozen_input

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        fusion_file = root / "fusion.csv"
        fusion_file.write_text("gene_1,gene_2,protein_length,recurrence_count\nA,B,100,1\n", encoding="utf-8")
        frozen_root = root / "frozen_inputs"
        frozen_root.mkdir(parents=True, exist_ok=True)

        path1 = ensure_frozen_input(fusion_file, frozen_root)
        path2 = ensure_frozen_input(fusion_file, frozen_root)

        assert path1 == path2
        assert path1.is_dir()
        assert (path1 / "manifest.json").exists()
        assert (path1 / ".frozen").exists()


# -----------------------------------------------------------------------------
# Test 2 — Schema Validation Failure
# -----------------------------------------------------------------------------
def test_schema_validation_failure_missing_required_column():
    """Missing required column must cause validation to fail."""
    schema = _load_fusion_schema()
    # DataFrame missing gene_1 (and gene_2) — only protein_length, recurrence_count
    df = pd.DataFrame({
        "protein_length": [100],
        "recurrence_count": [1],
    })
    with pytest.raises(ValueError, match="Missing required column"):
        _validate_dataframe_against_schema(df, schema)


# -----------------------------------------------------------------------------
# Test 3 — SciPy Optional Mode
# -----------------------------------------------------------------------------
def test_scipy_optional_mode_diagnostics_still_run():
    """Simulate SciPy missing; diagnostics must still run (no crash)."""
    script = """
import sys
import types
# Simulate SciPy missing: mock that raises ImportError on attribute access
def _raise_import(name):
    raise ImportError("scipy not available")
scipy_mock = types.ModuleType("scipy")
scipy_mock.__getattr__ = _raise_import
sys.modules["scipy"] = scipy_mock
# Import after blocking scipy so module sees scipy as unavailable
from week2_validation.benford.diagnostics import run_benford_diagnostics
import numpy as np
data = np.array([1.0, 2.0, 3.0] * 20)
result = run_benford_diagnostics(data, skip_applicability_check=True)
assert result is not None
# Diagnostic ran without crash; scipy optional mode exercised
print("ok")
"""
    proc = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=str(Path(__file__).resolve().parent.parent.parent),
        env={**__import__("os").environ, "PYTHONPATH": str(Path(__file__).resolve().parent.parent.parent)},
    )
    assert proc.returncode == 0, f"stderr: {proc.stderr!r} stdout: {proc.stdout!r}"
    assert "ok" in proc.stdout


# -----------------------------------------------------------------------------
# Test 4 — JSON Summary Structure
# -----------------------------------------------------------------------------
def test_json_summary_structure_contains_required_keys():
    """Summary JSON must contain required keys."""
    from week2_validation.reporting.json_summary import write_week2_summary_json

    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        run_metadata = {
            "week2_version": "0.1.0",
            "run_timestamp_utc": "2025-02-01T12:00:00Z",
            "dataset_hash": "abc123",
        }
        diagnostic_results = {
            "diagnostics_run": ["benford"],
            "results": {},
            "notes": [],
        }
        out_path = write_week2_summary_json(output_dir, run_metadata, diagnostic_results)
        assert out_path.exists()

        with open(out_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        required_keys = {"week2_version", "run_timestamp_utc", "dataset_hash", "diagnostics_run", "results", "notes"}
        for key in required_keys:
            assert key in data, f"Missing required key: {key}"
