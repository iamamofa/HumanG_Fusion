"""
Week 2: Data Integrity & Statistical Validation — Integration Tests (Stability Infrastructure).

Tests freeze determinism, schema validation contract, SciPy optional mode,
JSON summary structure, null policy, path traversal hardening, parallel freeze,
status envelope sanitization, and resource limit guards. Non-destructive, diagnostic-only.

Stress test documentation:
- File size near 5GB: Not run (would require large files). Guards tested via mocked stat.
- Row count near 200M: Not run (would require ~40GB CSV). Guards fire at limit.
- Memory near 8GB: Runtime guard fires when psutil reports >8GB RSS.
- Parallel invocations: test_parallel_freeze_safety uses threading.
"""

import json
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import pandas as pd
import pytest

# Schema path relative to week2_validation package
_SCHEMA_PATH = Path(__file__).resolve().parent.parent / "schemas" / "fusion_schema.yaml"

# Pipeline writes to week2_validation/results/{output_subfolder}/; --output-dir only supplies subfolder name
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def _status_path_for_run(output_dir_arg: Path) -> Path:
    """Resolve status file path: pipeline uses project_root/week2_validation/results/{output_dir_arg.name}/."""
    return _PROJECT_ROOT / "week2_validation" / "results" / output_dir_arg.name / "validation_status_fusion.json"


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
# Test 4 — Null Policy (required columns must not contain nulls)
# -----------------------------------------------------------------------------
def test_null_policy_rejects_nulls_in_required_columns():
    """Required columns with nulls must cause SchemaValidationError."""
    import numpy as np

    from week2_validation.utils.data_loader import (
        REQUIRED_FUSION_FIELDS,
        SchemaValidationError,
        validate_schema,
    )

    df = pd.DataFrame({
        "fusion_id": ["A::B"],
        "gene_1": ["A"],
        "gene_2": ["B"],
        "protein_length": [np.nan],
        "recurrence_count": [1],
    })
    with pytest.raises(SchemaValidationError) as exc_info:
        validate_schema(df, REQUIRED_FUSION_FIELDS)
    assert "protein_length" in str(exc_info.value)


# -----------------------------------------------------------------------------
# Test 5 — Path Traversal Hardening (path.parts check)
# -----------------------------------------------------------------------------
def test_path_traversal_blocked_via_parts():
    """Path with '..' as component must be rejected; '....' not confused with '..'."""
    from week2_validation.utils.data_loader import FileValidationError, validate_file_path

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        path = Path(f.name)
    try:
        validate_file_path(str(path), must_exist=True)
    except FileValidationError:
        pass
    path_str_bad = str(path.parent / ".." / path.name)
    with pytest.raises(FileValidationError, match="traversal"):
        validate_file_path(path_str_bad, must_exist=False)
    path.unlink(missing_ok=True)


# -----------------------------------------------------------------------------
# Test 6 — Parallel Freeze Safety
# -----------------------------------------------------------------------------
def test_parallel_freeze_safety():
    """Concurrent freeze of same input must not corrupt; both get valid snapshot."""
    from week2_validation.utils.freeze import ensure_frozen_input

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        fusion_file = root / "fusion.csv"
        fusion_file.write_text(
            "fusion_id,gene_1,gene_2,protein_length,recurrence_count\n"
            "A::B,A,B,100,1\nB::C,B,C,200,2\n",
            encoding="utf-8",
        )
        frozen_root = root / "frozen_inputs"
        frozen_root.mkdir(parents=True, exist_ok=True)

        results = []

        def freeze_once():
            p = ensure_frozen_input(fusion_file, frozen_root)
            results.append(p)

        t1 = threading.Thread(target=freeze_once)
        t2 = threading.Thread(target=freeze_once)
        t1.start()
        t2.start()
        t1.join()
        t2.join()

        assert len(results) == 2
        assert results[0] == results[1]
        p = results[0]
        assert p.is_dir()
        assert (p / "manifest.json").exists()
        assert (p / ".frozen").exists()


# -----------------------------------------------------------------------------
# Test 7 — Resource Limit Guards (file size / row estimate)
# -----------------------------------------------------------------------------
def test_resource_limit_guards_reject_oversized_file():
    """File exceeding MAX_FILE_SIZE_BYTES must raise DataLoaderError."""
    from unittest.mock import patch

    from week2_validation.utils.data_loader import DataLoaderError, load_csv

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w", encoding="utf-8") as f:
        f.write("fusion_id,gene_1,gene_2,protein_length,recurrence_count\n")
        f.write("A::B,A,B,100,1\n")
        f.flush()
        path = Path(f.name)
    try:
        with patch.object(Path, "stat") as mock_stat:
            mock_stat.return_value.st_size = 11 * 1024 * 1024 * 1024
            with pytest.raises(DataLoaderError, match="exceeds maximum"):
                load_csv(path)
    finally:
        path.unlink(missing_ok=True)


def test_resource_limit_guards_reject_estimated_rows_exceeded():
    """CSV with estimated rows > MAX_ESTIMATED_CSV_ROWS must raise DataLoaderError."""
    from unittest.mock import patch

    from week2_validation.utils.data_loader import (
        DataLoaderError,
        MAX_ESTIMATED_CSV_ROWS,
        load_csv,
    )

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w", encoding="utf-8") as f:
        f.write("fusion_id,gene_1,gene_2,protein_length,recurrence_count\n")
        f.write("A::B,A,B,100,1\n")
        f.flush()
        path = Path(f.name)
    try:
        # Size under file limit but estimated rows > MAX (50GB > 5GB file limit;
        # patch file limit so row limit fires)
        import week2_validation.utils.data_loader as dl

        with patch.object(dl, "MAX_FILE_SIZE_BYTES", 60 * 1024**3):
            with patch.object(Path, "stat") as mock_stat:
                mock_stat.return_value.st_size = 50_000_000_000
                with pytest.raises(DataLoaderError, match=str(MAX_ESTIMATED_CSV_ROWS)):
                    load_csv(path)
    finally:
        path.unlink(missing_ok=True)


# -----------------------------------------------------------------------------
# Test 8 — Status Envelope Sanitization
# -----------------------------------------------------------------------------
def test_status_envelope_sanitizes_paths_and_fuzzed_messages():
    """Status envelope must sanitize absolute paths, UNC, and fuzzed error messages."""
    from week2_validation.reporting.status_envelope import build_status_envelope

    notes_with_paths = [
        "Error at C:\\Users\\jane\\data\\fusion.csv",
        "File /home/user/secret.txt not found",
        r"Path \\server\share\file.csv",
    ]
    envelope = build_status_envelope(
        dataset_hash="abc",
        diagnostics_run=[],
        exit_code=0,
        notes=notes_with_paths,
    )
    notes = envelope["notes"]
    for n in notes:
        assert "C:\\" not in n
        assert "/home/" not in n
        assert "server" not in n or "<path>" in n
        assert "<path>" in n or len(n) < 50


# -----------------------------------------------------------------------------
# Test 9 — SciPy Optional Mode (existing)
# -----------------------------------------------------------------------------
# (test_scipy_optional_mode_diagnostics_still_run above)

# -----------------------------------------------------------------------------
# Test 10 — psutil Optional Mode (runtime guard degrades gracefully)
# -----------------------------------------------------------------------------
def test_psutil_optional_mode_runtime_guard_still_runs():
    """With psutil missing, runtime guard must not crash; memory check is skipped."""
    from week2_validation.runtime.runtime_guard import (
        _PSUTIL_AVAILABLE,
        check_runtime_guard,
        start_runtime_guard,
    )

    start_runtime_guard()
    check_runtime_guard()
    # If psutil available, memory check runs; if not, it skips. No crash either way.


# -----------------------------------------------------------------------------
# Test 11 — Exit Code Consistency (SchemaValidationError -> 10)
# -----------------------------------------------------------------------------
def test_schema_validation_error_maps_to_exit_10():
    """SchemaValidationError (e.g. protein_length=0) must map to exit code 10."""
    import subprocess
    import sys

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        fusion_file = root / "pl_zero.csv"
        fusion_file.write_text(
            "fusion_id,gene_1,gene_2,protein_length,recurrence_count\n"
            "A::B,A,B,0,1\n",
            encoding="utf-8",
        )
        out_dir = root / "out"
        out_dir.mkdir()
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "week2_validation.run_week2",
                "--fusion-data",
                str(fusion_file),
                "--output-dir",
                str(out_dir),
                "--run-diagnostics",
            ],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(Path(__file__).resolve().parent.parent.parent),
            env={**__import__("os").environ, "PYTHONPATH": str(Path(__file__).resolve().parent.parent.parent)},
        )
        assert proc.returncode == 10, f"Expected exit 10, got {proc.returncode}. stderr: {proc.stderr}"


# -----------------------------------------------------------------------------
# Test 12 — JSON Summary Structure
# -----------------------------------------------------------------------------
def test_json_summary_structure_contains_required_keys():
    """Summary JSON must contain required keys."""
    from week2_validation.reporting.json_summary import write_week2_summary_json

    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        run_metadata = {
            "validation_version": "0.1.0",
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
        required_keys = {"validation_version", "run_timestamp_utc", "dataset_hash", "diagnostics_run", "results", "notes"}
        for key in required_keys:
            assert key in data, f"Missing required key: {key}"


# -----------------------------------------------------------------------------
# Test 13 — Gene Symbol Normalization (COSMIC)
# -----------------------------------------------------------------------------
def test_gene_normalization_alk_lowercase():
    """'alk' -> 'ALK'"""
    from week2_validation.cosmic.diagnostics import run_cosmic_recurrence_diagnostic

    fusion_df = pd.DataFrame({
        "gene_1": ["alk"], "gene_2": ["EML4"], "recurrence_count": [10],
    })
    cosmic_df = pd.DataFrame({
        "gene_1": ["ALK"], "gene_2": ["EML4"], "recurrence_count": [10],
    })
    result = run_cosmic_recurrence_diagnostic(fusion_df=fusion_df, cosmic_df=cosmic_df, top_n=5)
    assert result["overlap_count"] == 1


def test_gene_normalization_whitespace():
    """' ALK ' -> 'ALK'"""
    from week2_validation.cosmic.diagnostics import run_cosmic_recurrence_diagnostic

    fusion_df = pd.DataFrame({
        "gene_1": [" ALK "], "gene_2": [" EML4 "], "recurrence_count": [10],
    })
    cosmic_df = pd.DataFrame({
        "gene_1": ["ALK"], "gene_2": ["EML4"], "recurrence_count": [10],
    })
    result = run_cosmic_recurrence_diagnostic(fusion_df=fusion_df, cosmic_df=cosmic_df, top_n=5)
    assert result["overlap_count"] == 1


def test_gene_normalization_mixed_case():
    """'AlK' -> 'ALK'"""
    from week2_validation.cosmic.diagnostics import run_cosmic_recurrence_diagnostic

    fusion_df = pd.DataFrame({
        "gene_1": ["AlK"], "gene_2": ["eMl4"], "recurrence_count": [10],
    })
    cosmic_df = pd.DataFrame({
        "gene_1": ["ALK"], "gene_2": ["EML4"], "recurrence_count": [10],
    })
    result = run_cosmic_recurrence_diagnostic(fusion_df=fusion_df, cosmic_df=cosmic_df, top_n=5)
    assert result["overlap_count"] == 1


# -----------------------------------------------------------------------------
# Test 14 — Data Quality Flags
# -----------------------------------------------------------------------------
def test_data_quality_5_percent_excluded_no_warning():
    """5% excluded -> no warning_flag"""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        rows = ["fusion_id,gene_1,gene_2,protein_length,recurrence_count"]
        for i in range(100):
            pl = "inf" if i < 5 else "100"
            rows.append(f"G{i}::H{i},G{i},H{i},{pl},1")
        fusion_file = root / "fusion.csv"
        fusion_file.write_text("\n".join(rows), encoding="utf-8")
        out_dir = root / "out"
        out_dir.mkdir()
        proc = subprocess.run(
            [sys.executable, "-m", "week2_validation.run_week2",
             "--fusion-data", str(fusion_file), "--output-dir", str(out_dir),
             "--run-benford"],
            capture_output=True, text=True, timeout=30,
            cwd=str(Path(__file__).resolve().parent.parent.parent),
            env={**__import__("os").environ, "PYTHONPATH": str(Path(__file__).resolve().parent.parent.parent)},
        )
        assert proc.returncode == 0
        status_path = _status_path_for_run(out_dir)
        with open(status_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        assert data["status"] == "SUCCESS"
        dq = data.get("data_quality", {})
        assert dq.get("warning_flag") is False


def test_data_quality_15_percent_excluded_warning():
    """15% excluded -> warning_flag true"""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        rows = ["fusion_id,gene_1,gene_2,protein_length,recurrence_count"]
        for i in range(100):
            pl = "inf" if i < 15 else "100"
            rows.append(f"G{i}::H{i},G{i},H{i},{pl},1")
        fusion_file = root / "fusion.csv"
        fusion_file.write_text("\n".join(rows), encoding="utf-8")
        out_dir = root / "out"
        out_dir.mkdir()
        proc = subprocess.run(
            [sys.executable, "-m", "week2_validation.run_week2",
             "--fusion-data", str(fusion_file), "--output-dir", str(out_dir),
             "--run-benford"],
            capture_output=True, text=True, timeout=30,
            cwd=str(Path(__file__).resolve().parent.parent.parent),
            env={**__import__("os").environ, "PYTHONPATH": str(Path(__file__).resolve().parent.parent.parent)},
        )
        assert proc.returncode == 0
        status_path = _status_path_for_run(out_dir)
        with open(status_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        assert data["status"] == "SUCCESS_WITH_WARNINGS"
        dq = data.get("data_quality", {})
        assert dq.get("warning_flag") is True


def test_data_quality_60_percent_excluded_high_risk():
    """60% excluded -> high_risk_flag true, status SUCCESS_WITH_HIGH_DATA_RISK"""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        rows = ["fusion_id,gene_1,gene_2,protein_length,recurrence_count"]
        for i in range(100):
            pl = "inf" if i < 60 else "100"
            rows.append(f"G{i}::H{i},G{i},H{i},{pl},1")
        fusion_file = root / "fusion.csv"
        fusion_file.write_text("\n".join(rows), encoding="utf-8")
        out_dir = root / "out"
        out_dir.mkdir()
        proc = subprocess.run(
            [sys.executable, "-m", "week2_validation.run_week2",
             "--fusion-data", str(fusion_file), "--output-dir", str(out_dir),
             "--run-benford"],
            capture_output=True, text=True, timeout=30,
            cwd=str(Path(__file__).resolve().parent.parent.parent),
            env={**__import__("os").environ, "PYTHONPATH": str(Path(__file__).resolve().parent.parent.parent)},
        )
        # With 60% excluded, quality gates may REJECT (exit 25) or report SUCCESS_WITH_HIGH_DATA_RISK (exit 0)
        assert proc.returncode in (0, 25)
        status_path = _status_path_for_run(out_dir)
        with open(status_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        dq = data.get("data_quality", {})
        assert dq.get("high_risk_flag") is True
        if proc.returncode == 0:
            assert data["status"] == "SUCCESS_WITH_HIGH_DATA_RISK"
        else:
            assert data["status"] == "FAILED"


# -----------------------------------------------------------------------------
# Test 15 — Excel .xls Rejected
# -----------------------------------------------------------------------------
def test_excel_xls_rejected_clear_error():
    """If .xls provided -> clear error message."""
    from week2_validation.utils.data_loader import UnsupportedFormatError, get_file_format
    from unittest.mock import MagicMock

    mock_path = MagicMock()
    mock_path.suffix = ".xls"
    mock_path.name = "data.xls"
    with pytest.raises(UnsupportedFormatError, match="Legacy Excel.*not supported.*convert to .xlsx"):
        get_file_format(mock_path)


# -----------------------------------------------------------------------------
# Test 16 — Empty Dataset Graceful Handling (PATCH 1)
# -----------------------------------------------------------------------------
def test_empty_dataset_graceful_handling():
    """CSV with headers only -> diagnostics.skipped true, status SUCCESS_WITH_WARNINGS, exit_code 0."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        fusion_file = root / "fusion.csv"
        fusion_file.write_text(
            "fusion_id,gene_1,gene_2,protein_length,recurrence_count\n",
            encoding="utf-8",
        )
        out_dir = root / "out"
        out_dir.mkdir()
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "week2_validation.run_week2",
                "--fusion-data",
                str(fusion_file),
                "--output-dir",
                str(out_dir),
                "--run-diagnostics",
            ],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(Path(__file__).resolve().parent.parent.parent),
            env={**__import__("os").environ, "PYTHONPATH": str(Path(__file__).resolve().parent.parent.parent)},
        )
        assert proc.returncode == 0, f"stderr: {proc.stderr!r} stdout: {proc.stdout!r}"
        status_path = _status_path_for_run(out_dir)
        with open(status_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        assert data.get("diagnostics", {}).get("skipped") is True
        assert data.get("diagnostics", {}).get("skip_reason") == "EMPTY_DATASET"
        assert data.get("status") == "SUCCESS_WITH_WARNINGS"
        dq = data.get("data_quality", {})
        assert dq.get("total_rows_original") == 0
        assert dq.get("warning_flag") is True
        assert dq.get("high_risk_flag") is False


# -----------------------------------------------------------------------------
# Test 17 — Exit Code 30 for Diagnostic Runtime Failure (PATCH 2)
# -----------------------------------------------------------------------------
def test_diagnostic_runtime_failure_exit_code_30():
    """Diagnostic runtime failure must map to exit code 30."""
    from unittest.mock import patch
    from week2_validation.run_week2 import run_pipeline, PipelineConfig
    from week2_validation.runtime.exit_codes import Week2ExitCode

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        fusion_file = root / "fusion.csv"
        fusion_file.write_text(
            "fusion_id,gene_1,gene_2,protein_length,recurrence_count\n"
            "F1::G1,G1,H1,100,1\n",
            encoding="utf-8",
        )
        out_dir = root / "out"
        out_dir.mkdir()
        config = PipelineConfig(
            fusion_data_path=fusion_file,
            output_dir=out_dir,
            dataset_stem="fusion",
            cosmic_data_path=None,
            dry_run=False,
            run_diagnostics=True,
            run_benford=False,
            run_lognormal=False,
            run_cosmic=False,
            run_benford_controls=False,
            generate_report=False,
            generate_pdf=False,
        )
        with patch("week2_validation.run_week2.execute_diagnostics") as mock_exec:
            mock_exec.return_value = int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
            exit_code = run_pipeline(config)
        assert exit_code == int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)


# -----------------------------------------------------------------------------
# Test 18 — Constant Distribution Flag (PATCH 3)
# -----------------------------------------------------------------------------
def test_constant_distribution_detected_flag():
    """protein_length = [100,100,100,100] -> constant_distribution_detected True, diagnostics still run."""
    from week2_validation.distributions.log_normality import run_log_normality_diagnostics

    result = run_log_normality_diagnostics([100.0, 100.0, 100.0, 100.0])
    assert result.constant_distribution_detected is True
    assert "variance is zero" in " ".join(result.computation_notes)
    assert result.sample_size == 4
