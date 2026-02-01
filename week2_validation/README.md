# Week 2 Validation Pipeline — Technical Documentation

## 1. Overview — What Week 2 Is (As Implemented)

Week 2 is a **diagnostic-only** data integrity and statistical validation stage for cancer genomics fusion datasets. It sits between Week 1 (upstream RNA-Seq/fusion processing) and Week 3 (downstream analysis).

### What Week 2 Actually Does

- **Validates input files** — Checks that fusion and optional COSMIC reference files exist, are readable, and conform to required schema.
- **Enforces defensive freeze** — Creates or validates an immutable snapshot of the input dataset before any diagnostics run.
- **Runs requested diagnostics** — When the dataset is frozen and the user explicitly requests diagnostics, computes descriptive statistics and structural characterizations.
- **Reports results to stdout** — Prints computed statistics; does not save plots or reports to disk.

### Guarantees

- **No inference** — No hypothesis testing, p-value interpretation, or statistical conclusions.
- **No gating** — Does not approve or reject datasets; does not halt based on statistical thresholds.
- **No biological claims** — All outputs are descriptive characterizations for exploratory inspection.
- **No file writes from diagnostics** — Histograms and statistics are computed in memory; figures are not saved.

### Why It Exists

Week 2 provides a reproducible checkpoint for data integrity. It ensures that downstream stages (Week 3) operate on a defensively frozen, auditable snapshot. It does not depend on Week 1 internals; it enforces its own freeze logic.

---

## 2. Project Layout — File Structure

```
week2_validation/
│
├── __init__.py           # Package root; exports __version__
├── run_week2.py          # Main program — run this to start the pipeline
├── README.md             # This documentation
├── requirements.txt      # Python packages needed to run
├── requirements.lock.txt # Locked dependency versions (reproducible installs)
│
├── adapters/             # Upstream format compatibility
│   └── week1_adapter.py  # Converts Week 1 output to Week 2 schema (geneA/geneB → gene_1/gene_2)
│
├── config/               # Settings and rules
│   ├── thresholds.yaml   # Configuration (sample sizes, test settings, etc.)
│   └── loader.py         # Reads and validates the config file
│
├── schemas/              # Input schema contracts
│   └── fusion_schema.yaml   # Required columns and constraints for fusion data
│
├── utils/                # Shared helpers
│   ├── __init__.py       # Package init
│   ├── data_loader.py    # Loads fusion and COSMIC data from files
│   ├── freeze.py         # Creates the locked copy of your data
│   └── state.py          # Checks if data is marked as "ready for analysis"
│
├── distributions/        # Distribution checks
│   ├── visualize.py      # Charts and summary stats (--run-diagnostics)
│   └── log_normality.py  # Log-scale spread check (--run-lognormal)
│
├── benford/              # Benford's Law checks
│   ├── diagnostics.py    # First-digit pattern analysis (--run-benford, --run-benford-controls)
│   ├── benford_controls.py   # Test data generators (standalone)
│   └── benford_test.py   # Chi-squared test utilities (standalone)
│
├── cosmic/               # COSMIC comparison
│   └── diagnostics.py    # Compare your fusion list to COSMIC (--run-cosmic)
│
├── reporting/            # Structured output
│   ├── json_summary.py   # Optional: writes week2_integrity_summary.json when invoked (not called by default pipeline)
│   ├── status_envelope.py   # Builds machine-readable status envelope (version, exit_code, status)
│   └── approved_dataset_writer.py   # Writes week2_cleaned_dataset.csv and week2_dataset_certification.json
│
├── runtime/              # Survivability layer
│   ├── exit_codes.py     # Deterministic exit code contract (Week2ExitCode)
│   ├── runtime_guard.py  # Lightweight runtime limits (duration, memory via psutil)
│   └── safe_runner.py    # Atomic failure containment wrapper
│
├── tests/                # Integration tests
│   └── test_week2_integration.py   # Freeze, schema, SciPy optional, JSON summary
│
├── outputs/              # Runtime: pipeline output (user-specified via --output-dir)
│
└── frozen_inputs/        # Runtime: immutable dataset snapshots (created on run)
    └── <hash_prefix>/    # One dir per frozen input (hash-based)
```

| Folder / File | Purpose |
|---------------|---------|
| `run_week2.py` | Entry point — the program you run |
| `adapters/` | Week 1 format compatibility (geneA/geneB → gene_1/gene_2) |
| `config/` | Configuration and thresholds |
| `schemas/` | Input schema contracts (fusion_schema.yaml) |
| `utils/` | File loading, freeze logic, state checks |
| `distributions/` | Distribution and log-normality diagnostics |
| `benford/` | Benford's Law diagnostics and test utilities |
| `cosmic/` | COSMIC database comparison |
| `reporting/` | Structured JSON output, status envelope, approved dataset writer |
| `runtime/` | Exit codes, runtime guard, safe execution wrapper |
| `tests/` | Integration tests (pytest) |
| `outputs/` (or user-specified) | Runtime pipeline output directory (week2_status.json, week2_cleaned_dataset.csv, week2_dataset_certification.json) |
| `frozen_inputs/` | Runtime directory for frozen dataset snapshots (created on first run; safe to delete to clear old snapshots) |

**Using real data:** Point `--fusion-data` to your fusion CSV and `--output-dir` to where you want outputs. The pipeline creates `frozen_inputs/` automatically when it runs (one snapshot per input file). No simulation or test artifacts are required.

### Dependencies (`requirements.txt` / `requirements.lock.txt`)

- **Required:** pandas, numpy, matplotlib, PyYAML
- **Optional (statistical):** scipy (graceful degradation if missing; Benford p-value, log-normality Anderson-Darling, bias-corrected skewness)
- **Optional (file formats):** openpyxl (Excel), pyarrow (Parquet)
- **Optional (runtime):** psutil (memory limit checks in `runtime_guard`; if missing, runtime guard skips memory check)
- **Development:** pytest

**Optional dependency behavior:**
- **scipy missing:** Benford still computes chi-squared; p-value omitted. Log-normality uses manual KS; Anderson-Darling omitted. Visualization uses NumPy skewness fallback. Status envelope reports `scipy_available: false`.
- **psutil missing:** Runtime guard checks only elapsed time (no memory limit). Status envelope reports `psutil_available: false`. No crash; pipeline runs normally.

---

## 2.1 Operational Envelope (Safe Operating Zone)

| Parameter | Safe Zone | Caution Zone | Hard Limit |
|-----------|-----------|--------------|------------|
| **File size** | < 1 GB | 1–5 GB (warning) | 5 GB (rejected) |
| **CSV/TSV rows (est.)** | < 10M | 10M–100M (warning) | 200M (rejected) |
| **Runtime** | < 100 s typical | 100–3600 s | 3600 s (RuntimeError) |
| **Memory (RSS)** | < 2 GB typical | 2–8 GB | 8 GB (RuntimeError, requires psutil) |

When inputs approach limits, the pipeline logs non-blocking warnings. Exceeding hard limits raises `DataLoaderError` or `RuntimeError` with deterministic exit codes. See `utils/data_loader.py` (MAX_FILE_SIZE_BYTES, MAX_ESTIMATED_CSV_ROWS) and `runtime/runtime_guard.py` (MAX_RUNTIME_SECONDS, MAX_MEMORY_MB).

---

## 3. High-Level Workflow Diagram

### Simple Overview

```
    [START]
       |
       v
    [1] Run the program (you provide: data file, output folder)
       |
       v
    [2] Check that files exist and are readable
       |
       v
    [3] Make a safe, locked copy of your data (cannot be changed during analysis)
       |
       v
    [4] Ask: "Is the data marked as ready for analysis?"
       |
       +---> No  ---> Stop. Message: "Please freeze the dataset first."
       |
       v
       Yes
       |
       v
    [5] Ask: "Did you request any diagnostic checks?"
       |
       +---> No  ---> Stop. Message: "Add --run-diagnostics (or similar) to run checks."
       |
       v
       Yes
       |
       v
    [6] Run the requested checks on the safe copy
       |
       v
    [7] Show results on screen (no files saved)
       |
       v
    [END]
```

### Defensive Freeze: Why a Copy Is Used

```
    Your original data file
           |
           |  (Make a safe copy)
           v
    [LOCKED COPY]  <-- Analysis always uses this. The original is never touched again.
           |
           v
    All diagnostic checks read from the copy
```

This ensures that if someone changes the original file during analysis, the results stay consistent because they are based on the locked copy.

---

## 4. Step-by-Step Execution Trace

When a user runs:

```bash
python -m week2_validation.run_week2 [flags]
```

the following occurs in order:

1. **CLI parsing** — `create_argument_parser()` defines required (`--fusion-data`, `--output-dir`) and optional arguments. If no arguments are provided, help is printed and the process exits with code 1.

2. **Benford controls shortcut** — If `--run-benford-controls` is present and no other diagnostic flags (`--run-diagnostics`, `--run-benford`, `--run-lognormal`, `--run-cosmic`) are present, `execute_benford_controls()` runs immediately. This uses synthetic data only and does not require `--fusion-data` or `--output-dir`. The process then exits.

3. **Required argument enforcement** — `validate_arguments()` ensures `--fusion-data` and `--output-dir` exist, are valid paths, and (for files) that they exist on disk. Raises `ConfigurationError` on failure.

4. **Config loading** — `load_config()` reads `thresholds.yaml` from `week2_validation/config/thresholds.yaml`. All sections (dataset_frozen_required, benford, distribution_tests, visualization, etc.) are validated. Raises `ConfigError` on failure.

5. **Defensive freeze** — `ensure_frozen_input()` is called with the fusion data path and `frozen_root = week2_validation/frozen_inputs/`. If the input hash matches an existing frozen snapshot, that path is reused. Otherwise, a new snapshot is created: copy file, write manifest, create `.frozen` marker. All subsequent steps use the frozen path.

6. **Freeze-state logic** — `resolve_freeze_state(flag_file_path)` checks whether `week2_validation/.frozen` exists. If it exists, `is_frozen=True`. If not, `is_frozen=False` (given `cli_frozen=False` as default). The flag file is authoritative; there is no CLI override for "frozen."

7. **Input validation** — `validate_inputs()` loads fusion data (and COSMIC if provided) from the frozen path, validates schema (required columns: `fusion_id`, `gene_1`, `gene_2`, `protein_length`, `recurrence_count`; see `schemas/fusion_schema.yaml`), and confirms the output directory is writable.

8. **Exit conditions**:
   - `--dry-run`: Exit 0 after validation.
   - Not frozen: Print message, exit 0.
   - Frozen but no diagnostic flags: Print message and suggested flags, exit 0.
   - Frozen and diagnostic flags: Run each requested diagnostic, collect exit codes, return the last non-zero if any.

9. **Diagnostic dispatch** — Each diagnostic is lazy-imported when its flag is set. Data is loaded from the frozen path. The `protein_length` column is extracted for distribution, Benford, and log-normality diagnostics. COSMIC uses the full fusion DataFrame (requires `gene_1`, `gene_2`, `recurrence_count` at runtime).

---

## 5. Defensive Freeze Logic (Critical Section)

### Why Week 2 Owns Freeze Logic

Week 2 cannot assume that Week 1 froze the data. To ensure reproducibility and auditability, Week 2 implements its own freeze mechanism. All diagnostics operate on an immutable snapshot; the original input is never trusted after the freeze step.

### What "Frozen" Means

A dataset is frozen when:

1. A copy of the input file exists under `frozen_inputs/<hash_prefix>/`
2. A `manifest.json` records source path, hash, timestamp, and week2 version
3. A `.frozen` marker file is present

The hash is SHA-256 of the file contents. The directory name uses the first 12 characters of the hash. Same input always yields the same frozen directory (deterministic and idempotent).

### How `ensure_frozen_input()` Works

1. Validate that the input path exists and is a file.
2. Search `frozen_root` for an existing directory whose manifest hash matches the current input hash.
3. If found and valid: return that directory.
4. If not found: compute hash, create `frozen_root/<hash_prefix>/`, copy the file, write manifest, create `.frozen`, validate the snapshot.
5. On any failure: raise `RuntimeError`; diagnostics do not run.

### Idempotency

- Same input file → same hash → same directory. Running twice does not create duplicates.
- If the frozen snapshot is corrupted, it is removed and recreated.

### What Diagnostics Are Allowed to See

Diagnostics receive the path to the frozen data file via `get_frozen_data_path(frozen_dir)`. They load data from that path only. The original `--fusion-data` path is not used for analysis after the freeze.

### Why Original Input Is Never Trusted

The original path may change (e.g., file overwritten) between validation and analysis. The freeze creates a point-in-time snapshot. All diagnostic results are tied to the frozen copy, which is verifiable via the manifest and hash.

---

## 6. Diagnostic Modules

### 6.1 Visualization (`distributions/visualize.py`)

**What it computes:** Descriptive statistics (count, min, max, mean, median, std, skewness) and histogram figures (linear and optional log-scale) for protein length data.

**Required inputs:** Array-like numeric data (list, numpy array, or pandas Series) of strictly positive, finite values.

**Optional dependencies:** SciPy (for bias-corrected skewness; otherwise NumPy fallback), pandas (for Series input).

**Real vs synthetic:** Operates on whatever data is passed; used with real frozen data when `--run-diagnostics` is set.

**Writes files:** No. Returns `DiagnosticResult` with in-memory `Figure` objects.

**Performs inference:** No.

### 6.2 Log-Normality (`distributions/log_normality.py`)

**What it computes:** Kolmogorov-Smirnov statistic (and optionally Anderson-Darling via SciPy) describing deviation from a fitted log-normal distribution. KS p-value is intentionally not reported (parameters estimated from same data).

**Required inputs:** Numeric array of strictly positive, finite values.

**Optional dependencies:** SciPy (for Anderson-Darling and KS via `kstest`; without SciPy, KS is computed manually).

**Real vs synthetic:** Used with real frozen data when `--run-lognormal` is set.

**Writes files:** No.

**Performs inference:** No.

### 6.3 Benford Diagnostics (`benford/diagnostics.py`)

**What it computes:** First significant digit distribution, observed vs expected (Benford) frequencies, chi-squared statistic, optional p-value. Heuristic applicability assessment (sample size, scale span).

**Required inputs:** Numeric array; only strictly positive finite values are used.

**Optional dependencies:** SciPy (for chi-squared p-value), pandas (for Series input).

**Real vs synthetic:** Used with real frozen data when `--run-benford` is set.

**Writes files:** No.

**Performs inference:** No.

### 6.4 Benford Controls (`benford/diagnostics.py`)

**What it computes:** Implementation self-test using synthetic Benford-compliant (log-uniform) and non-Benford (uniform) data. Verifies that the Benford pipeline behaves as expected on known inputs.

**Required inputs:** None (synthetic data generated internally).

**Optional dependencies:** NumPy (for `default_rng`), SciPy (for p-value in controls).

**Real vs synthetic:** Synthetic only. Triggered by `--run-benford-controls`.

**Writes files:** No.

**Performs inference:** No.

### 6.5 COSMIC Diagnostics (`cosmic/diagnostics.py`)

**What it computes:** Descriptive rank-order comparison between fusion and COSMIC recurrence data. Counts overlapping pairs, rank discrepancies; reports top N discrepancies by absolute rank difference.

**Required inputs:** `fusion_df` (must conform to fusion schema: `fusion_id`, `gene_1`, `gene_2`, `protein_length`, `recurrence_count`) and optionally `cosmic_df` with columns `gene_1`, `gene_2`, `recurrence_count`.

**Optional dependencies:** None beyond pandas (provided by data loader).

**Real vs synthetic:** Used with real fusion data and optional COSMIC reference when `--run-cosmic` is set.

**Writes files:** No.

**Performs inference:** No. No statistical tests, correlations, or validation conclusions.

### 6.6 Reporting (`reporting/`)

**`json_summary.py`** — Provides `write_week2_summary_json()` to write a structured JSON file (`week2_integrity_summary.json`) when explicitly invoked. The main pipeline does not call this; it is available for optional use.

**`status_envelope.py`** — Builds the machine-readable status envelope with `week2_version`, `run_timestamp_utc`, `dataset_hash`, `diagnostics_run`, `exit_code`, and `status` (SUCCESS/FAILED/PARTIAL). Written as `week2_status.json` to the output directory on every run.

**`approved_dataset_writer.py`** — When validation succeeds, writes `week2_cleaned_dataset.csv` (copy of frozen validated dataset) and `week2_dataset_certification.json` (certification record for power-law modeling) to the output directory. Additive only; does not affect pipeline success/failure.

**Required inputs (json_summary, when invoked):** `output_dir` (Path), `run_metadata` (dict), `diagnostic_results` (dict).

**Optional dependencies:** None (uses stdlib `json` and `datetime`).

**Writes files:** `status_envelope` and `approved_dataset_writer` are invoked by the pipeline. `json_summary` writes only when explicitly called.

**Performs inference:** No.

---

## 7. Diagnostic Flag Matrix

| Flag | What it does | Needs locked data? | Uses real data? |
|------|--------------|--------------------|-----------------|
| `--run-diagnostics` | Show distribution charts and summary stats | Yes | Yes |
| `--run-benford` | Check first-digit pattern in numbers | Yes | Yes |
| `--run-lognormal` | Check how values are spread (log-scale) | Yes | Yes |
| `--run-cosmic` | Compare your fusion list to COSMIC database | Yes | Yes |
| `--run-benford-controls` | Self-test using fake data (no real data used) | No | No (test data only) |
| `--dry-run` | Only check that files exist; do not run any analysis | No | N/A |

All checks print results to the screen. No files are saved.

---

## 8. What Week 2 Will Never Do (By Design)

- **Does NOT** approve or reject datasets based on statistics.
- **Does NOT** halt or fail based on Benford, log-normality, or COSMIC results.
- **Does NOT** interpret biology or draw conclusions about data quality.
- **Does NOT** fit models for inference or prediction.
- **Does NOT** generate final reports or save plots to disk.
- **Does NOT** prepare Week 3 inputs automatically (e.g., filtered files, approval flags).

Week 2 is a diagnostic layer. All interpretation and decision-making is explicitly deferred to downstream stages and human reviewers.

---

## 9. Relationship to Week 1 and Week 3

### Week 1

Week 2 does not depend on Week 1 internals. It receives a fusion data file path and optionally a COSMIC reference path. The defensive freeze exists because Week 1’s behavior (including whether it freezes data) is unknown. Week 2 creates or validates its own frozen snapshot.

### Week 2 Position

How Week 2 fits into the overall pipeline:

```
    [Raw RNA-Seq data]
            |
            v
    [Week 1]  Process data → Find fusion genes
            |
            v
    [Week 2]  Check data → Lock a copy → Run diagnostics  <-- You are here
            |
            v
    [Week 3]  Final analysis and results
```

### Week 3

Week 3 must not bypass Week 2. The defensive freeze and freeze-state check ensure that analysis uses a reproducible snapshot. Week 2 does not produce formal "approved" outputs; it provides diagnostic information. Week 3 (or the operator) decides how to use that information.

---

## 10. Current Status (Fact-Only)

### Implemented

- CLI with `--fusion-data`, `--output-dir`, `--cosmic-data`, `--run-diagnostics`, `--run-benford`, `--run-lognormal`, `--run-cosmic`, `--run-benford-controls`, `--dry-run`
- Week 1 adapter (`adapters/week1_adapter.py`) — maps geneA/geneB, recurrence_frequency/samples_detected to Week 2 schema
- Configuration loading from `thresholds.yaml`
- Schema contract in `schemas/fusion_schema.yaml` (required columns: `fusion_id`, `gene_1`, `gene_2`, `protein_length`, `recurrence_count`)
- Defensive freeze via `ensure_frozen_input()` and `get_frozen_data_path()`
- Freeze state resolution from flag file `week2_validation/.frozen`
- Input validation (path checks, schema validation for fusion data)
- Distribution diagnostics (visualization, descriptive statistics)
- Log-normality diagnostics (KS, optional Anderson-Darling)
- Benford diagnostics (FSD distribution, chi-squared, applicability heuristic)
- Benford controls (synthetic positive/negative)
- Benford test utilities (`benford_test.py`) — chi-squared statistic and p-value
- COSMIC rank-order diagnostic
- Reporting: `reporting/json_summary.py` (structured JSON output), `reporting/status_envelope.py` (status envelope), `reporting/approved_dataset_writer.py` (cleaned dataset and certification JSON)
- Runtime layer: `runtime/exit_codes.py`, `runtime/runtime_guard.py`, `runtime/safe_runner.py`
- Integration tests (`tests/test_week2_integration.py`) for freeze, schema, SciPy optional mode, JSON summary

### Runnable

- The pipeline runs via `python -m week2_validation.run_week2`
- All diagnostic modules support lazy import; missing optional deps cause graceful degradation with notes.
- Run integration tests via `python -m pytest week2_validation/tests/` (requires pytest).

### Executed

- Execution is user-driven. No automated execution or scheduling is implemented in the codebase.

### Deliverables Produced

- Stdout output: validation messages, diagnostic statistics, disclaimers.
- Frozen snapshot directory under `week2_validation/frozen_inputs/<hash_prefix>/` when freeze is performed.
- Output directory (user-specified via `--output-dir`): `week2_status.json` (every run), `week2_cleaned_dataset.csv` and `week2_dataset_certification.json` (on success), `week2_adapted_fusion.csv` (when Week 1 format is detected).
- No files written by diagnostic modules (no plots, no reports).
