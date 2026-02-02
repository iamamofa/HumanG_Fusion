# Data Integrity & Statistical Validation — Documentation

*This is Week 2: Data Integrity & Statistical Validation — part of the pipeline between Week 1 (Pipeline Execution & Data Generation) and Week 3 (Zipf's Law Modeling & Interpretation).*

## 1. Overview — What Is This Layer?

**Objective:** Validate authenticity and reliability of fusion protein length data.

**In plain terms:** This layer is a data-checking step for cancer fusion gene datasets. It takes your data file, makes a safe copy, runs checks you request, and tells you what it found. It does not approve or reject your data—it only describes it.

### What This Layer Does

- **Checks your data file** — Confirms the file exists, can be read, and has the required columns (fusion_id, gene_1, gene_2, protein_length, recurrence_count).
- **Makes a locked copy** — Creates an unchanged snapshot of your data before any analysis. This way, if someone edits the original file, your results stay the same.
- **Runs the checks you ask for** — Distribution stats, Benford’s Law, log-normality, or COSMIC comparison. Results are shown on screen and written to output files.
- **Shows results on screen** — Prints statistics; does not save plots or charts to disk.

### What This Layer Does NOT Do

- **Does not pass or fail your data** — It does not approve or reject based on statistics.
- **Does not draw conclusions** — No hypothesis tests, p-value interpretation, or biological claims.
- **Does not change your data** — It only reads and reports; it does not modify values.

### Why It Exists

This Data Integrity & Statistical Validation layer gives you a clear, reproducible checkpoint. It locks your data so downstream analysis uses the exact same copy, and it records what was checked and what was found.

---

## 2. Project Layout — What’s in the week2_validation Folder

*For developers: the main file you run is `run_week2.py`. The folders below contain the code that loads data, runs checks, and writes outputs.*

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
│   └── week1_adapter.py  # Converts Week 1 (Pipeline Execution & Data Generation) output to this layer's schema
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
| `run_week2.py` | The main program — this is what you run |
| `adapters/` | Converts Week 1 (Pipeline Execution & Data Generation) column names to the format this layer expects |
| `config/` | Settings and rules (e.g. thresholds.yaml) |
| `schemas/` | Defines required columns and constraints for your data |
| `utils/` | Loads files, creates the locked copy, checks state |
| `distributions/` | Distribution and log-normality checks |
| `benford/` | Benford’s Law first-digit checks |
| `cosmic/` | COSMIC database comparison |
| `reporting/` | Writes status JSON, cleaned dataset, and certification files |
| `runtime/` | Exit codes and safety checks |
| `tests/` | Integration tests (for developers) |
| `outputs/` | Where results are saved (week2_status.json, week2_cleaned_dataset.csv, etc.) — you choose the folder with `--output-dir` |
| `frozen_inputs/` | Locked copies of your data (created automatically; you can delete to clear old snapshots) |

### Dependencies (`requirements.txt` / `requirements.lock.txt`)

- **Required:** pandas, numpy, matplotlib, PyYAML (install these to run the pipeline)
- **Optional:** scipy (for some statistics), openpyxl (for Excel files), pyarrow (for Parquet), psutil (for memory checks)
- **Development:** pytest (for running tests)

*If optional packages are missing, the pipeline still runs but may omit some statistics or features.*

---

## 2.1 Limits — File Size, Rows, and Runtime

*The pipeline has limits to avoid crashing or running too long. Very large files or very long runs may be rejected or warned about.*

| Parameter | Safe | Warning | Hard Limit (rejected) |
|-----------|-----------|--------------|------------|
| **File size** | < 1 GB | 1–5 GB (warning) | 5 GB (rejected) |
| **CSV/TSV rows (est.)** | < 10M | 10M–100M (warning) | 200M (rejected) |
| **Runtime** | < 100 s typical | 100–3600 s | 3600 s (RuntimeError) |
| **Memory (RSS)** | < 2 GB typical | 2–8 GB | 8 GB (RuntimeError, requires psutil) |

When you approach these limits, the pipeline may log warnings. If you exceed the hard limits, it will stop and report an error. The 5 GB file size limit aligns with the 8 GB RAM runtime guard to prevent memory exhaustion.

---

## 2.2 Running with a Real-World Dataset

This section describes how to run the Data Integrity & Statistical Validation pipeline on your own fusion dataset, what happens step by step, and what output to expect.

### Input Requirements

- **Format:** CSV, TSV, JSON, Parquet, or Excel (.xlsx)
- **Required columns:** Your spreadsheet or file must have these column names: `fusion_id`, `gene_1`, `gene_2`, `protein_length`, `recurrence_count`
- **Valid values:** `protein_length` must be greater than 0; `recurrence_count` must be 0 or more
- **No empty cells** in required columns

- Be in a supported format: CSV, TSV, JSON, Parquet, or Excel (.xlsx)
- Include these required columns: `fusion_id`, `gene_1`, `gene_2`, `protein_length`, `recurrence_count`
- Have `protein_length` > 0 and `recurrence_count` >= 0 for all rows
- Have no nulls in required columns

If your data uses Week 1 (Pipeline Execution & Data Generation) column names (`geneA`, `geneB`, `samples_detected` or `recurrence_frequency`), the pipeline adapts them automatically.

### How to Run

1. **Install required packages** (run once, from the main project folder):

   ```bash
   pip install -r week2_validation/requirements.txt
   ```

2. **Run the pipeline** (from the main project folder, e.g. `HumanG_Fusion/`):

   ```bash
   python -m week2_validation.run_week2 --fusion-data /path/to/your_fusion.csv --output-dir /path/to/outputs --run-diagnostics
   ```

   Replace `/path/to/your_fusion.csv` with your actual data file path, and `/path/to/outputs` with where you want the results saved.

3. **Optional extra checks** — You can add any of these flags:
   - `--run-benford` — First-digit pattern check (Benford’s Law)
   - `--run-lognormal` — Log-scale distribution check
   - `--run-cosmic` — Compare with COSMIC database (add `--cosmic-data /path/to/cosmic.tsv` if you have it)
   - `--dry-run` — Only check that files exist; do not run any analysis

### What Happens Step by Step

| Step | What happens |
|------|--------------|
| 1 | **Read your command** — Check that the data file and output folder paths are valid |
| 2 | **Load settings** — Read configuration from `thresholds.yaml` |
| 3 | **Make a locked copy** — Create a snapshot of your data under `frozen_inputs/`. All checks use this copy, not your original file. |
| 4 | **Convert column names if needed** — If your file uses Week 1 (Pipeline Execution & Data Generation) names (geneA, geneB, etc.), convert to the expected format |
| 5 | **Validate your data** — Load the file, check that all required columns exist, and confirm the output folder can be written to |
| 6 | **Check data quality** — Count how many rows have valid `protein_length`. If the file is empty, skip checks and report a warning. |
| 7 | **Run the checks you requested** — Run distribution, Benford, log-normality, or COSMIC checks. Results are printed on screen. |
| 8 | **Write output files** — Save `week2_status.json` (always), and on success: `week2_cleaned_dataset.csv` and `week2_dataset_certification.json` |

### What You See on Screen

```
============================================================
Data Integrity & Statistical Validation
============================================================

Loading configuration...
  Configuration loaded from: thresholds.yaml
  ...
Ensuring dataset is frozen for Data Integrity & Statistical Validation diagnostics...
Dataset not frozen — freezing input for diagnostics
Freeze complete — diagnostics may now proceed
  Frozen data location: .../frozen_inputs/<hash>/your_fusion.csv

Validating inputs...
Validating fusion data: your_fusion.csv
  Loaded N records with 5 fields
  Schema validation: PASSED
Output directory validated: /path/to/outputs

============================================================
RUNNING DIAGNOSTICS
============================================================
...
Starting diagnostic analysis...
Running diagnostics on N protein length values...

Diagnostic Results:
  Sample size: N
  Min: ...
  Max: ...
  Mean: ...
  ...
Status written to /path/to/outputs/week2_status.json
```

### Output Files (in your chosen output folder)

| File | When created | What it contains |
|------|--------------|------------------|
| `week2_status.json` | Every run | Run summary: version, dataset ID, exit code, status, which checks ran, data quality info, run time |
| `week2_cleaned_dataset.csv` | On success | A copy of your validated data (unchanged) |
| `week2_dataset_certification.json` | On success | A record that the dataset passed validation and is certified for downstream use |
| `week2_adapted_fusion.csv` | When Week 1 (Pipeline Execution & Data Generation) format is detected | Your data after column names were converted |

**Frozen snapshot** — Created under `week2_validation/frozen_inputs/<hash>/`:

The pipeline automatically creates a locked copy of your input under `week2_validation/frozen_inputs/`. This copy is never changed. It includes:

- A copy of your input file
- A manifest file (source path, hash, timestamp)
- A marker file indicating the data is locked

---

## 3. High-Level Workflow Diagram

**In plain terms:** You give the program your data file and an output folder. It checks the files, makes a locked copy, runs the checks you asked for, and shows results on screen and writes output files.

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
    [4] Is the data ready for analysis?
       |
       +---> No  ---> Stop. Message: "Please freeze the dataset first."
       |
       v
       Yes
       |
       v
    [5] Did you request any checks? (e.g. --run-diagnostics)
       |
       +---> No  ---> Stop. Message: "Add --run-diagnostics (or similar) to run checks."
       |
       v
       Yes
       |
       v
    [6] Run the requested checks on the locked copy
       |
       v
    [7] Show results on screen and write output files
       |
       v
    [END]
```

### Why a Copy Is Made

```
    Your original data file
           |
           |  (Make a locked copy)
           v
    [LOCKED COPY]  <-- Analysis always uses this. The original is never touched again.
           |
           v
    All checks read from the copy
```

If someone edits the original file while the pipeline runs, your results stay correct because they are based on the locked copy.

---

## 4. Step-by-Step Execution Trace (Technical Detail)

*For a simpler overview, see Section 2.2 “Running with a Real-World Dataset” and Section 3 “High-Level Workflow.”*

When a user runs (from the project root, e.g. `HumanG_Fusion/`):

```bash
python -m week2_validation.run_week2 --fusion-data <path/to/fusion.csv> --output-dir <path/to/output> [--run-diagnostics] [--run-benford] [--run-lognormal] [--run-cosmic] [--dry-run]
```

Example:

```bash
python -m week2_validation.run_week2 --fusion-data ./my_fusion.csv --output-dir ./outputs --run-diagnostics
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

## 5. How the Locked Copy Works (Technical Detail)

### Why This Layer Makes Its Own Copy

The Data Integrity & Statistical Validation layer does not assume that any earlier step locked the data. To keep results reproducible and auditable, it makes its own locked copy. All checks use that copy; the original file is not used after the copy is made.

### What "Frozen" Means

A dataset is frozen when:

1. A copy of your input file exists under `frozen_inputs/<hash>/`
2. A `manifest.json` records the source path, file hash, timestamp, and version
3. A `.frozen` marker file is present

Same input file always produces the same frozen folder. Running twice does not create duplicates. If the copy is corrupted, it is removed and recreated.

### Why the Original File Is Not Used After the Copy

The original file could be changed (e.g. overwritten) while the pipeline runs. The locked copy is a snapshot in time. All results are tied to that copy, which can be verified using the manifest and hash.

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

## 7. Diagnostic Flags — What Each One Does

| Flag | What it does |
|------|--------------|
| `--run-diagnostics` | Run distribution checks and show summary stats (min, max, mean, etc.) |
| `--run-benford` | Check first-digit pattern in numbers (Benford’s Law) |
| `--run-lognormal` | Check how values are spread (log-scale distribution) |
| `--run-cosmic` | Compare your fusion list to the COSMIC database |
| `--run-benford-controls` | Internal test with fake data (for developers; does not use your data) |
| `--dry-run` | Only check that files exist; do not run any analysis |

All checks print results to the screen. Output files (e.g. `week2_status.json`) are written to your chosen output folder.

---

## 8. What This Layer Will Never Do (By Design)

- **Does NOT** approve or reject your data based on statistics
- **Does NOT** pass or fail based on Benford, log-normality, or COSMIC results
- **Does NOT** interpret biology or draw conclusions about data quality
- **Does NOT** fit models or make predictions
- **Does NOT** save plots or charts to disk
- **Does NOT** prepare or filter files for the next stage automatically

This layer only describes your data. Week 3 (Zipf's Law Modeling & Interpretation) or the user decides what to do with that information.

---

## 9. Where This Layer Fits in the Pipeline

The Data Integrity & Statistical Validation layer receives a fusion data file and optionally a COSMIC reference file. It does not depend on how Week 1 (Pipeline Execution & Data Generation) works; it makes its own locked copy.

How this layer fits into the overall pipeline:

```
    [Raw RNA-Seq data]
            |
            v
    [Week 1: Pipeline Execution & Data Generation]  Process data → Find fusion genes
            |
            v
    [Week 2: Data Integrity & Statistical Validation]  Check data → Lock a copy → Run diagnostics  <-- You are here
            |
            v
    [Week 3: Zipf's Law Modeling & Interpretation]  Final analysis and results
```

This layer does not produce formal "approved" outputs. It provides diagnostic information. Week 3 (Zipf's Law Modeling & Interpretation) or the user decides how to use that information.

---

## 10. Current Status (Fact-Only)

### Implemented

- CLI with `--fusion-data`, `--output-dir`, `--cosmic-data`, `--run-diagnostics`, `--run-benford`, `--run-lognormal`, `--run-cosmic`, `--run-benford-controls`, `--dry-run`
- Week 1 (Pipeline Execution & Data Generation) adapter — maps geneA/geneB, recurrence_frequency/samples_detected to this layer's schema
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
- Output directory (user-specified via `--output-dir`): `week2_status.json` (every run), `week2_cleaned_dataset.csv` and `week2_dataset_certification.json` (on success), `week2_adapted_fusion.csv` (when Week 1 (Pipeline Execution & Data Generation) format is detected).
- No files written by diagnostic modules (no plots, no reports).
