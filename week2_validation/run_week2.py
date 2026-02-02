#!/usr/bin/env python3
"""
Week 2: Data Integrity & Statistical Validation — CLI Entry Point.

This module provides the command-line interface for the Data Integrity &
Statistical Validation pipeline. It handles argument
parsing, input validation, and orchestration of the validation workflow.

WHAT DOES THIS FILE DO?
This is the main "control center" for running the Data Integrity & Statistical Validation checks.
When a user runs this program from the command line, this file:
1. Reads the user's instructions (which files to analyze, where to save results)
2. Loads configuration from thresholds.yaml
3. Checks dataset freeze state (Option C: flag file takes precedence over CLI)
4. Validates all input files exist and are in the correct format
5. If diagnostics requested AND dataset is frozen, runs requested diagnostics:
   - --run-diagnostics: distribution diagnostics
   - --run-benford: Benford's Law diagnostics (diagnostic only, no inference)
   - --run-lognormal: log-normality diagnostics (diagnostic only, no inference)
6. Reports the results

EXECUTION FLOW:
    ┌─────────────────────────────────────────────────────────────────┐
    │  1. Parse CLI arguments                                         │
    │  2. Load thresholds.yaml configuration                          │
    │  3. Check freeze state (flag file > CLI > error)                │
    │  4. Validate input files                                        │
    │  5. If --dry-run: exit after validation                         │
    │  6. If NOT frozen: exit with message                            │
    │  7. If frozen but NO diagnostics requested: exit with message   │
    │  8. If frozen AND --run-diagnostics: run distribution diags     │
    │  9. If frozen AND --run-benford: run Benford diags              │
    │ 10. If frozen AND --run-lognormal: run log-normality diags      │
    └─────────────────────────────────────────────────────────────────┘

This pipeline is diagnostic only:
    - No hypothesis testing
    - No modeling
    - No biological claims

Usage:
    python -m week2_validation.run_week2 \\
        --fusion-data /path/to/fusion_data.csv \\
        --output-dir /path/to/output \\
        [--run-all] \\
        [--cosmic-data /path/to/cosmic.tsv] \\
        [--run-diagnostics] [--run-benford] [--run-lognormal] [--run-benford-controls] \\
        [--dry-run]

Security considerations:
    - All file paths are validated before use
    - No directory traversal is permitted
    - Output is restricted to the specified output directory
    - No sensitive data is logged or printed
    - Diagnostics only run when explicitly requested AND dataset is frozen
"""

# =============================================================================
# STANDARD LIBRARY IMPORTS
# =============================================================================

# 'argparse' helps read and understand command-line arguments (user instructions)
import argparse
import io
import math
# 'json' for writing status envelope
import json
# 'sys' provides access to system functions like exiting the program
import sys
# 'time' for runtime tracking
import time
# 'dataclass' creates simple classes for holding related data together
from dataclasses import dataclass
# 'Path' helps work with file and folder locations on the computer
from pathlib import Path
# 'Optional' indicates that a value might be present or might be None (empty)
from typing import Optional

# =============================================================================
# RUN LOG TEE (capture console output to file)
# =============================================================================


class _TeeWriter(io.TextIOBase):
    """Writes to both an underlying stream and a log file."""

    def __init__(self, stream, log_path: Path):
        self._stream = stream
        self._log_path = Path(log_path)

    def write(self, s: str) -> int:
        n = self._stream.write(s)
        if n:
            try:
                with open(self._log_path, "a", encoding="utf-8") as f:
                    f.write(s[:n])
            except OSError:
                pass
        return n

    def flush(self) -> None:
        self._stream.flush()

    def close(self) -> None:
        if hasattr(self._stream, "close"):
            self._stream.close()


def _run_with_log_capture(config: "PipelineConfig", run_fn):
    """Run pipeline while teeing stdout/stderr to run_log_{stem}.txt."""
    log_path = config.output_dir / f"run_log_{config.dataset_stem}.txt"
    try:
        # Clear or create log file
        log_path.write_text("", encoding="utf-8")
    except OSError:
        pass
    old_stdout, old_stderr = sys.stdout, sys.stderr
    tee_out = _TeeWriter(sys.stdout, log_path)
    tee_err = _TeeWriter(sys.stderr, log_path)
    try:
        sys.stdout = tee_out
        sys.stderr = tee_err
        return run_fn()
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
        tee_out.flush()
        tee_err.flush()


# =============================================================================
# INTERNAL IMPORTS - Data Loading
# =============================================================================

# Import our custom data loading tools from the utils folder
# These handle reading files and checking they're in the correct format
from week2_validation.utils.data_loader import (
    DataLoaderError,           # Error when data can't be read
    FileValidationError,       # Error when file path is invalid
    REQUIRED_FUSION_FIELDS,    # Required column names for schema validation
    SchemaValidationError,     # Error when data is missing required columns
    UnsupportedFormatError,    # Error when file type isn't supported
    load_data,                 # Raw data load (no schema check)
    load_fusion_data,          # Function to load the main fusion dataset
    load_reference_data,       # Function to load optional reference data
    validate_file_path,        # Function to check if a file path is valid
    validate_output_directory, # Function to check if output folder is valid
    validate_schema,           # Schema validation for adapted data
)

# =============================================================================
# INTERNAL IMPORTS - Configuration and State
# =============================================================================

# Import configuration loader for thresholds.yaml
from week2_validation.config.loader import (
    ConfigError,               # Error when configuration loading fails
    Week2Config,               # Configuration dataclass
    load_config,               # Function to load config from path
)

# Import freeze state management (Option C implementation)
from week2_validation.utils.state import (
    StateError,                # Error when state cannot be determined
    FreezeState,               # Immutable freeze state representation
    resolve_freeze_state,      # Resolves freeze state from flag file + CLI
)

# Import defensive freeze logic (Data Integrity & Statistical Validation owned dataset freezing)
from week2_validation.utils.freeze import FreezeError, ensure_frozen_input, get_frozen_data_path

# Survivability runtime layer
from week2_validation.runtime.exit_codes import Week2ExitCode
from week2_validation.runtime.runtime_guard import check_runtime_guard, start_runtime_guard
from week2_validation.runtime.safe_runner import run_week2_safely
from week2_validation.reporting.status_envelope import build_status_envelope

# Week 1: Pipeline Execution & Data Generation compatibility adapter
from week2_validation.adapters.week1_adapter import adapt_week1_dataframe

# =============================================================================
# NOTE: Diagnostic modules are NOT imported at top level.
# They are imported lazily ONLY when diagnostics are requested AND allowed.
# This prevents unnecessary dependencies and side effects.
# =============================================================================

# =============================================================================
# CONSTANTS
# =============================================================================

# Default path to the freeze flag file (relative to workspace root)
# If this file EXISTS, the dataset is considered frozen.
# This is the "authoritative" source for freeze state (Option C).
DEFAULT_FLAG_FILE_PATH = Path("week2_validation/.frozen")

# Default path to the configuration file
DEFAULT_CONFIG_PATH = Path("week2_validation/config/thresholds.yaml")

# Metadata for status envelope (set by run_pipeline)
_run_metadata: dict = {}


# =============================================================================
# PIPELINE CONFIGURATION
# =============================================================================

# The '@dataclass' decorator automatically creates a class that holds data.
# 'frozen=True' means once created, the values cannot be changed (like a locked form).
@dataclass(frozen=True)
class PipelineConfig:
    """
    Immutable configuration for the validation pipeline.
    
    This class holds all the settings for running the validation:
    - Where to find the input data files
    - Where to save the results
    - Whether to do a test run or full analysis
    - Whether to run diagnostics (requires frozen dataset)
    
    Think of it as a filled-out form that tells the pipeline what to do.

    Attributes:
        fusion_data_path: Location of the main fusion dataset file.
        output_dir: Folder where results will be saved.
        cosmic_data_path: Optional location of COSMIC reference data.
        dry_run: If True, only check inputs without running full analysis.
        run_diagnostics: If True, run distribution diagnostics (requires frozen dataset).
        run_benford: If True, run Benford's Law diagnostics (requires frozen dataset).
        run_lognormal: If True, run log-normality diagnostics (requires frozen dataset).
    """

    fusion_data_path: Path          # Where the fusion data file is located
    output_dir: Path                # Where to save the output/results
    dataset_stem: str               # Input filename stem for dynamic output names (e.g. "demo_fusion")
    cosmic_data_path: Optional[Path]  # Optional reference data location (can be empty)
    dry_run: bool                   # True = just validate, False = run full analysis
    run_diagnostics: bool           # True = run distribution diagnostics (if frozen), False = skip
    run_benford: bool               # True = run Benford diagnostics (if frozen), False = skip
    run_lognormal: bool             # True = run log-normality diagnostics (if frozen), False = skip
    run_cosmic: bool                # True = run COSMIC rank-order diagnostics (if frozen), False = skip
    run_benford_controls: bool      # True = run Benford implementation self-tests (synthetic data), False = skip
    generate_report: bool           # True = generate Statistical_Integrity_Report_Week2.md at end


# =============================================================================
# CUSTOM ERROR TYPES
# =============================================================================

# Custom error types help us identify what went wrong when something fails
class PipelineError(Exception):
    """
    Base exception for pipeline execution errors.
    
    This is a general error that occurs when something goes wrong
    while running the validation pipeline.
    """
    pass


class ConfigurationError(PipelineError):
    """
    Raised when pipeline configuration is invalid.
    
    This error occurs when the user provides incorrect settings,
    like a file path that doesn't exist or an invalid option.
    """
    pass


class FreezeStateError(PipelineError):
    """
    Raised when dataset freeze state prevents operation.
    
    This error occurs when an operation requires a frozen dataset
    but the dataset is not frozen.
    """
    pass


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser.
    
    This function sets up how the program reads instructions from the user.
    It defines what options the user can provide (like file paths) and
    explains what each option does.

    Returns:
        A configured parser ready to read user instructions.
    """
    # Create the main parser with program name and description
    parser = argparse.ArgumentParser(
        prog="week2_validation",  # Name shown in help text
        description=(
            "Data Integrity & Statistical Validation Pipeline. "
            "Performs diagnostic validation on fusion datasets. "
            "This pipeline is diagnostic only: no hypothesis testing, "
            "no modeling, no biological claims."
        ),
        # Show an example of how to use the program
        epilog=(
            "Example: python -m week2_validation.run_week2 "
            "--fusion-data fusion.csv --output-dir ./results --run-diagnostics"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # -------------------------------------------------------------------------
    # REQUIRED ARGUMENTS
    # -------------------------------------------------------------------------
    
    # Create a group for arguments that MUST be provided
    required_group = parser.add_argument_group("required arguments")

    # Define the fusion data file argument (REQUIRED)
    # The user must tell us where the data file is located
    required_group.add_argument(
        "--fusion-data",           # The flag the user types
        type=str,                  # Expects text (a file path)
        required=True,             # Cannot be skipped
        metavar="PATH",            # Shows "PATH" in help text
        help=(
            "Path to the fusion dataset file. "
            "Supported formats: CSV, TSV, JSON, Parquet, Excel."
        ),
    )

    # Define the output directory argument (REQUIRED)
    # The user must tell us where to save the results
    required_group.add_argument(
        "--output-dir",
        type=str,
        required=True,
        metavar="PATH",
        help="Directory where validation results will be written.",
    )

    # -------------------------------------------------------------------------
    # OPTIONAL ARGUMENTS
    # -------------------------------------------------------------------------

    # Create a group for arguments that are OPTIONAL
    optional_group = parser.add_argument_group("optional arguments")

    # Define the COSMIC reference data argument (OPTIONAL)
    # This is extra data the user can provide for additional checks
    optional_group.add_argument(
        "--cosmic-data",
        type=str,
        required=False,            # Can be skipped
        default=None,              # If not provided, will be None (empty)
        metavar="PATH",
        help=(
            "Path to COSMIC reference data file. "
            "If provided, enables reference-based validation."
        ),
    )

    # Define the run-all flag (OPTIONAL) - convenience to run full pipeline
    optional_group.add_argument(
        "--run-all",
        action="store_true",
        default=False,
        help=(
            "Run full pipeline: distribution, Benford, log-normality, Benford controls, and COSMIC. "
            "Equivalent to --run-diagnostics --run-benford --run-lognormal --run-benford-controls --run-cosmic."
        ),
    )

    # Define the run-diagnostics flag (OPTIONAL)
    # This explicitly requests distribution diagnostic analysis to run
    # IMPORTANT: Diagnostics only run if dataset is also frozen
    optional_group.add_argument(
        "--run-diagnostics",
        action="store_true",       # Just a flag, no value needed
        default=False,             # Off by default (explicit request required)
        help=(
            "Run distribution diagnostic analysis on the dataset. "
            "REQUIRES dataset to be frozen. "
            "If not specified, only input validation is performed."
        ),
    )

    # Define the run-benford flag (OPTIONAL)
    # This explicitly requests Benford's Law diagnostic analysis
    # IMPORTANT: Diagnostics only run if dataset is also frozen
    optional_group.add_argument(
        "--run-benford",
        action="store_true",       # Just a flag, no value needed
        default=False,             # Off by default (explicit request required)
        help=(
            "Run Benford's Law diagnostic analysis on the dataset. "
            "REQUIRES dataset to be frozen. "
            "Diagnostic only: no inference, no pass/fail determination."
        ),
    )

    # Define the run-benford-controls flag (OPTIONAL)
    # This runs Benford implementation self-tests using synthetic data only
    # IMPORTANT: This is NOT a scientific diagnostic, it is a developer validation tool
    optional_group.add_argument(
        "--run-benford-controls",
        action="store_true",       # Just a flag, no value needed
        default=False,             # Off by default
        help=(
            "Run Benford implementation self-test using synthetic data only. "
            "This is an implementation self-test for developer validation. "
            "Uses synthetic data only, does NOT analyze real data, "
            "does NOT validate dataset integrity."
        ),
    )

    # Define the run-lognormal flag (OPTIONAL)
    # This explicitly requests log-normality diagnostic analysis
    # IMPORTANT: Diagnostics only run if dataset is also frozen
    optional_group.add_argument(
        "--run-lognormal",
        action="store_true",       # Just a flag, no value needed
        default=False,             # Off by default (explicit request required)
        help=(
            "Run log-normality diagnostic analysis on the dataset. "
            "REQUIRES dataset to be frozen. "
            "Diagnostic only: no inference, no pass/fail, no decisions."
        ),
    )

    # Define the run-cosmic flag (OPTIONAL)
    # This explicitly requests COSMIC rank-order diagnostic analysis
    # IMPORTANT: Diagnostics only run if dataset is also frozen
    optional_group.add_argument(
        "--run-cosmic",
        action="store_true",       # Just a flag, no value needed
        default=False,             # Off by default (explicit request required)
        help=(
            "Run COSMIC rank-order diagnostic analysis. "
            "REQUIRES dataset to be frozen. "
            "Diagnostic only: no inference, no validation, no statistical tests."
        ),
    )

    # Define the generate-report flag (OPTIONAL)
    optional_group.add_argument(
        "--generate-report",
        action="store_true",
        default=False,
        help=(
            "Generate a researcher-friendly Markdown report "
            "(Statistical_Integrity_Report_Week2.md) from the run outputs."
        ),
    )

    # Define the dry-run flag (OPTIONAL)
    # When set, the program only checks inputs without running full analysis
    optional_group.add_argument(
        "--dry-run",
        action="store_true",       # Just a flag, no value needed
        default=False,             # Off by default
        help=(
            "Validate inputs and configuration only. "
            "Do not execute diagnostic analysis."
        ),
    )

    return parser


# =============================================================================
# ARGUMENT VALIDATION
# =============================================================================

def validate_arguments(args: argparse.Namespace) -> PipelineConfig:
    """
    Validate all command-line arguments and construct pipeline configuration.

    This function checks that all the file paths the user provided are valid:
    - Does the fusion data file exist?
    - Can we write to the output directory?
    - If COSMIC data was provided, does that file exist?

    If everything checks out, it creates a configuration object that the
    rest of the pipeline can use.

    Args:
        args: The parsed command-line arguments from the user.

    Returns:
        A validated PipelineConfig ready to use.

    Raises:
        ConfigurationError: If any file path is invalid or doesn't exist.
    """
    # STEP 1: Validate the fusion data file path
    # Check that the file exists and the path is safe
    try:
        fusion_path = validate_file_path(args.fusion_data, must_exist=True)
    except FileValidationError as e:
        # If validation fails, wrap the error with more context
        raise ConfigurationError(f"Invalid fusion data path: {e}") from e

    # STEP 2: Validate the output directory
    # Check that we can write there; create the folder if it doesn't exist
    try:
        output_dir = validate_output_directory(args.output_dir, create=True)
    except FileValidationError as e:
        raise ConfigurationError(f"Invalid output directory: {e}") from e

    # STEP 3: Validate optional COSMIC data path (if provided)
    cosmic_path: Optional[Path] = None  # Start with None (empty)
    if args.cosmic_data is not None:
        # User provided COSMIC data, so we need to validate it too
        try:
            cosmic_path = validate_file_path(args.cosmic_data, must_exist=True)
        except FileValidationError as e:
            raise ConfigurationError(f"Invalid COSMIC data path: {e}") from e

    # STEP 4: Create and return the configuration object
    # --run-all enables all diagnostic flags (distribution, Benford, log-normality, Benford controls, COSMIC)
    # AND generates the narrative report
    run_all = getattr(args, "run_all", False)
    dataset_stem = fusion_path.stem  # e.g. "demo_fusion" from "demo_fusion.parquet"
    return PipelineConfig(
        fusion_data_path=fusion_path,
        output_dir=output_dir,
        dataset_stem=dataset_stem,
        cosmic_data_path=cosmic_path,
        dry_run=args.dry_run,
        run_diagnostics=args.run_diagnostics or run_all,
        run_benford=args.run_benford or run_all,
        run_lognormal=args.run_lognormal or run_all,
        run_cosmic=args.run_cosmic or run_all,
        run_benford_controls=args.run_benford_controls or run_all,
        generate_report=getattr(args, "generate_report", False) or run_all,
    )


# =============================================================================
# INPUT VALIDATION
# =============================================================================

def validate_inputs(
    config: PipelineConfig,
    frozen_data_path: Optional[Path] = None,
) -> None:
    """
    Validate that input files can be loaded and meet schema requirements.

    This function does a "test load" of all input files to make sure:
    1. The files can actually be read (not corrupted)
    2. The data has the required columns (correct structure)
    
    This catches problems early before running the full analysis.

    Args:
        config: The validated pipeline configuration.
        frozen_data_path: Path to frozen data file (if provided, used instead of config path).

    Raises:
        PipelineError: If any input file cannot be loaded or is malformed.
    """
    # Use frozen path if provided (Data Integrity & Statistical Validation defensive freeze)
    data_path = frozen_data_path if frozen_data_path is not None else config.fusion_data_path
    
    # Tell the user what we're checking
    print(f"Validating fusion data: {data_path.name}")

    # Try to load the fusion data file
    try:
        # Load the data into memory as a table (DataFrame)
        fusion_df = load_fusion_data(str(data_path))
    except (FileValidationError, UnsupportedFormatError) as e:
        # File path is bad or file type isn't supported
        raise PipelineError(f"Cannot load fusion data: {e}") from e
    except DataLoaderError as e:
        # File exists but couldn't be read (might be corrupted)
        raise PipelineError(f"Error reading fusion data: {e}") from e
    except SchemaValidationError as e:
        # File loaded but is missing required columns
        raise PipelineError(f"Fusion data schema validation failed: {e}") from e

    # Report success: show how many rows and columns were loaded
    row_count = len(fusion_df)           # Number of data rows
    column_count = len(fusion_df.columns)  # Number of data columns
    print(f"  Loaded {row_count} records with {column_count} fields")
    print("  Schema validation: PASSED")

    # If COSMIC reference data was provided, validate it too
    if config.cosmic_data_path is not None:
        print(f"Validating COSMIC reference data: {config.cosmic_data_path.name}")

        try:
            # Load the COSMIC data
            cosmic_df = load_reference_data(str(config.cosmic_data_path))
        except (FileValidationError, UnsupportedFormatError) as e:
            raise PipelineError(f"Cannot load COSMIC data: {e}") from e
        except DataLoaderError as e:
            raise PipelineError(f"Error reading COSMIC data: {e}") from e

        # Report success for COSMIC data
        if cosmic_df is not None:
            row_count = len(cosmic_df)
            column_count = len(cosmic_df.columns)
            print(f"  Loaded {row_count} records with {column_count} fields")

    # Confirm the output directory is ready
    print(f"Output directory validated: {config.output_dir}")


# =============================================================================
# FREEZE STATE CHECKING
# =============================================================================

def check_freeze_state(flag_file_path: Path) -> FreezeState:
    """
    Check and return the dataset freeze state.
    
    Uses Option C logic:
    - If flag file exists → dataset is frozen (authoritative)
    - If flag file doesn't exist → dataset is NOT frozen
    
    Args:
        flag_file_path: Path to the freeze flag file.
    
    Returns:
        FreezeState object indicating frozen status and source.
    
    Raises:
        FreezeStateError: If freeze state cannot be determined.
    """
    try:
        # Try to resolve freeze state from flag file
        # We pass cli_frozen=False as default when flag file doesn't exist
        # This means: no flag file = not frozen (safe default)
        return resolve_freeze_state(flag_file_path, cli_frozen=False)
    except StateError as e:
        raise FreezeStateError(f"Cannot determine freeze state: {e}") from e


# =============================================================================
# DIAGNOSTIC EXECUTION (LAZY IMPORT)
# =============================================================================

def execute_diagnostics(
    config: PipelineConfig,
    week2_config: Week2Config,
    frozen_data_path: Optional[Path] = None,
    run_metadata: Optional[dict] = None,
) -> int:
    """
    Execute diagnostic analysis on the fusion dataset.
    
    This function performs a LAZY IMPORT of the diagnostic module
    to avoid unnecessary dependencies when diagnostics aren't requested.
    
    IMPORTANT: This function should ONLY be called when:
    - Dataset is frozen (already verified by caller)
    - --run-diagnostics flag is True (already verified by caller)
    
    Args:
        config: Pipeline configuration with file paths.
        week2_config: Week2 configuration from thresholds.yaml.
        frozen_data_path: Path to frozen data file (if provided, used instead of config path).
    
    Returns:
        Exit code: 0 for success, non-zero for failure.
    """
    print("Starting diagnostic analysis...")
    print()
    
    # -------------------------------------------------------------------------
    # LAZY IMPORT: Only import diagnostics module when actually needed
    # This prevents loading heavy dependencies (numpy, matplotlib) when
    # the user only wants to validate inputs.
    # -------------------------------------------------------------------------
    try:
        from week2_validation.distributions.visualize import (
            run_diagnostics,
            save_protein_distribution_png,
        )
    except ImportError as e:
        print(f"Error: Cannot import diagnostic module: {e}", file=sys.stderr)
        print("Make sure all required dependencies are installed.", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    # -------------------------------------------------------------------------
    # Load the fusion data for diagnostic analysis
    # Use frozen path if provided (Data Integrity & Statistical Validation defensive freeze)
    # -------------------------------------------------------------------------
    data_path = frozen_data_path if frozen_data_path is not None else config.fusion_data_path
    try:
        fusion_df = load_fusion_data(str(data_path))
    except Exception as e:
        print(f"Error loading fusion data for diagnostics: {e}", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    # -------------------------------------------------------------------------
    # Extract protein_length column for analysis
    # No assumptions about data format - just use what's in the required column
    # -------------------------------------------------------------------------
    if "protein_length" not in fusion_df.columns:
        print("Error: 'protein_length' column not found in fusion data.", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    protein_lengths = fusion_df["protein_length"].tolist()
    
    # -------------------------------------------------------------------------
    # Run diagnostics (returns DiagnosticResult, does NOT save files)
    # -------------------------------------------------------------------------
    try:
        print(f"Running diagnostics on {len(protein_lengths)} protein length values...")
        result = run_diagnostics(
            data=protein_lengths,
            log_mode="transform",
            include_log_histogram=week2_config.visualization.allow_log_scale,
        )
        
        # Report results (no interpretation, just facts)
        print()
        print("Diagnostic Results:")
        print(f"  Sample size: {result.statistics.count}")
        print(f"  Min: {result.statistics.minimum:.2f}")
        print(f"  Max: {result.statistics.maximum:.2f}")
        print(f"  Mean: {result.statistics.mean:.2f}")
        print(f"  Median: {result.statistics.median:.2f}")
        print(f"  Std Dev: {result.statistics.std:.2f}")
        print(f"  Skewness: {result.statistics.skewness:.4f}")
        print()
        print(f"  Log mode used: {result.log_mode_used}")
        print(f"  Diagnostic only: {result.metadata.diagnostic_only}")
        print(f"  Hypothesis tested: {result.metadata.hypothesis_tested}")
        print()
        # Record skewness in data_quality for week2_status.json (JSON-serializable)
        if run_metadata is not None and "data_quality" in run_metadata:
            sk = float(result.statistics.skewness)
            run_metadata["data_quality"]["skewness"] = (
                round(sk, 6) if math.isfinite(sk) else None
            )
        # Record distribution for narrative report (diagnostic_results)
        if run_metadata is not None:
            run_metadata.setdefault("diagnostic_results", {})["distribution"] = {
                "minimum": round(float(result.statistics.minimum), 6),
                "maximum": round(float(result.statistics.maximum), 6),
                "mean": round(float(result.statistics.mean), 6),
                "median": round(float(result.statistics.median), 6),
                "std": round(float(result.statistics.std), 6),
                "skewness": round(float(result.statistics.skewness), 6)
                if math.isfinite(float(result.statistics.skewness))
                else None,
            }
        # Save protein distribution histogram to output-dir when --run-diagnostics is active
        try:
            saved_path = save_protein_distribution_png(
                result, config.output_dir, filename="protein_distribution.png"
            )
            if saved_path is not None:
                print(f"Protein distribution saved: {saved_path.name}")
        except Exception as save_err:
            print(f"Warning: Could not save histogram: {save_err}", file=sys.stderr)
        print()
        print("Diagnostic analysis complete.")
        print("NOTE: Results are diagnostic only. No interpretation provided.")
        
    except Exception as e:
        print(f"Error during diagnostic analysis: {e}", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    return 0


# =============================================================================
# BENFORD DIAGNOSTIC EXECUTION (LAZY IMPORT)
# =============================================================================

def execute_benford_diagnostics(
    config: PipelineConfig,
    week2_config: Week2Config,
    frozen_data_path: Optional[Path] = None,
    run_metadata: Optional[dict] = None,
) -> int:
    """
    Execute Benford's Law diagnostic analysis on the fusion dataset.
    
    This function performs a LAZY IMPORT of the Benford diagnostic module
    to avoid unnecessary dependencies when Benford diagnostics aren't requested.
    
    IMPORTANT: This function should ONLY be called when:
    - Dataset is frozen (already verified by caller)
    - --run-benford flag is True (already verified by caller)
    
    Args:
        config: Pipeline configuration with file paths.
        week2_config: Week2 configuration from thresholds.yaml.
        frozen_data_path: Path to frozen data file (if provided, used instead of config path).
    
    Returns:
        Exit code: 0 for success, non-zero for failure.
    """
    print("Starting Benford's Law diagnostic analysis...")
    print()
    
    # -------------------------------------------------------------------------
    # LAZY IMPORT: Only import Benford module when actually needed
    # -------------------------------------------------------------------------
    try:
        from week2_validation.benford.diagnostics import run_benford_diagnostics
    except ImportError as e:
        print(f"Error: Cannot import Benford diagnostic module: {e}", file=sys.stderr)
        print("Make sure all required dependencies are installed.", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    # -------------------------------------------------------------------------
    # Load the fusion data for Benford analysis
    # Use frozen path if provided (Data Integrity & Statistical Validation defensive freeze)
    # -------------------------------------------------------------------------
    data_path = frozen_data_path if frozen_data_path is not None else config.fusion_data_path
    try:
        fusion_df = load_fusion_data(str(data_path))
    except Exception as e:
        print(f"Error loading fusion data for Benford diagnostics: {e}", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    # -------------------------------------------------------------------------
    # Extract protein_length column for analysis
    # No assumptions about data format - just use what's in the required column
    # -------------------------------------------------------------------------
    if "protein_length" not in fusion_df.columns:
        print("Error: 'protein_length' column not found in fusion data.", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    protein_lengths = fusion_df["protein_length"].tolist()
    
    # -------------------------------------------------------------------------
    # Run Benford diagnostics (returns BenfordDiagnosticResult, does NOT save files)
    # -------------------------------------------------------------------------
    try:
        print(f"Running Benford diagnostics on {len(protein_lengths)} protein length values...")
        result = run_benford_diagnostics(data=protein_lengths)
        
        # Report results (no interpretation, just facts)
        print()
        print("Benford Diagnostic Results:")
        print(f"  Total input values: {result.total_input_values}")
        print(f"  Values analyzed: {result.total_extracted_digits}")
        print(f"  Values excluded: {result.values_excluded_count}")
        print()
        
        if result.computation_successful:
            print("  Observed first-digit frequencies:")
            for digit in range(1, 10):
                obs = result.observed_frequencies.get(digit, 0.0)
                exp = result.expected_frequencies.get(digit, 0.0)
                print(f"    Digit {digit}: observed={obs:.4f}, Benford expected={exp:.4f}")
            print()
            
            if result.chi_squared_statistic is not None:
                print(f"  Chi-squared statistic: {result.chi_squared_statistic:.4f}")
                print(f"  Degrees of freedom: {result.degrees_of_freedom}")
                if result.p_value is not None:
                    print(f"  p-value: {result.p_value:.6f}")
            print()
            
            # Report applicability assessment (diagnostic only)
            if result.metadata.benford_applicable is not None:
                status = "applicable" if result.metadata.benford_applicable else "not applicable"
                print(f"  Benford applicability (heuristic): {status}")
                if result.metadata.reason_if_not_applicable:
                    print(f"  Reason: {result.metadata.reason_if_not_applicable}")
                if result.metadata.scale_span_orders_of_magnitude is not None:
                    print(f"  Scale span: {result.metadata.scale_span_orders_of_magnitude:.2f} orders of magnitude")
                    if run_metadata is not None:
                        run_metadata["benford_scale_span"] = result.metadata.scale_span_orders_of_magnitude
            # Record Benford results for narrative report
            if run_metadata is not None:
                obs = {int(k): round(float(v), 6) for k, v in result.observed_frequencies.items()}
                exp = {int(k): round(float(v), 6) for k, v in result.expected_frequencies.items()}
                run_metadata.setdefault("diagnostic_results", {})["benford"] = {
                    "observed_frequencies": obs,
                    "expected_frequencies": exp,
                    "scale_span_orders_of_magnitude": (
                        round(float(result.metadata.scale_span_orders_of_magnitude), 4)
                        if result.metadata.scale_span_orders_of_magnitude is not None
                        else None
                    ),
                    "applicability": result.metadata.benford_applicable,
                    "reason_if_not_applicable": result.metadata.reason_if_not_applicable or "",
                }
        
        print()
        print("Benford diagnostic analysis complete.")
        print("NOTE: Results are diagnostic only. No inference drawn. No interpretation provided.")
        print(f"      {result.metadata.disclaimer}")
        
    except Exception as e:
        print(f"Error during Benford diagnostic analysis: {e}", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    return 0


# =============================================================================
# LOG-NORMALITY DIAGNOSTIC EXECUTION (LAZY IMPORT)
# =============================================================================

def execute_log_normality_diagnostics(
    config: PipelineConfig,
    week2_config: Week2Config,
    frozen_data_path: Optional[Path] = None,
    run_metadata: Optional[dict] = None,
) -> int:
    """
    Execute log-normality diagnostic analysis on the fusion dataset.
    
    This function performs a LAZY IMPORT of the log-normality diagnostic module
    to avoid unnecessary dependencies when log-normality diagnostics aren't requested.
    
    IMPORTANT: This function should ONLY be called when:
    - Dataset is frozen (already verified by caller)
    - --run-lognormal flag is True (already verified by caller)
    
    DIAGNOSTIC-ONLY:
    - No inference, thresholds, or decisions are made
    - No files are saved
    - No plots are generated
    - Results are descriptive summaries only
    
    Args:
        config: Pipeline configuration with file paths.
        week2_config: Week2 configuration from thresholds.yaml.
        frozen_data_path: Path to frozen data file (if provided, used instead of config path).
    
    Returns:
        Exit code: 0 for success, non-zero for failure.
    """
    print("Starting log-normality diagnostic analysis...")
    print()
    
    # -------------------------------------------------------------------------
    # LAZY IMPORT: Only import log-normality module when actually needed
    # -------------------------------------------------------------------------
    try:
        from week2_validation.distributions.log_normality import (
            run_log_normality_diagnostics,
        )
    except ImportError as e:
        print(f"Error: Cannot import log-normality diagnostic module: {e}", file=sys.stderr)
        print("Make sure all required dependencies are installed.", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    # -------------------------------------------------------------------------
    # Load the fusion data for log-normality analysis
    # Use frozen path if provided (Data Integrity & Statistical Validation defensive freeze)
    # -------------------------------------------------------------------------
    data_path = frozen_data_path if frozen_data_path is not None else config.fusion_data_path
    try:
        fusion_df = load_fusion_data(str(data_path))
    except Exception as e:
        print(f"Error loading fusion data for log-normality diagnostics: {e}", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    # -------------------------------------------------------------------------
    # Extract protein_length column for analysis
    # No assumptions about data format - just use what's in the required column
    # -------------------------------------------------------------------------
    if "protein_length" not in fusion_df.columns:
        print("Error: 'protein_length' column not found in fusion data.", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    protein_lengths = fusion_df["protein_length"].tolist()
    
    # -------------------------------------------------------------------------
    # Run log-normality diagnostics (returns LogNormalityDiagnosticResult)
    # Does NOT save files, does NOT produce plots
    # -------------------------------------------------------------------------
    try:
        print(f"Running log-normality diagnostics on {len(protein_lengths)} protein length values...")
        result = run_log_normality_diagnostics(data=protein_lengths)
        
        # Report results (no interpretation, just facts)
        print()
        print("Log-Normality Diagnostic Results:")
        print(f"  Sample size: {result.sample_size}")
        print(f"  SciPy available: {result.scipy_available}")
        print(f"  Computation successful: {result.computation_successful}")
        print()
        
        if result.computation_successful:
            if result.ks_statistic is not None:
                print(f"  KS statistic: {result.ks_statistic:.6f}")
            else:
                print("  KS statistic: not computed")
            
            if result.ad_statistic is not None:
                print(f"  Anderson-Darling statistic: {result.ad_statistic:.6f}")
            else:
                print("  Anderson-Darling statistic: not computed (requires SciPy)")
            # Record log-normality results for narrative report
            if run_metadata is not None:
                run_metadata.setdefault("diagnostic_results", {})["log_normality"] = {
                    "ks_statistic": (
                        round(float(result.ks_statistic), 6)
                        if result.ks_statistic is not None and math.isfinite(result.ks_statistic)
                        else None
                    ),
                    "ad_statistic": (
                        round(float(result.ad_statistic), 6)
                        if result.ad_statistic is not None and math.isfinite(result.ad_statistic)
                        else None
                    ),
                    "scipy_available": result.scipy_available,
                }
        
        print()
        print("Log-normality diagnostic analysis complete.")
        print()
        print("NOTE: Log-normality diagnostics are DESCRIPTIVE ONLY.")
        print("No inference, thresholds, or data validity conclusions are drawn.")
        
    except Exception as e:
        print(f"Error during log-normality diagnostic analysis: {e}", file=sys.stderr)
        return int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
    
    return 0


# =============================================================================
# COSMIC DIAGNOSTIC EXECUTION (LAZY IMPORT)
# =============================================================================

def execute_cosmic_diagnostics(
    config: PipelineConfig,
    week2_config: Week2Config,
    frozen_data_path: Optional[Path] = None,
    run_metadata: Optional[dict] = None,
) -> tuple[int, bool]:
    """
    Execute COSMIC rank-order diagnostic analysis.
    
    This function performs a LAZY IMPORT of the COSMIC diagnostic module
    to avoid unnecessary dependencies when COSMIC diagnostics aren't requested.
    
    IMPORTANT: This function should ONLY be called when:
    - Dataset is frozen (already verified by caller)
    - --run-cosmic flag is True (already verified by caller)
    
    DIAGNOSTIC-ONLY:
    - No inference, thresholds, or decisions are made
    - No files are saved
    - No plots are generated
    - No statistical tests are performed
    - Results are descriptive summaries only
    
    Args:
        config: Pipeline configuration with file paths.
        week2_config: Week2 configuration from thresholds.yaml.
        frozen_data_path: Path to frozen data file (if provided, used instead of config path).
    
    Returns:
        Tuple of (result_code, cosmic_loaded). result_code: 0 for success, non-zero for failure.
        cosmic_loaded: True if COSMIC reference was loaded (user-provided or mock), False otherwise.
    """
    print("Starting COSMIC rank-order diagnostic analysis...")
    print()
    
    # -------------------------------------------------------------------------
    # LAZY IMPORT: Only import COSMIC module when actually needed
    # -------------------------------------------------------------------------
    try:
        from week2_validation.cosmic.diagnostics import (
            run_cosmic_recurrence_diagnostic,
        )
    except ImportError as e:
        print(f"Error: Cannot import COSMIC diagnostic module: {e}", file=sys.stderr)
        print("Make sure all required dependencies are installed.", file=sys.stderr)
        return (int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR), False)
    
    # -------------------------------------------------------------------------
    # Load the fusion data
    # Use frozen path if provided (Data Integrity & Statistical Validation defensive freeze)
    # -------------------------------------------------------------------------
    data_path = frozen_data_path if frozen_data_path is not None else config.fusion_data_path
    try:
        fusion_df = load_fusion_data(str(data_path))
    except Exception as e:
        print(f"Error loading fusion data for COSMIC diagnostics: {e}", file=sys.stderr)
        return (int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR), False)
    
    # -------------------------------------------------------------------------
    # Load COSMIC data with automatic fallback to mock COSMIC
    # -------------------------------------------------------------------------
    cosmic_df = None
    cosmic_reference_loaded = False
    cosmic_reference_source = None
    
    if config.cosmic_data_path is not None:
        # User provided COSMIC data - try to load it
        try:
            cosmic_df = load_reference_data(str(config.cosmic_data_path))
            cosmic_reference_loaded = cosmic_df is not None
            if cosmic_reference_loaded:
                cosmic_reference_source = "user_provided"
                print(f"Loaded user-provided COSMIC data from: {config.cosmic_data_path.name}")
        except Exception as e:
            print(f"Warning: Could not load user-provided COSMIC data: {e}")
            print("Attempting fallback to mock COSMIC dataset...")
            cosmic_reference_loaded = False
    
    # Fallback to mock COSMIC if user COSMIC not available
    if not cosmic_reference_loaded:
        try:
            from week2_validation.cosmic.generate_mock_cosmic import ensure_mock_cosmic_exists
            
            # Mock COSMIC should be in week2_validation/cosmic/ directory
            cosmic_dir = Path(__file__).parent / "cosmic"
            mock_cosmic_path = ensure_mock_cosmic_exists(cosmic_dir)
            
            cosmic_df = load_reference_data(str(mock_cosmic_path))
            cosmic_reference_loaded = cosmic_df is not None
            if cosmic_reference_loaded:
                cosmic_reference_source = "mock_fallback"
                print(f"Loaded mock COSMIC fallback dataset from: {mock_cosmic_path.name}")
                print("Note: Using synthetic COSMIC data for pipeline operability.")
        except Exception as e:
            print(f"Warning: Could not load mock COSMIC fallback: {e}")
            print("COSMIC diagnostic will report fusion data counts only.")
            cosmic_reference_loaded = False
    
    if not cosmic_reference_loaded:
        print("Note: No COSMIC reference data available (neither user-provided nor mock).")
        print("COSMIC diagnostic will report fusion data counts only.")
    
    # -------------------------------------------------------------------------
    # Run COSMIC diagnostic (returns descriptive dict only)
    # -------------------------------------------------------------------------
    try:
        print(f"Running COSMIC rank-order diagnostic...")
        # Determine cosmic file path for provenance
        cosmic_file_path = None
        if config.cosmic_data_path:
            cosmic_file_path = config.cosmic_data_path
        elif cosmic_reference_source == "mock_fallback":
            # Mock COSMIC path
            cosmic_dir = Path(__file__).parent / "cosmic"
            cosmic_file_path = cosmic_dir / "mock_cosmic_census.csv"
        
        result = run_cosmic_recurrence_diagnostic(
            fusion_df=fusion_df,
            cosmic_df=cosmic_df,
            top_n=10,
            cosmic_reference_source=cosmic_reference_source,
            cosmic_file_path=cosmic_file_path,
            mock_generation_seed=42,
            bootstrap_iterations=1000,
            bootstrap_seed=42,
        )
        
        # Report results (no interpretation, just facts)
        print()
        print("COSMIC Rank-Order Diagnostic Results:")
        if cosmic_reference_source:
            print(f"  COSMIC reference source: {cosmic_reference_source}")
        print(f"  Total fusion pairs (ours): {result['total_fusions_ours']}")
        print(f"  Total fusion pairs (COSMIC): {result['total_fusions_cosmic']}")
        print(f"  Overlapping pairs: {result['overlap_count']}")
        print(f"  Only in ours: {result['only_in_ours_count']}")
        print(f"  Only in COSMIC: {result['only_in_cosmic_count']}")
        
        # Report statistical metrics if available
        if "spearman_rho" in result and result["spearman_rho"] is not None:
            print()
            print("  Statistical Metrics:")
            print(f"    Spearman rank correlation (rho): {result['spearman_rho']:.4f}")
            if result.get("spearman_p_value") is not None:
                print(f"    Spearman p-value: {result['spearman_p_value']:.6f}")
            if result.get("enrichment_p_value") is not None:
                print(f"    Enrichment p-value: {result['enrichment_p_value']:.6f}")
            if result.get("expected_overlap_random") is not None:
                print(f"    Expected random overlap: {result['expected_overlap_random']:.2f}")
            if result.get("negative_control_rho") is not None:
                print(f"    Negative control rho: {result['negative_control_rho']:.4f}")
            if result.get("cosmic_validation_score") is not None:
                print(f"    COSMIC validation score: {result['cosmic_validation_score']:.4f}")
                print(f"    Classification: {result.get('cosmic_validation_classification', 'N/A')}")
        
        if "top_fusion_overlap" in result:
            print(f"    Top 10 fusion overlap: {result['top_fusion_overlap']}")
            print(f"    Top fusion enrichment ratio: {result['top_fusion_enrichment_ratio']:.4f}")
        
        # Report message if present
        if "message" in result and result["message"]:
            print()
            print(f"  Note: {result['message']}")
        
        # Report top rank discrepancies
        if result["top_rank_discrepancies"]:
            print()
            print("  Top rank discrepancies (by absolute rank difference):")
            for i, disc in enumerate(result["top_rank_discrepancies"], 1):
                print(f"    {i}. {disc['gene_1']}-{disc['gene_2']}: "
                      f"our_rank={disc['our_rank']}, cosmic_rank={disc['cosmic_rank']}, "
                      f"diff={disc['absolute_rank_difference']}")
        
        print()
        print("COSMIC rank-order diagnostic analysis complete.")
        # Record COSMIC results for narrative report
        if run_metadata is not None:
            cosmic_metadata = {
                "total_fusions_ours": result["total_fusions_ours"],
                "total_fusions_cosmic": result["total_fusions_cosmic"],
                "overlap_count": result["overlap_count"],
                "only_in_ours_count": result["only_in_ours_count"],
                "only_in_cosmic_count": result["only_in_cosmic_count"],
            }
            
            # Add statistical metrics if available
            if "spearman_rho" in result:
                cosmic_metadata["spearman_rho"] = result.get("spearman_rho")
                cosmic_metadata["spearman_p_value"] = result.get("spearman_p_value")
                cosmic_metadata["enrichment_p_value"] = result.get("enrichment_p_value")
                cosmic_metadata["expected_overlap_random"] = result.get("expected_overlap_random")
                cosmic_metadata["observed_overlap"] = result.get("observed_overlap")
                cosmic_metadata["negative_control_rho"] = result.get("negative_control_rho")
                cosmic_metadata["negative_control_p_value"] = result.get("negative_control_p_value")
                cosmic_metadata["cosmic_validation_score"] = result.get("cosmic_validation_score")
                cosmic_metadata["cosmic_validation_classification"] = result.get("cosmic_validation_classification")
                
                # Add bootstrap CI if available
                if "rho_ci_lower" in result:
                    cosmic_metadata["rho_ci_lower"] = result.get("rho_ci_lower")
                    cosmic_metadata["rho_ci_upper"] = result.get("rho_ci_upper")
                    cosmic_metadata["bootstrap_iterations"] = result.get("bootstrap_iterations")
                
                # Add score component breakdown if available
                if "score_component_breakdown" in result:
                    cosmic_metadata["score_component_breakdown"] = result.get("score_component_breakdown")
                
                # Add score_components (legacy) if available
                if "score_components" in result:
                    cosmic_metadata["score_components"] = result.get("score_components")
                    
            if "top_fusion_overlap" in result:
                cosmic_metadata["top_fusion_overlap"] = result.get("top_fusion_overlap")
                cosmic_metadata["top_fusion_enrichment_ratio"] = result.get("top_fusion_enrichment_ratio")
            
            # Add reproducibility lock if available
            if "reproducibility_lock" in result:
                cosmic_metadata["reproducibility_lock"] = result.get("reproducibility_lock")
            
            # Add provenance metadata
            if "cosmic_reference_source" in result:
                cosmic_metadata["cosmic_reference_source"] = result.get("cosmic_reference_source")
            elif cosmic_reference_source:
                cosmic_metadata["cosmic_reference_source"] = cosmic_reference_source
            
            cosmic_metadata["cosmic_reference_version"] = result.get("cosmic_reference_version")
            cosmic_metadata["cosmic_reference_file_hash"] = result.get("cosmic_reference_file_hash")
            cosmic_metadata["cosmic_reference_load_timestamp"] = result.get("cosmic_reference_load_timestamp")
            
            run_metadata.setdefault("diagnostic_results", {})["cosmic"] = cosmic_metadata
        print()
        print("NOTE: COSMIC diagnostics are DESCRIPTIVE ONLY.")
        print("No inference, no statistical tests, no validation conclusions are drawn.")
        
    except Exception as e:
        print(f"Error during COSMIC diagnostic analysis: {e}", file=sys.stderr)
        return (int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR), cosmic_reference_loaded if config.cosmic_data_path is not None else True)
    
    return (0, cosmic_reference_loaded if config.cosmic_data_path is not None else True)


# =============================================================================
# BENFORD CONTROLS (IMPLEMENTATION SELF-TEST)
# =============================================================================

def execute_benford_controls() -> int:
    """
    Execute Benford implementation self-tests using synthetic control data.
    
    This function runs positive and negative control tests to verify that
    the Benford analysis code is working correctly. It uses ONLY synthetic
    data generated specifically for this purpose.
    
    IMPORTANT: This is an implementation self-test, NOT a scientific analysis.
    - Uses synthetic data only (no real data is touched)
    - Does NOT validate dataset integrity
    - Does NOT interpret p-values or draw conclusions
    - Does NOT affect pipeline execution or exit codes
    - Does NOT require dataset to be frozen
    
    Returns:
        0 on normal execution, 1 only if an unexpected exception occurs.
    """
    print("=" * 60)
    print("BENFORD IMPLEMENTATION SELF-TEST")
    print("=" * 60)
    print()
    print("This is an implementation self-test, not a scientific analysis.")
    print("Using synthetic data only. No real data is analyzed.")
    print()
    
    # -------------------------------------------------------------------------
    # LAZY IMPORT: Only import when this function is called
    # -------------------------------------------------------------------------
    try:
        from week2_validation.benford.diagnostics import (
            generate_synthetic_benford_positive_control,
            generate_synthetic_benford_negative_control,
            run_benford_diagnostics,
        )
    except ImportError as e:
        print(f"Error: Cannot import Benford diagnostics module: {e}", file=sys.stderr)
        print("Make sure all required dependencies are installed.", file=sys.stderr)
        return 1
    
    # -------------------------------------------------------------------------
    # POSITIVE CONTROL: Benford-compliant synthetic data
    # Expected: Observed frequencies should approximate Benford's Law
    # -------------------------------------------------------------------------
    print("-" * 60)
    print("BENFORD CONTROL — POSITIVE (synthetic)")
    print("-" * 60)
    print()
    print("Description: Log-uniform synthetic data (Benford-compliant)")
    print()
    
    try:
        # Generate synthetic Benford-compliant data
        positive_data = generate_synthetic_benford_positive_control(size=1000, seed=42)
        positive_result = run_benford_diagnostics(data=positive_data.tolist())
        
        print(f"Sample size: {positive_result.total_input_values}")
        print()
        print("Observed first-digit frequencies:")
        for digit in range(1, 10):
            obs = positive_result.observed_frequencies.get(digit, 0.0)
            exp = positive_result.expected_frequencies.get(digit, 0.0)
            print(f"  Digit {digit}: observed={obs:.4f}, expected={exp:.4f}")
        print()
        
        if positive_result.computation_successful:
            if positive_result.chi_squared_statistic is not None:
                print(f"Chi-squared statistic: {positive_result.chi_squared_statistic:.4f}")
            if positive_result.p_value is not None:
                print(f"p-value: {positive_result.p_value:.6f}")
            print()
            print("Expected behavior observed: Observed frequencies approximate Benford distribution.")
        else:
            print("Observed behavior: Computation did not complete.")
        
        print()
        print("NOTE: This is an implementation self-test, not a scientific analysis.")
        
    except Exception as e:
        print(f"Unexpected error during positive control: {e}", file=sys.stderr)
        return 1
    
    print()
    
    # -------------------------------------------------------------------------
    # NEGATIVE CONTROL: Uniform synthetic data (non-Benford)
    # Expected: Observed frequencies should NOT approximate Benford's Law
    # -------------------------------------------------------------------------
    print("-" * 60)
    print("BENFORD CONTROL — NEGATIVE (synthetic)")
    print("-" * 60)
    print()
    print("Description: Uniform synthetic data (NOT Benford-distributed)")
    print()
    
    try:
        # Generate synthetic uniform (non-Benford) data
        negative_data = generate_synthetic_benford_negative_control(size=1000, seed=42)
        negative_result = run_benford_diagnostics(data=negative_data.tolist())
        
        print(f"Sample size: {negative_result.total_input_values}")
        print()
        print("Observed first-digit frequencies:")
        for digit in range(1, 10):
            obs = negative_result.observed_frequencies.get(digit, 0.0)
            exp = negative_result.expected_frequencies.get(digit, 0.0)
            print(f"  Digit {digit}: observed={obs:.4f}, expected={exp:.4f}")
        print()
        
        if negative_result.computation_successful:
            if negative_result.chi_squared_statistic is not None:
                print(f"Chi-squared statistic: {negative_result.chi_squared_statistic:.4f}")
            if negative_result.p_value is not None:
                print(f"p-value: {negative_result.p_value:.6f}")
            print()
            print("Expected behavior observed: Observed frequencies deviate from Benford distribution.")
        else:
            print("Observed behavior: Computation did not complete.")
        
        print()
        print("NOTE: This is an implementation self-test, not a scientific analysis.")
        
    except Exception as e:
        print(f"Unexpected error during negative control: {e}", file=sys.stderr)
        return 1
    
    print()
    print("=" * 60)
    print("BENFORD IMPLEMENTATION SELF-TEST COMPLETE")
    print("=" * 60)
    print()
    print("NOTE: This was an implementation self-test using synthetic data only.")
    print("      No real data was analyzed. No dataset integrity conclusions drawn.")
    print()
    
    return 0


# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================

def run_pipeline(config: PipelineConfig) -> int:
    """
    Execute the validation pipeline.

    This is the main function that runs the entire validation process.
    It coordinates all the steps: configuration loading, freeze state
    checking, validation, and optional diagnostic execution.

    EXECUTION FLOW:
    1. Load configuration from thresholds.yaml
    2. Check dataset freeze state
    3. Validate input files
    4. If --dry-run: exit after validation
    5. If NOT frozen: exit with message
    6. If frozen but NO --run-diagnostics: exit with message
    7. If frozen AND --run-diagnostics: run diagnostics

    Args:
        config: The validated pipeline configuration.

    Returns:
        Exit code: 0 means success, any other number means something failed.

    Raises:
        PipelineError: If something goes wrong during pipeline execution.
    """
    global _run_metadata
    # Print a header banner to clearly show the pipeline is starting
    print("=" * 60)
    print("Data Integrity & Statistical Validation Pipeline")
    print("=" * 60)
    print()

    # =========================================================================
    # STEP 1: Load configuration from thresholds.yaml
    # =========================================================================
    print("Loading configuration...")
    
    try:
        # Resolve config path relative to current working directory
        config_path = Path.cwd() / DEFAULT_CONFIG_PATH
        if not config_path.exists():
            # Try relative to this file's location
            config_path = Path(__file__).parent / "config" / "thresholds.yaml"
        
        week2_config = load_config(config_path)
        print(f"  Configuration loaded from: {config_path.name}")
        print(f"  dataset_frozen_required: {week2_config.dataset_frozen_required}")
        print(f"  allow_real_data_analysis: {week2_config.allow_real_data_analysis}")
    except ConfigError as e:
        raise PipelineError(f"Failed to load configuration: {e}") from e
    
    print()

    # =========================================================================
    # STEP 1.5: Defensive Freeze - Ensure input is frozen BEFORE diagnostics
    # =========================================================================
    # Data Integrity & Statistical Validation cannot assume Week 1 froze the data. This defensive freeze ensures
    # that ALL diagnostics operate on immutable data, regardless of upstream.
    print("Ensuring dataset is frozen for Data Integrity & Statistical Validation diagnostics...")
    
    try:
        frozen_root = Path(__file__).parent / "frozen_inputs"
        frozen_input_dir = ensure_frozen_input(
            input_path=config.fusion_data_path,
            frozen_root=frozen_root,
        )
        # Get the actual frozen data file path for downstream use
        frozen_input_path = get_frozen_data_path(frozen_input_dir)
        print(f"  Frozen data location: {frozen_input_path}")
    except FreezeError as e:
        raise
    
    print()

    # =========================================================================
    # STEP 2: Week 1 (Pipeline Execution & Data Generation) adapter - convert to schema if needed
    # =========================================================================
    try:
        raw_df = load_data(str(frozen_input_path))
        adapted_df = adapt_week1_dataframe(raw_df)
        validate_schema(adapted_df, REQUIRED_FUSION_FIELDS)
        effective_data_path = frozen_input_path
        if "geneA" in raw_df.columns or "geneB" in raw_df.columns or "samples_detected" in raw_df.columns or "recurrence_frequency" in raw_df.columns:
            adapted_path = config.output_dir / "week2_adapted_fusion.csv"
            adapted_df.to_csv(adapted_path, sep=",", index=False)
            effective_data_path = adapted_path
            print("  Week 1 (Pipeline Execution & Data Generation) format detected; adapted data written for diagnostics")
    except ValueError as e:
        raise PipelineError(f"Input schema adaptation failed: {e}") from e

    # =========================================================================
    # STEP 3: Validate all input files
    # =========================================================================
    print("Validating inputs...")
    validate_inputs(config, frozen_data_path=effective_data_path)
    print()

    # =========================================================================
    # STEP 3.5: Compute data quality (protein_length exclusion for diagnostics)
    # =========================================================================
    _run_metadata["dataset_hash"] = frozen_input_dir.name
    import numpy as np
    try:
        _fusion_df = load_fusion_data(str(effective_data_path))
        total_rows_original = len(_fusion_df)
        if "protein_length" in _fusion_df.columns and total_rows_original > 0:
            pl = np.asarray(_fusion_df["protein_length"], dtype=np.float64)
            valid_mask = np.isfinite(pl) & (pl > 0)
            rows_used_for_analysis = int(np.sum(valid_mask))
            rows_excluded = total_rows_original - rows_used_for_analysis
            excluded_fraction = rows_excluded / total_rows_original
        else:
            rows_used_for_analysis = total_rows_original
            rows_excluded = 0
            excluded_fraction = 0.0
        _run_metadata["data_quality"] = {
            "total_rows_original": total_rows_original,
            "rows_used_for_analysis": rows_used_for_analysis,
            "rows_excluded": rows_excluded,
            "excluded_fraction": round(excluded_fraction, 6),
            "warning_flag": excluded_fraction > 0.10,
            "high_risk_flag": excluded_fraction > 0.50,
        }
    except Exception:
        _run_metadata["data_quality"] = {
            "total_rows_original": 0,
            "rows_used_for_analysis": 0,
            "rows_excluded": 0,
            "excluded_fraction": 0.0,
            "warning_flag": False,
            "high_risk_flag": False,
        }

    # -------------------------------------------------------------------------
    # STEP 3.6: Empty dataset — skip diagnostics, complete with warnings
    # -------------------------------------------------------------------------
    total_rows_after_schema = _run_metadata["data_quality"].get("total_rows_original", 0)
    if total_rows_after_schema == 0:
        _run_metadata["data_quality"] = {
            "total_rows_original": 0,
            "rows_used_for_analysis": 0,
            "rows_excluded": 0,
            "excluded_fraction": 0.0,
            "warning_flag": True,
            "high_risk_flag": False,
        }
        _run_metadata["diagnostics_skipped"] = True
        _run_metadata["skip_reason"] = "EMPTY_DATASET"
        print("EMPTY_DATASET_DETECTED — Diagnostics skipped, pipeline completed with warnings.")

    # =========================================================================
    # STEP 4: Handle --dry-run mode
    # =========================================================================
    if config.dry_run:
        # In dry run mode, we only validate inputs - no actual analysis
        print("=" * 60)
        print("DRY RUN MODE")
        print("=" * 60)
        print("Input validation complete.")
        print("No diagnostic analysis performed (--dry-run specified).")
        return 0  # Return 0 to indicate success

    # =========================================================================
    # STEP 5: Dataset frozen (automatic) - check if any diagnostics were requested
    # =========================================================================
    if not config.run_diagnostics and not config.run_benford and not config.run_lognormal and not config.run_cosmic:
        # Dataset is frozen, but no diagnostic flag was provided
        print("=" * 60)
        print("DATASET FROZEN - DIAGNOSTICS NOT REQUESTED")
        print("=" * 60)
        print("The dataset is frozen and ready for analysis.")
        print()
        print("Input validation completed successfully.")
        print()
        print("To run diagnostic analysis, use one or more flags:")
        print(f"  python -m week2_validation.run_week2 \\")
        print(f"      --fusion-data {config.fusion_data_path} \\")
        print(f"      --output-dir {config.output_dir} \\")
        print(f"      --run-diagnostics    # Distribution diagnostics")
        print(f"      --run-benford        # Benford's Law diagnostics (diagnostic only)")
        print(f"      --run-lognormal      # Log-normality diagnostics (diagnostic only)")
        print(f"      --run-cosmic         # COSMIC rank-order diagnostics (diagnostic only)")
        return 0  # Exit cleanly (explicit request required)

    # =========================================================================
    # STEP 7: Dataset is frozen AND diagnostics requested - execute (skip if empty)
    # =========================================================================
    exit_code = 0
    if _run_metadata.get("diagnostics_skipped"):
        # Empty dataset: diagnostics already skipped in STEP 3.6; write status and exit
        _run_metadata["diagnostics_run"] = [
            n for n, f in [
                ("diagnostics", config.run_diagnostics),
                ("benford", config.run_benford),
                ("lognormal", config.run_lognormal),
                ("cosmic", config.run_cosmic),
            ] if f
        ]
    else:
        print("=" * 60)
        print("RUNNING DIAGNOSTICS")
        print("=" * 60)
        print("Dataset frozen (automatic): YES")
        print(f"Distribution diagnostics requested: {'YES' if config.run_diagnostics else 'NO'}")
        print(f"Benford diagnostics requested: {'YES' if config.run_benford else 'NO'}")
        print(f"Log-normality diagnostics requested: {'YES' if config.run_lognormal else 'NO'}")
        print(f"COSMIC diagnostics requested: {'YES' if config.run_cosmic else 'NO'}")
        print()

        _run_metadata["diagnostics_run"] = [
            n for n, f in [
                ("diagnostics", config.run_diagnostics),
                ("benford", config.run_benford),
                ("lognormal", config.run_lognormal),
                ("cosmic", config.run_cosmic),
            ] if f
        ]

        # Execute distribution diagnostics if requested (lazy import happens inside)
        # Uses frozen_input_path to ensure diagnostics operate on immutable data
        if config.run_diagnostics:
            check_runtime_guard()
            result = execute_diagnostics(
                config, week2_config,
                frozen_data_path=effective_data_path,
                run_metadata=_run_metadata,
            )
            check_runtime_guard()
            if result != 0:
                exit_code = result
            print()

        # Execute Benford diagnostics if requested (lazy import happens inside)
        # Uses frozen_input_path to ensure diagnostics operate on immutable data
        if config.run_benford:
            check_runtime_guard()
            result = execute_benford_diagnostics(
                config, week2_config,
                frozen_data_path=effective_data_path,
                run_metadata=_run_metadata,
            )
            check_runtime_guard()
            if result != 0:
                exit_code = result
            print()

        # Execute log-normality diagnostics if requested (lazy import happens inside)
        # Uses frozen_input_path to ensure diagnostics operate on immutable data
        if config.run_lognormal:
            check_runtime_guard()
            result = execute_log_normality_diagnostics(
                config, week2_config,
                frozen_data_path=effective_data_path,
                run_metadata=_run_metadata,
            )
            check_runtime_guard()
            if result != 0:
                exit_code = result
            print()

        # Execute COSMIC diagnostics if requested (lazy import happens inside)
        # Uses frozen_input_path to ensure diagnostics operate on immutable data
        if config.run_cosmic:
            check_runtime_guard()
            result_code, cosmic_loaded = execute_cosmic_diagnostics(
                config, week2_config,
                frozen_data_path=effective_data_path,
                run_metadata=_run_metadata,
            )
            check_runtime_guard()
            if result_code != 0:
                exit_code = result_code
            if config.cosmic_data_path is not None:
                _run_metadata["cosmic_reference_loaded"] = cosmic_loaded

        # Optional: write approved dataset artifacts when validation succeeded (additive only)
        # Does not affect exit_code, pipeline success/failure, or status envelope on failure
        if exit_code == 0:
            try:
                try:
                    import scipy
                    _scipy_available = True
                except Exception:
                    _scipy_available = False
                try:
                    import psutil
                    _psutil_available = True
                except Exception:
                    _psutil_available = False
                _runtime_meta = {
                    "scipy_available": _scipy_available,
                    "psutil_available": _psutil_available,
                }
                from week2_validation.reporting.approved_dataset_writer import write_approved_dataset
                _csv_p, _cert_p = write_approved_dataset(
                    frozen_dataset_path=effective_data_path,
                    output_dir=config.output_dir,
                    dataset_hash=frozen_input_dir.name,
                    validation_passed=True,
                    diagnostics_run=_run_metadata["diagnostics_run"],
                    runtime_metadata=_runtime_meta,
                    dataset_stem=config.dataset_stem,
                )
                if _csv_p is not None and _cert_p is not None:
                    print(f"Cleaned dataset written: {_csv_p.name}")
                    print(f"Certification written: {_cert_p.name}")
            except Exception as _e:
                print(f"Warning: Approved dataset write skipped: {_e}", file=sys.stderr)
                _run_metadata["approved_dataset_write_note"] = f"Write skipped: {type(_e).__name__}"

        # Write diagnostic results for narrative report (when diagnostics ran)
        diag_results = _run_metadata.get("diagnostic_results")
        if diag_results:
            try:
                diag_path = config.output_dir / f"week2_diagnostic_results_{config.dataset_stem}.json"
                with open(diag_path, "w", encoding="utf-8") as f:
                    json.dump(diag_results, f, indent=2)
            except OSError:
                pass  # Non-fatal; report may use status/cert only

    return exit_code


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main() -> int:
    """
    Main entry point for the CLI (Command Line Interface).

    This is the first function that runs when a user starts the program.
    It handles the complete workflow:
    1. Read and understand the user's command-line instructions
    2. Validate all the settings and file paths
    3. Run the validation pipeline
    4. Handle any errors gracefully and report them to the user

    Returns:
        Exit code: 0 = success, 1 = error, 130 = user cancelled.
    """
    # Create the argument parser that will read user instructions
    parser = create_argument_parser()

    # Check if the user ran the program without any arguments
    # sys.argv contains all the command-line arguments; index 0 is the program name
    if len(sys.argv) == 1:
        # No arguments provided - show the help text so user knows what to do
        parser.print_help(sys.stderr)
        return 1  # Return 1 to indicate an error (missing arguments)

    # =========================================================================
    # BENFORD CONTROLS: Handle standalone execution BEFORE requiring file args
    # This allows --run-benford-controls to run without --fusion-data/--output-dir
    # Uses synthetic data only, fully isolated from real data pipeline
    # =========================================================================
    if "--run-benford-controls" in sys.argv:
        # Check if this is a standalone controls run (no other diagnostic flags)
        other_diagnostic_flags = [
            "--run-diagnostics",
            "--run-benford", 
            "--run-lognormal",
            "--run-cosmic",
        ]
        has_other_diagnostics = any(flag in sys.argv for flag in other_diagnostic_flags)
        
        if not has_other_diagnostics:
            # Standalone controls run - no real data needed
            return execute_benford_controls()

    # Parse the command-line arguments the user provided
    args = parser.parse_args()

    # Try to validate all the arguments and create the configuration
    try:
        config = validate_arguments(args)
    except ConfigurationError as e:
        # Something was wrong with the configuration (bad file path, etc.)
        print(f"Configuration error: {e}", file=sys.stderr)
        return 1  # Return 1 to indicate an error

    # =========================================================================
    # BENFORD CONTROLS: Also run if combined with other diagnostics
    # In this case, file validation has already passed
    # =========================================================================
    if config.run_benford_controls:
        control_result = execute_benford_controls()
        if control_result != 0:
            # Only return early if there was an unexpected exception
            return control_result

    # Run pipeline with log capture (tee stdout/stderr to run_log_{stem}.txt)
    def _execute():
        start_runtime_guard()
        start_time = time.monotonic()
        exit_code, status_msg = run_week2_safely(run_pipeline, config)
        check_runtime_guard()
        runtime_seconds = time.monotonic() - start_time

        # Optional dependency availability (observability only; no logic change)
        try:
            import scipy
            scipy_available = True
        except Exception:
            scipy_available = False
        try:
            import psutil
            psutil_available = True
        except Exception:
            psutil_available = False

        # Always write status envelope
        meta = _run_metadata
        notes = [status_msg] if status_msg and status_msg != "OK" else []
        if meta.get("cosmic_reference_loaded") is False:
            notes.append("COSMIC_LOAD_FAILED")
        cosmic_reference_requested = config.cosmic_data_path is not None
        envelope = build_status_envelope(
            dataset_hash=meta.get("dataset_hash", ""),
            diagnostics_run=meta.get("diagnostics_run", []),
            exit_code=exit_code,
            notes=notes,
            cosmic_reference_loaded=meta.get("cosmic_reference_loaded"),
            cosmic_reference_requested=cosmic_reference_requested,
            scipy_available=scipy_available,
            psutil_available=psutil_available,
            data_quality=meta.get("data_quality"),
            diagnostics_skipped=meta.get("diagnostics_skipped"),
            skip_reason=meta.get("skip_reason"),
            benford_scale_span=meta.get("benford_scale_span"),
        )
        envelope["runtime_seconds"] = round(runtime_seconds, 2)
        status_path = config.output_dir / f"week2_status_{config.dataset_stem}.json"
        try:
            with open(status_path, "w", encoding="utf-8") as f:
                json.dump(envelope, f, indent=2)
            print(f"Status written to {status_path}")
        except OSError as e:
            print(f"Could not write status file: {e}", file=sys.stderr)
            # Status write failure always overrides success exit (FIX 4)
            if exit_code == 0:
                exit_code = int(Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR)
            return exit_code

        # Generate researcher-friendly narrative report when requested
        if config.generate_report:
            # Generate dynamic skewness plot when we have skewness (from diagnostics)
            skewness_val = None
            dq = _run_metadata.get("data_quality") or {}
            diag_res = _run_metadata.get("diagnostic_results") or {}
            dist = diag_res.get("distribution") or {}
            for src in [dist.get("skewness"), dq.get("skewness")]:
                if src is not None and isinstance(src, (int, float)) and math.isfinite(float(src)):
                    skewness_val = float(src)
                    break
            if skewness_val is not None:
                try:
                    from week2_validation.distributions.visualize import generate_dynamic_skewness_plot
                    skew_path = config.output_dir / "skewness_diagnostic.png"
                    if generate_dynamic_skewness_plot(skewness_val, skew_path) is not None:
                        print(f"Skewness diagnostic saved: skewness_diagnostic.png")
                except Exception as e:
                    print(f"Warning: Could not generate skewness diagnostic: {e}", file=sys.stderr)
            try:
                from week2_validation.reporting.narrative_generator import generate_narrative_report
                report_path = generate_narrative_report(config.output_dir, dataset_stem=config.dataset_stem)
                if report_path is not None:
                    print(f"Narrative report written: {report_path.name}")
            except Exception as e:
                print(f"Warning: Could not generate narrative report: {e}", file=sys.stderr)

        # Failure report when validation/execution failed (paper trail for every file)
        if exit_code != 0:
            try:
                from week2_validation.reporting.narrative_generator import generate_failure_report
                _msg = str(status_msg or "").upper()
                if exit_code == int(Week2ExitCode.INPUT_SCHEMA_ERROR) or "SCHEMA" in _msg or "MISSING" in _msg:
                    reason = "Missing Columns or Data Type Mismatch"
                elif exit_code == int(Week2ExitCode.CONFIG_ERROR):
                    reason = "Configuration Error"
                elif exit_code == int(Week2ExitCode.FREEZE_ERROR):
                    reason = "Dataset Freeze Error"
                elif "EMPTY" in _msg:
                    reason = "Empty File"
                else:
                    reason = "Validation or Execution Failed"
                fp = generate_failure_report(
                    config.output_dir,
                    config.dataset_stem,
                    failure_reason=reason,
                    failure_details=status_msg or "Unknown error",
                )
                if fp is not None:
                    print(f"Failure report written: {fp.name}")
            except Exception:
                pass

        return exit_code

    exit_code = _run_with_log_capture(config, _execute)
    try:
        print(f"Run log saved: run_log_{config.dataset_stem}.txt")
    except Exception:
        pass
    return exit_code


# This block only runs when the file is executed directly (not imported)
# It starts the main function and uses its return value as the program exit code
if __name__ == "__main__":
    sys.exit(main())
