#!/usr/bin/env python3
"""
Week 2 Validation Pipeline - CLI Entry Point.

This module provides the command-line interface for the Week 2 data
integrity and statistical validation pipeline. It handles argument
parsing, input validation, and orchestration of the validation workflow.

WHAT DOES THIS FILE DO?
This is the main "control center" for running the Week 2 validation checks.
When a user runs this program from the command line, this file:
1. Reads the user's instructions (which files to analyze, where to save results)
2. Loads configuration from thresholds.yaml
3. Checks dataset freeze state (Option C: flag file takes precedence over CLI)
4. Validates all input files exist and are in the correct format
5. If --run-diagnostics is requested AND dataset is frozen, runs diagnostics
6. Reports the results

EXECUTION FLOW:
    ┌─────────────────────────────────────────────────────────────────┐
    │  1. Parse CLI arguments                                         │
    │  2. Load thresholds.yaml configuration                          │
    │  3. Check freeze state (flag file > CLI > error)                │
    │  4. Validate input files                                        │
    │  5. If --dry-run: exit after validation                         │
    │  6. If NOT frozen: exit with message                            │
    │  7. If frozen but NO --run-diagnostics: exit with message       │
    │  8. If frozen AND --run-diagnostics: run diagnostics            │
    └─────────────────────────────────────────────────────────────────┘

This pipeline is diagnostic only:
    - No hypothesis testing
    - No modeling
    - No biological claims

Usage:
    python -m week2_validation.run_week2 \\
        --fusion-data /path/to/fusion_data.csv \\
        --output-dir /path/to/output \\
        [--cosmic-data /path/to/cosmic.tsv] \\
        [--run-diagnostics] \\
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
# 'sys' provides access to system functions like exiting the program
import sys
# 'dataclass' creates simple classes for holding related data together
from dataclasses import dataclass
# 'Path' helps work with file and folder locations on the computer
from pathlib import Path
# 'Optional' indicates that a value might be present or might be None (empty)
from typing import Optional

# =============================================================================
# INTERNAL IMPORTS - Data Loading
# =============================================================================

# Import our custom data loading tools from the utils folder
# These handle reading files and checking they're in the correct format
from week2_validation.utils.data_loader import (
    DataLoaderError,           # Error when data can't be read
    FileValidationError,       # Error when file path is invalid
    SchemaValidationError,     # Error when data is missing required columns
    UnsupportedFormatError,    # Error when file type isn't supported
    load_fusion_data,          # Function to load the main fusion dataset
    load_reference_data,       # Function to load optional reference data
    validate_file_path,        # Function to check if a file path is valid
    validate_output_directory, # Function to check if output folder is valid
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
        run_diagnostics: If True, run diagnostic analysis (requires frozen dataset).
    """

    fusion_data_path: Path          # Where the fusion data file is located
    output_dir: Path                # Where to save the output/results
    cosmic_data_path: Optional[Path]  # Optional reference data location (can be empty)
    dry_run: bool                   # True = just validate, False = run full analysis
    run_diagnostics: bool           # True = run diagnostics (if frozen), False = skip


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
            "Week 2 Data Integrity & Statistical Validation Pipeline. "
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

    # Define the run-diagnostics flag (OPTIONAL)
    # This explicitly requests diagnostic analysis to run
    # IMPORTANT: Diagnostics only run if dataset is also frozen
    optional_group.add_argument(
        "--run-diagnostics",
        action="store_true",       # Just a flag, no value needed
        default=False,             # Off by default (explicit request required)
        help=(
            "Run diagnostic analysis on the dataset. "
            "REQUIRES dataset to be frozen. "
            "If not specified, only input validation is performed."
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
    # This packages all validated settings together
    return PipelineConfig(
        fusion_data_path=fusion_path,
        output_dir=output_dir,
        cosmic_data_path=cosmic_path,
        dry_run=args.dry_run,
        run_diagnostics=args.run_diagnostics,
    )


# =============================================================================
# INPUT VALIDATION
# =============================================================================

def validate_inputs(config: PipelineConfig) -> None:
    """
    Validate that input files can be loaded and meet schema requirements.

    This function does a "test load" of all input files to make sure:
    1. The files can actually be read (not corrupted)
    2. The data has the required columns (correct structure)
    
    This catches problems early before running the full analysis.

    Args:
        config: The validated pipeline configuration.

    Raises:
        PipelineError: If any input file cannot be loaded or is malformed.
    """
    # Tell the user what we're checking
    print(f"Validating fusion data: {config.fusion_data_path.name}")

    # Try to load the fusion data file
    try:
        # Load the data into memory as a table (DataFrame)
        fusion_df = load_fusion_data(str(config.fusion_data_path))
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

def execute_diagnostics(config: PipelineConfig, week2_config: Week2Config) -> int:
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
        from week2_validation.distributions.visualize import run_diagnostics
    except ImportError as e:
        print(f"Error: Cannot import diagnostic module: {e}", file=sys.stderr)
        print("Make sure all required dependencies are installed.", file=sys.stderr)
        return 1
    
    # -------------------------------------------------------------------------
    # Load the fusion data for diagnostic analysis
    # -------------------------------------------------------------------------
    try:
        fusion_df = load_fusion_data(str(config.fusion_data_path))
    except Exception as e:
        print(f"Error loading fusion data for diagnostics: {e}", file=sys.stderr)
        return 1
    
    # -------------------------------------------------------------------------
    # Extract protein_length column for analysis
    # No assumptions about data format - just use what's in the required column
    # -------------------------------------------------------------------------
    if "protein_length" not in fusion_df.columns:
        print("Error: 'protein_length' column not found in fusion data.", file=sys.stderr)
        return 1
    
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
        print()
        print(f"  Log mode used: {result.log_mode_used}")
        print(f"  Diagnostic only: {result.metadata.diagnostic_only}")
        print(f"  Hypothesis tested: {result.metadata.hypothesis_tested}")
        print()
        print("Diagnostic analysis complete.")
        print("NOTE: Results are diagnostic only. No files saved. No interpretation provided.")
        
    except Exception as e:
        print(f"Error during diagnostic analysis: {e}", file=sys.stderr)
        return 1
    
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
    # Print a header banner to clearly show the pipeline is starting
    print("=" * 60)
    print("Week 2 Validation Pipeline")
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
    # STEP 2: Check dataset freeze state
    # =========================================================================
    print("Checking dataset freeze state...")
    
    # Resolve flag file path relative to current working directory
    flag_file_path = Path.cwd() / DEFAULT_FLAG_FILE_PATH
    if not flag_file_path.parent.exists():
        # Try relative to this file's location
        flag_file_path = Path(__file__).parent / ".frozen"
    
    freeze_state = check_freeze_state(flag_file_path)
    
    print(f"  Flag file path: {flag_file_path}")
    print(f"  Freeze state: {'FROZEN' if freeze_state.is_frozen else 'NOT FROZEN'}")
    print(f"  State source: {freeze_state.source}")
    print()

    # =========================================================================
    # STEP 3: Validate all input files
    # =========================================================================
    print("Validating inputs...")
    validate_inputs(config)
    print()

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
    # STEP 5: Check if dataset is frozen (required for any analysis)
    # =========================================================================
    if not freeze_state.is_frozen:
        # Dataset is NOT frozen - cannot proceed with any analysis
        print("=" * 60)
        print("DATASET NOT FROZEN")
        print("=" * 60)
        print("The dataset is not yet frozen.")
        print()
        print("To run diagnostics, the dataset must first be frozen by:")
        print(f"  1. Creating the flag file: {flag_file_path}")
        print("  2. Or completing Week 1 data freeze process")
        print()
        print("Input validation completed successfully.")
        print("No diagnostic analysis performed (dataset not frozen).")
        return 0  # Exit cleanly (not an error, just a gate)

    # =========================================================================
    # STEP 6: Dataset is frozen - check if diagnostics were requested
    # =========================================================================
    if not config.run_diagnostics:
        # Dataset is frozen, but --run-diagnostics flag was not provided
        print("=" * 60)
        print("DATASET FROZEN - DIAGNOSTICS NOT REQUESTED")
        print("=" * 60)
        print("The dataset is frozen and ready for analysis.")
        print()
        print("Input validation completed successfully.")
        print()
        print("To run diagnostic analysis, use the --run-diagnostics flag:")
        print(f"  python -m week2_validation.run_week2 \\")
        print(f"      --fusion-data {config.fusion_data_path} \\")
        print(f"      --output-dir {config.output_dir} \\")
        print(f"      --run-diagnostics")
        return 0  # Exit cleanly (explicit request required)

    # =========================================================================
    # STEP 7: Dataset is frozen AND --run-diagnostics requested
    # =========================================================================
    print("=" * 60)
    print("RUNNING DIAGNOSTICS")
    print("=" * 60)
    print("Dataset is frozen: YES")
    print("Diagnostics requested: YES")
    print()
    
    # Execute diagnostics (lazy import happens inside)
    return execute_diagnostics(config, week2_config)


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

    # Parse the command-line arguments the user provided
    args = parser.parse_args()

    # Try to validate all the arguments and create the configuration
    try:
        config = validate_arguments(args)
    except ConfigurationError as e:
        # Something was wrong with the configuration (bad file path, etc.)
        print(f"Configuration error: {e}", file=sys.stderr)
        return 1  # Return 1 to indicate an error

    # Try to run the actual pipeline
    try:
        return run_pipeline(config)
    except PipelineError as e:
        # Something went wrong during pipeline execution
        print(f"Pipeline error: {e}", file=sys.stderr)
        return 1  # Return 1 to indicate an error
    except FreezeStateError as e:
        # Freeze state issue
        print(f"Freeze state error: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        # User pressed Ctrl+C to stop the program
        print("\nPipeline interrupted by user.", file=sys.stderr)
        return 130  # Standard exit code for keyboard interrupt


# This block only runs when the file is executed directly (not imported)
# It starts the main function and uses its return value as the program exit code
if __name__ == "__main__":
    sys.exit(main())
