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
2. Checks that all the specified files exist and are in the correct format
3. Runs the validation analysis
4. Reports the results

This pipeline is diagnostic only:
    - No hypothesis testing
    - No modeling
    - No biological claims

Usage:
    python -m week2_validation.run_week2 \\
        --fusion-data /path/to/fusion_data.csv \\
        --output-dir /path/to/output \\
        [--cosmic-data /path/to/cosmic.tsv] \\
        [--dry-run]

Security considerations:
    - All file paths are validated before use
    - No directory traversal is permitted
    - Output is restricted to the specified output directory
    - No sensitive data is logged or printed
"""

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
    
    Think of it as a filled-out form that tells the pipeline what to do.

    Attributes:
        fusion_data_path: Location of the main fusion dataset file.
        output_dir: Folder where results will be saved.
        cosmic_data_path: Optional location of COSMIC reference data.
        dry_run: If True, only check inputs without running full analysis.
    """

    fusion_data_path: Path          # Where the fusion data file is located
    output_dir: Path                # Where to save the output/results
    cosmic_data_path: Optional[Path]  # Optional reference data location (can be empty)
    dry_run: bool                   # True = just validate, False = run full analysis


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
            "--fusion-data fusion.csv --output-dir ./results"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

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
    )


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


def run_pipeline(config: PipelineConfig) -> int:
    """
    Execute the validation pipeline.

    This is the main function that runs the entire validation process.
    It coordinates all the steps: validation, analysis, and reporting.

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

    # STEP 1: Validate all input files before doing any analysis
    # This ensures we catch problems early
    validate_inputs(config)
    print()

    # STEP 2: Check if this is a dry run (test mode)
    if config.dry_run:
        # In dry run mode, we only validate inputs - no actual analysis
        print("Dry run mode: Input validation complete.")
        print("No diagnostic analysis performed.")
        return 0  # Return 0 to indicate success

    # STEP 3: If not a dry run, proceed with the full pipeline
    print("Input validation complete.")
    print("Pipeline is ready for diagnostic execution.")
    print()
    # Note: The actual analysis modules will be added in future development
    print("Note: Diagnostic analysis modules not yet implemented.")
    print("This entry point validates orchestration readiness only.")

    return 0  # Return 0 to indicate success


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
    except KeyboardInterrupt:
        # User pressed Ctrl+C to stop the program
        print("\nPipeline interrupted by user.", file=sys.stderr)
        return 130  # Standard exit code for keyboard interrupt


# This block only runs when the file is executed directly (not imported)
# It starts the main function and uses its return value as the program exit code
if __name__ == "__main__":
    sys.exit(main())