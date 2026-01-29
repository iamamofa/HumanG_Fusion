#!/usr/bin/env python3
"""
Week 2 Validation Pipeline - CLI Entry Point.

This module provides the command-line interface for the Week 2 data
integrity and statistical validation pipeline. It handles argument
parsing, input validation, and orchestration of the validation workflow.

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

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from week2_validation.utils.data_loader import (
    DataLoaderError,
    FileValidationError,
    SchemaValidationError,
    UnsupportedFormatError,
    load_fusion_data,
    load_reference_data,
    validate_file_path,
    validate_output_directory,
)


@dataclass(frozen=True)
class PipelineConfig:
    """
    Immutable configuration for the validation pipeline.

    Attributes:
        fusion_data_path: Validated path to the fusion dataset.
        output_dir: Validated path to the output directory.
        cosmic_data_path: Optional validated path to COSMIC reference data.
        dry_run: If True, validate inputs only without running diagnostics.
    """

    fusion_data_path: Path
    output_dir: Path
    cosmic_data_path: Optional[Path]
    dry_run: bool


class PipelineError(Exception):
    """Base exception for pipeline execution errors."""

    pass


class ConfigurationError(PipelineError):
    """Raised when pipeline configuration is invalid."""

    pass


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser.

    Returns:
        Configured ArgumentParser instance.
    """
    parser = argparse.ArgumentParser(
        prog="week2_validation",
        description=(
            "Week 2 Data Integrity & Statistical Validation Pipeline. "
            "Performs diagnostic validation on fusion datasets. "
            "This pipeline is diagnostic only: no hypothesis testing, "
            "no modeling, no biological claims."
        ),
        epilog=(
            "Example: python -m week2_validation.run_week2 "
            "--fusion-data fusion.csv --output-dir ./results"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    required_group = parser.add_argument_group("required arguments")

    required_group.add_argument(
        "--fusion-data",
        type=str,
        required=True,
        metavar="PATH",
        help=(
            "Path to the fusion dataset file. "
            "Supported formats: CSV, TSV, JSON, Parquet, Excel."
        ),
    )

    required_group.add_argument(
        "--output-dir",
        type=str,
        required=True,
        metavar="PATH",
        help="Directory where validation results will be written.",
    )

    optional_group = parser.add_argument_group("optional arguments")

    optional_group.add_argument(
        "--cosmic-data",
        type=str,
        required=False,
        default=None,
        metavar="PATH",
        help=(
            "Path to COSMIC reference data file. "
            "If provided, enables reference-based validation."
        ),
    )

    optional_group.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help=(
            "Validate inputs and configuration only. "
            "Do not execute diagnostic analysis."
        ),
    )

    return parser


def validate_arguments(args: argparse.Namespace) -> PipelineConfig:
    """
    Validate all command-line arguments and construct pipeline configuration.

    Performs comprehensive validation of all provided paths and options,
    ensuring they meet security and format requirements.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Validated PipelineConfig instance.

    Raises:
        ConfigurationError: If any argument fails validation.
    """
    try:
        fusion_path = validate_file_path(args.fusion_data, must_exist=True)
    except FileValidationError as e:
        raise ConfigurationError(f"Invalid fusion data path: {e}") from e

    try:
        output_dir = validate_output_directory(args.output_dir, create=True)
    except FileValidationError as e:
        raise ConfigurationError(f"Invalid output directory: {e}") from e

    cosmic_path: Optional[Path] = None
    if args.cosmic_data is not None:
        try:
            cosmic_path = validate_file_path(args.cosmic_data, must_exist=True)
        except FileValidationError as e:
            raise ConfigurationError(f"Invalid COSMIC data path: {e}") from e

    return PipelineConfig(
        fusion_data_path=fusion_path,
        output_dir=output_dir,
        cosmic_data_path=cosmic_path,
        dry_run=args.dry_run,
    )


def validate_inputs(config: PipelineConfig) -> None:
    """
    Validate that input files can be loaded and meet schema requirements.

    This function attempts to load all input files and validates their
    structure without performing any analysis.

    Args:
        config: Validated pipeline configuration.

    Raises:
        PipelineError: If any input file cannot be loaded or fails validation.
    """
    print(f"Validating fusion data: {config.fusion_data_path.name}")

    try:
        fusion_df = load_fusion_data(str(config.fusion_data_path))
    except (FileValidationError, UnsupportedFormatError) as e:
        raise PipelineError(f"Cannot load fusion data: {e}") from e
    except DataLoaderError as e:
        raise PipelineError(f"Error reading fusion data: {e}") from e
    except SchemaValidationError as e:
        raise PipelineError(f"Fusion data schema validation failed: {e}") from e

    row_count = len(fusion_df)
    column_count = len(fusion_df.columns)
    print(f"  Loaded {row_count} records with {column_count} fields")
    print("  Schema validation: PASSED")

    if config.cosmic_data_path is not None:
        print(f"Validating COSMIC reference data: {config.cosmic_data_path.name}")

        try:
            cosmic_df = load_reference_data(str(config.cosmic_data_path))
        except (FileValidationError, UnsupportedFormatError) as e:
            raise PipelineError(f"Cannot load COSMIC data: {e}") from e
        except DataLoaderError as e:
            raise PipelineError(f"Error reading COSMIC data: {e}") from e

        if cosmic_df is not None:
            row_count = len(cosmic_df)
            column_count = len(cosmic_df.columns)
            print(f"  Loaded {row_count} records with {column_count} fields")

    print(f"Output directory validated: {config.output_dir}")


def run_pipeline(config: PipelineConfig) -> int:
    """
    Execute the validation pipeline.

    This function orchestrates the validation workflow based on the
    provided configuration.

    Args:
        config: Validated pipeline configuration.

    Returns:
        Exit code (0 for success, non-zero for failure).

    Raises:
        PipelineError: If pipeline execution fails.
    """
    print("=" * 60)
    print("Week 2 Validation Pipeline")
    print("=" * 60)
    print()

    validate_inputs(config)
    print()

    if config.dry_run:
        print("Dry run mode: Input validation complete.")
        print("No diagnostic analysis performed.")
        return 0

    print("Input validation complete.")
    print("Pipeline is ready for diagnostic execution.")
    print()
    print("Note: Diagnostic analysis modules not yet implemented.")
    print("This entry point validates orchestration readiness only.")

    return 0


def main() -> int:
    """
    Main entry point for the CLI.

    Parses arguments, validates configuration, and executes the pipeline.
    All errors are caught and reported with appropriate exit codes.

    Returns:
        Exit code (0 for success, non-zero for failure).
    """
    parser = create_argument_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    try:
        config = validate_arguments(args)
    except ConfigurationError as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        return 1

    try:
        return run_pipeline(config)
    except PipelineError as e:
        print(f"Pipeline error: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())