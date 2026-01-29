"""
Data loader module for Week 2 validation pipeline.

This module provides secure, format-agnostic data loading capabilities
for cancer genomics fusion datasets. It supports multiple file formats
and enforces strict schema validation to ensure data integrity.

Supported formats:
    - CSV (.csv)
    - TSV (.tsv)
    - JSON (.json)
    - Parquet (.parquet)
    - Excel (.xlsx, .xls)

Security considerations:
    - No execution of arbitrary code
    - No unsafe deserialization
    - Path validation against directory traversal
    - No logging of raw data contents
"""

from pathlib import Path
from typing import Optional, Set

import pandas as pd


class DataLoaderError(Exception):
    """Base exception for data loader errors."""

    pass


class UnsupportedFormatError(DataLoaderError):
    """Raised when file format is not supported."""

    pass


class FileValidationError(DataLoaderError):
    """Raised when file path validation fails."""

    pass


class SchemaValidationError(DataLoaderError):
    """Raised when required schema fields are missing."""

    pass


SUPPORTED_EXTENSIONS: Set[str] = {".csv", ".tsv", ".json", ".parquet", ".xlsx", ".xls"}

REQUIRED_FUSION_FIELDS: Set[str] = {"fusion_id", "protein_length", "recurrence_count"}


def validate_file_path(file_path: str, must_exist: bool = True) -> Path:
    """
    Validate and resolve a file path securely.

    Ensures the path is safe from directory traversal attacks and
    optionally verifies the file exists.

    Args:
        file_path: The file path string to validate.
        must_exist: If True, verify the file exists on disk.

    Returns:
        A resolved Path object representing the validated path.

    Raises:
        FileValidationError: If the path is invalid, contains traversal
            attempts, or the file does not exist when must_exist is True.
    """
    if not file_path:
        raise FileValidationError("File path cannot be empty")

    if not isinstance(file_path, str):
        raise FileValidationError(
            f"File path must be a string, got {type(file_path).__name__}"
        )

    try:
        path = Path(file_path)
    except (TypeError, ValueError) as e:
        raise FileValidationError(f"Invalid file path format: {e}") from e

    try:
        resolved_path = path.resolve(strict=False)
    except (OSError, RuntimeError) as e:
        raise FileValidationError(f"Cannot resolve file path: {e}") from e

    path_str = str(path)
    if ".." in path_str.split("/") or ".." in path_str.split("\\"):
        raise FileValidationError(
            "Directory traversal patterns detected in path"
        )

    if must_exist:
        if not resolved_path.exists():
            raise FileValidationError(f"File does not exist: {resolved_path}")
        if not resolved_path.is_file():
            raise FileValidationError(f"Path is not a file: {resolved_path}")

    return resolved_path


def validate_output_directory(dir_path: str, create: bool = False) -> Path:
    """
    Validate an output directory path.

    Args:
        dir_path: The directory path string to validate.
        create: If True, create the directory if it does not exist.

    Returns:
        A resolved Path object representing the validated directory.

    Raises:
        FileValidationError: If the path is invalid or cannot be used
            as an output directory.
    """
    if not dir_path:
        raise FileValidationError("Output directory path cannot be empty")

    if not isinstance(dir_path, str):
        raise FileValidationError(
            f"Directory path must be a string, got {type(dir_path).__name__}"
        )

    try:
        path = Path(dir_path)
    except (TypeError, ValueError) as e:
        raise FileValidationError(f"Invalid directory path format: {e}") from e

    try:
        resolved_path = path.resolve(strict=False)
    except (OSError, RuntimeError) as e:
        raise FileValidationError(f"Cannot resolve directory path: {e}") from e

    path_str = str(path)
    if ".." in path_str.split("/") or ".." in path_str.split("\\"):
        raise FileValidationError(
            "Directory traversal patterns detected in path"
        )

    if resolved_path.exists():
        if not resolved_path.is_dir():
            raise FileValidationError(
                f"Path exists but is not a directory: {resolved_path}"
            )
    elif create:
        try:
            resolved_path.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            raise FileValidationError(
                f"Cannot create output directory: {e}"
            ) from e
    else:
        raise FileValidationError(f"Directory does not exist: {resolved_path}")

    return resolved_path


def get_file_format(file_path: Path) -> str:
    """
    Determine file format from extension.

    Args:
        file_path: Path object to examine.

    Returns:
        Lowercase file extension including the leading dot.

    Raises:
        UnsupportedFormatError: If the file extension is not supported.
    """
    extension = file_path.suffix.lower()

    if not extension:
        raise UnsupportedFormatError(
            f"File has no extension, cannot determine format: {file_path.name}"
        )

    if extension not in SUPPORTED_EXTENSIONS:
        raise UnsupportedFormatError(
            f"Unsupported file format '{extension}'. "
            f"Supported formats: {', '.join(sorted(SUPPORTED_EXTENSIONS))}"
        )

    return extension


def validate_schema(df: pd.DataFrame, required_fields: Set[str]) -> None:
    """
    Validate that a DataFrame contains all required fields.

    Args:
        df: The pandas DataFrame to validate.
        required_fields: Set of column names that must be present.

    Raises:
        SchemaValidationError: If any required fields are missing.
    """
    if df is None:
        raise SchemaValidationError("DataFrame is None")

    if not isinstance(df, pd.DataFrame):
        raise SchemaValidationError(
            f"Expected pandas DataFrame, got {type(df).__name__}"
        )

    actual_columns = set(df.columns)
    missing_fields = required_fields - actual_columns

    if missing_fields:
        raise SchemaValidationError(
            f"Missing required fields: {sorted(missing_fields)}. "
            f"Available fields: {sorted(actual_columns)}"
        )


def load_csv(file_path: Path) -> pd.DataFrame:
    """
    Load a CSV file into a DataFrame.

    Args:
        file_path: Path to the CSV file.

    Returns:
        pandas DataFrame containing the loaded data.

    Raises:
        DataLoaderError: If the file cannot be read.
    """
    try:
        return pd.read_csv(file_path, low_memory=False)
    except pd.errors.EmptyDataError as e:
        raise DataLoaderError(f"CSV file is empty: {file_path.name}") from e
    except pd.errors.ParserError as e:
        raise DataLoaderError(f"CSV parsing error: {e}") from e


def load_tsv(file_path: Path) -> pd.DataFrame:
    """
    Load a TSV file into a DataFrame.

    Args:
        file_path: Path to the TSV file.

    Returns:
        pandas DataFrame containing the loaded data.

    Raises:
        DataLoaderError: If the file cannot be read.
    """
    try:
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    except pd.errors.EmptyDataError as e:
        raise DataLoaderError(f"TSV file is empty: {file_path.name}") from e
    except pd.errors.ParserError as e:
        raise DataLoaderError(f"TSV parsing error: {e}") from e


def load_json(file_path: Path) -> pd.DataFrame:
    """
    Load a JSON file into a DataFrame.

    Only supports JSON structures that pandas can safely interpret
    as tabular data. Does not execute arbitrary code.

    Args:
        file_path: Path to the JSON file.

    Returns:
        pandas DataFrame containing the loaded data.

    Raises:
        DataLoaderError: If the file cannot be read or parsed.
    """
    try:
        return pd.read_json(file_path)
    except ValueError as e:
        raise DataLoaderError(f"JSON parsing error: {e}") from e


def load_parquet(file_path: Path) -> pd.DataFrame:
    """
    Load a Parquet file into a DataFrame.

    Args:
        file_path: Path to the Parquet file.

    Returns:
        pandas DataFrame containing the loaded data.

    Raises:
        DataLoaderError: If the file cannot be read.
    """
    try:
        return pd.read_parquet(file_path)
    except Exception as e:
        raise DataLoaderError(f"Parquet reading error: {e}") from e


def load_excel(file_path: Path) -> pd.DataFrame:
    """
    Load an Excel file into a DataFrame.

    Reads only the first sheet. Does not execute macros or
    embedded code.

    Args:
        file_path: Path to the Excel file.

    Returns:
        pandas DataFrame containing the loaded data.

    Raises:
        DataLoaderError: If the file cannot be read.
    """
    try:
        return pd.read_excel(file_path, engine="openpyxl")
    except ValueError as e:
        raise DataLoaderError(f"Excel reading error: {e}") from e


def load_data(file_path: str) -> pd.DataFrame:
    """
    Load data from a supported file format into a DataFrame.

    This is the primary entry point for loading data files.
    Format is inferred strictly from file extension.

    Args:
        file_path: Path to the data file.

    Returns:
        pandas DataFrame containing the loaded data.

    Raises:
        FileValidationError: If the file path is invalid or file does not exist.
        UnsupportedFormatError: If the file format is not supported.
        DataLoaderError: If the file cannot be read or parsed.
    """
    validated_path = validate_file_path(file_path, must_exist=True)
    file_format = get_file_format(validated_path)

    loaders = {
        ".csv": load_csv,
        ".tsv": load_tsv,
        ".json": load_json,
        ".parquet": load_parquet,
        ".xlsx": load_excel,
        ".xls": load_excel,
    }

    loader_func = loaders.get(file_format)
    if loader_func is None:
        raise UnsupportedFormatError(f"No loader implemented for {file_format}")

    return loader_func(validated_path)


def load_fusion_data(file_path: str) -> pd.DataFrame:
    """
    Load and validate fusion dataset.

    Loads data from the specified file and validates that all
    required fusion schema fields are present.

    Required fields:
        - fusion_id
        - protein_length
        - recurrence_count

    Args:
        file_path: Path to the fusion data file.

    Returns:
        pandas DataFrame containing validated fusion data.

    Raises:
        FileValidationError: If the file path is invalid.
        UnsupportedFormatError: If the file format is not supported.
        DataLoaderError: If the file cannot be read.
        SchemaValidationError: If required fields are missing.
    """
    df = load_data(file_path)
    validate_schema(df, REQUIRED_FUSION_FIELDS)
    return df


def load_reference_data(file_path: Optional[str]) -> Optional[pd.DataFrame]:
    """
    Load optional reference data (e.g., COSMIC).

    If no path is provided, returns None. Reference data does not
    undergo fusion schema validation.

    Args:
        file_path: Path to the reference data file, or None.

    Returns:
        pandas DataFrame containing reference data, or None if no path provided.

    Raises:
        FileValidationError: If the file path is invalid.
        UnsupportedFormatError: If the file format is not supported.
        DataLoaderError: If the file cannot be read.
    """
    if file_path is None:
        return None

    return load_data(file_path)