"""
Data loader module for Week 2 validation pipeline.

This module provides secure, format-agnostic data loading capabilities
for cancer genomics fusion datasets. It supports multiple file formats
and enforces strict schema validation to ensure data integrity.

WHAT DOES THIS FILE DO?
This is the "file reader" for the validation pipeline. It knows how to:
1. Read data from different file types (CSV, Excel, JSON, etc.)
2. Check that files are safe and valid before reading
3. Verify that the data has the required columns (schema validation)

Think of it as a careful librarian who checks books before lending them out.

Supported formats:
    - CSV (.csv) - Comma-separated values, like a simple spreadsheet
    - TSV (.tsv) - Tab-separated values, similar to CSV but uses tabs
    - JSON (.json) - JavaScript Object Notation, a structured text format
    - Parquet (.parquet) - A compressed format for large datasets
    - Excel (.xlsx) - Microsoft Excel spreadsheet files (.xlsx only)

Security considerations:
    - No execution of arbitrary code
    - No unsafe deserialization
    - Path validation against directory traversal
    - No logging of raw data contents
"""

# 'Path' helps work with file and folder locations on the computer
import logging
from pathlib import Path
# 'Optional' means a value can be present or None; 'Set' is a collection of unique items
from typing import Optional, Set

# 'pandas' is a powerful library for working with tabular data (like spreadsheets)
# We abbreviate it as 'pd' for convenience
import pandas as pd

_logger = logging.getLogger(__name__)


# ----- CUSTOM ERROR TYPES -----
# These help us identify exactly what went wrong when loading data fails

class DataLoaderError(Exception):
    """
    Base exception for data loader errors.
    
    This is a general error that occurs when data cannot be loaded.
    All other data loader errors inherit from this one.
    """

    pass


class UnsupportedFormatError(DataLoaderError):
    """
    Raised when file format is not supported.
    
    This error occurs when someone tries to load a file type we don't
    know how to read (like a .pdf or .docx file).
    """

    pass


class FileValidationError(DataLoaderError):
    """
    Raised when file path validation fails.
    
    This error occurs when:
    - The file doesn't exist
    - The path is empty
    - The path contains suspicious patterns (like ".." for directory traversal)
    """

    pass


class SchemaValidationError(DataLoaderError):
    """
    Raised when required schema fields are missing.
    
    This error occurs when the data file is missing one or more
    required columns. For example, fusion data must have 'fusion_id'.
    """

    pass


# ----- CONFIGURATION CONSTANTS -----
# These define what file types we support and what columns we require

# All the file extensions (types) that we know how to read
# Excel: .xlsx only (legacy .xls not supported)
SUPPORTED_EXTENSIONS: Set[str] = {".csv", ".tsv", ".json", ".parquet", ".xlsx"}

# The columns that MUST be present in fusion data files (matches schema fusion_schema.yaml)
# Without these columns, we cannot perform the validation analysis
REQUIRED_FUSION_FIELDS: Set[str] = {
    "fusion_id",
    "gene_1",
    "gene_2",
    "protein_length",
    "recurrence_count",
}

# Memory exhaustion protection: max file size before load (5GB; aligned with runtime guard)
MAX_FILE_SIZE_BYTES = 5 * 1024 * 1024 * 1024
# CSV/TSV row estimate: avg bytes per row for guard; reject if estimated rows exceed
AVG_CSV_ROW_BYTES = 200
MAX_ESTIMATED_CSV_ROWS = 200_000_000  # 200M estimated rows ceiling

# Soft warning thresholds (observability only; non-blocking)
SOFT_WARNING_FILE_SIZE_BYTES = 5 * 1024 * 1024 * 1024  # 5 GB
SOFT_WARNING_EST_ROWS = 100_000_000  # 100M

# Caution zone: approaching limits — early warning (non-blocking)
CAUTION_FILE_SIZE_BYTES = 1 * 1024 * 1024 * 1024  # 1 GB
CAUTION_EST_ROWS = 10_000_000  # 10M rows


def validate_file_path(file_path: str, must_exist: bool = True) -> Path:
    """
    Validate and resolve a file path securely.

    This function performs several safety checks on a file path:
    1. Makes sure the path is not empty
    2. Checks for security risks (like ".." which could access other folders)
    3. Optionally verifies the file actually exists on the computer

    Think of it as a security guard checking IDs before allowing entry.

    Args:
        file_path: The file location as text (e.g., "C:/data/myfile.csv").
        must_exist: If True, check that the file actually exists on disk.

    Returns:
        A validated Path object that's safe to use.

    Raises:
        FileValidationError: If the path is empty, dangerous, or file not found.
    """
    # CHECK 1: Make sure the path is not empty
    if not file_path:
        raise FileValidationError("File path cannot be empty")

    # CHECK 2: Make sure we received text, not a number or something else
    if not isinstance(file_path, str):
        raise FileValidationError(
            f"File path must be a string, got {type(file_path).__name__}"
        )

    # CHECK 3: Try to convert the text into a proper Path object
    try:
        path = Path(file_path)
    except (TypeError, ValueError) as e:
        raise FileValidationError(f"Invalid file path format: {e}") from e

    # CHECK 4: Resolve the path to its full, absolute form
    # This converts relative paths (like "../data/file.csv") to full paths
    try:
        resolved_path = path.resolve(strict=False)
    except (OSError, RuntimeError) as e:
        raise FileValidationError(f"Cannot resolve file path: {e}") from e

    # CHECK 5: Security check for directory traversal attempts
    # Block any path component exactly equal to ".." (blocks .... and other variants)
    try:
        parts = path.parts
    except (TypeError, ValueError):
        raise FileValidationError("Invalid path structure") from None
    if ".." in parts:
        raise FileValidationError(
            "Directory traversal patterns detected in path"
        )

    # CHECK 6: If requested, verify the file actually exists
    if must_exist:
        if not resolved_path.exists():
            raise FileValidationError(f"File does not exist: {resolved_path}")
        if not resolved_path.is_file():
            raise FileValidationError(f"Path is not a file: {resolved_path}")

    # All checks passed - return the validated path
    return resolved_path


def validate_output_directory(dir_path: str, create: bool = False) -> Path:
    """
    Validate an output directory path (where results will be saved).

    This function checks that we can safely write output files to the
    specified folder. It can optionally create the folder if it doesn't exist.

    Args:
        dir_path: The folder location as text (e.g., "C:/results/").
        create: If True, create the folder if it doesn't exist yet.

    Returns:
        A validated Path object representing the output folder.

    Raises:
        FileValidationError: If the path is invalid or we can't use it.
    """
    # CHECK 1: Make sure the path is not empty
    if not dir_path:
        raise FileValidationError("Output directory path cannot be empty")

    # CHECK 2: Make sure we received text
    if not isinstance(dir_path, str):
        raise FileValidationError(
            f"Directory path must be a string, got {type(dir_path).__name__}"
        )

    # CHECK 3: Convert text to a Path object
    try:
        path = Path(dir_path)
    except (TypeError, ValueError) as e:
        raise FileValidationError(f"Invalid directory path format: {e}") from e

    # CHECK 4: Resolve to full, absolute path
    try:
        resolved_path = path.resolve(strict=False)
    except (OSError, RuntimeError) as e:
        raise FileValidationError(f"Cannot resolve directory path: {e}") from e

    # CHECK 5: Security check for directory traversal
    try:
        parts = path.parts
    except (TypeError, ValueError):
        raise FileValidationError("Invalid path structure") from None
    if ".." in parts:
        raise FileValidationError(
            "Directory traversal patterns detected in path"
        )

    # CHECK 6: Handle the directory existence
    if resolved_path.exists():
        # Path exists - make sure it's actually a folder, not a file
        if not resolved_path.is_dir():
            raise FileValidationError(
                f"Path exists but is not a directory: {resolved_path}"
            )
    elif create:
        # Path doesn't exist but we're allowed to create it
        try:
            # 'parents=True' means create parent folders too if needed
            # 'exist_ok=True' means don't error if it already exists
            resolved_path.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            raise FileValidationError(
                f"Cannot create output directory: {e}"
            ) from e
    else:
        # Path doesn't exist and we're not allowed to create it
        raise FileValidationError(f"Directory does not exist: {resolved_path}")

    # All checks passed
    return resolved_path


def get_file_format(file_path: Path) -> str:
    """
    Determine the file format from its extension (the part after the dot).

    For example:
    - "data.csv" has extension ".csv" (CSV format)
    - "report.xlsx" has extension ".xlsx" (Excel format)

    Args:
        file_path: The Path object pointing to the file.

    Returns:
        The file extension in lowercase (e.g., ".csv", ".json").

    Raises:
        UnsupportedFormatError: If we don't know how to read this file type.
    """
    # Get the file extension (like ".csv" or ".xlsx") and make it lowercase
    extension = file_path.suffix.lower()

    # Check if there's no extension at all
    if not extension:
        raise UnsupportedFormatError(
            f"File has no extension, cannot determine format: {file_path.name}"
        )

    # Legacy Excel (.xls) is not supported
    if extension == ".xls":
        raise UnsupportedFormatError(
            "Legacy Excel (.xls) is not supported. Please convert to .xlsx."
        )

    # Check if we support this file type
    if extension not in SUPPORTED_EXTENSIONS:
        raise UnsupportedFormatError(
            f"Unsupported file format '{extension}'. "
            f"Supported formats: {', '.join(sorted(SUPPORTED_EXTENSIONS))}"
        )

    return extension


def validate_schema(df: pd.DataFrame, required_fields: Set[str]) -> None:
    """
    Validate that a data table contains all required columns.

    A "schema" is like a template that defines what columns should exist
    in the data. This function checks that the actual data matches our
    expected schema (has all the required columns).

    For example, fusion data MUST have 'fusion_id', 'protein_length', and
    'recurrence_count' columns. If any are missing, this function raises an error.

    Args:
        df: The data table (DataFrame) to check.
        required_fields: A set of column names that MUST be present.

    Raises:
        SchemaValidationError: If any required columns are missing.
    """
    # Check if we received no data at all
    if df is None:
        raise SchemaValidationError("DataFrame is None")

    # Check if we received the right type of object
    if not isinstance(df, pd.DataFrame):
        raise SchemaValidationError(
            f"Expected pandas DataFrame, got {type(df).__name__}"
        )

    # Get the list of columns that actually exist in the data
    actual_columns = set(df.columns)
    
    # Find which required columns are missing (required but not in actual)
    missing_fields = required_fields - actual_columns

    # If any columns are missing, report the error
    if missing_fields:
        raise SchemaValidationError(
            f"Missing required fields: {sorted(missing_fields)}. "
            f"Available fields: {sorted(actual_columns)}"
        )

    # Null policy: required columns must not contain nulls (matches fusion_schema.yaml)
    for col in required_fields:
        if col in actual_columns and df[col].isnull().any():
            null_count = int(df[col].isnull().sum())
            raise SchemaValidationError(
                f"Nulls not allowed in required column '{col}'. "
                f"Found {null_count} null value(s)."
            )


def _validate_fusion_numeric_constraints(df: pd.DataFrame) -> None:
    """
    Enforce runtime numeric constraints defined in fusion_schema.yaml
    WITHOUT changing schema loader behavior.

    NaN values are allowed (existing diagnostics decide). Only enforce
    when numeric value exists.
    """
    if "protein_length" in df.columns:
        # NaN <= 0 is False; only 0 and negative trigger
        invalid_mask = df["protein_length"] <= 0
        if invalid_mask.any():
            raise SchemaValidationError(
                "Invalid protein_length values detected (must be > 0)"
            )

    if "recurrence_count" in df.columns:
        # NaN < 0 is False; only negative triggers
        invalid_mask = df["recurrence_count"] < 0
        if invalid_mask.any():
            raise SchemaValidationError(
                "Invalid recurrence_count values detected (must be >= 0)"
            )


# ----- FILE FORMAT LOADERS -----
# Each function below reads a specific file type and returns the data as a table

def _check_file_size_and_csv_rows(file_path: Path, extension: str) -> None:
    """
    Enforce file size and (for CSV/TSV) row estimate guards before load.

    Raises:
        DataLoaderError: If file exceeds MAX_FILE_SIZE_BYTES or estimated rows exceed limit.
    """
    try:
        size = file_path.stat().st_size
    except OSError as e:
        raise DataLoaderError(f"Cannot stat file: {e}") from e

    # Caution zone: approaching limits (early warning)
    if size > CAUTION_FILE_SIZE_BYTES and size <= SOFT_WARNING_FILE_SIZE_BYTES:
        _logger.warning(
            "Input file >1GB: approaching safe limit (5GB). Monitor memory usage."
        )
    # Soft warning: high risk zone
    elif size > SOFT_WARNING_FILE_SIZE_BYTES:
        _logger.warning(
            "Large input file detected (>5GB). Runtime memory pressure possible."
        )

    if size > MAX_FILE_SIZE_BYTES:
        raise DataLoaderError(
            f"File size {size} exceeds maximum allowed {MAX_FILE_SIZE_BYTES} bytes"
        )
    if extension in (".csv", ".tsv"):
        estimated_rows = size // AVG_CSV_ROW_BYTES
        if estimated_rows > CAUTION_EST_ROWS and estimated_rows <= SOFT_WARNING_EST_ROWS:
            _logger.warning(
                "Estimated CSV/TSV rows >10M: approaching safe limit (200M). Monitor performance."
            )
        elif estimated_rows > SOFT_WARNING_EST_ROWS:
            _logger.warning(
                "Very large estimated row count (>100M). Runtime performance may degrade."
            )
        if estimated_rows > MAX_ESTIMATED_CSV_ROWS:
            raise DataLoaderError(
                f"Estimated CSV/TSV rows ({estimated_rows}) exceed limit {MAX_ESTIMATED_CSV_ROWS}"
            )


def load_csv(file_path: Path) -> pd.DataFrame:
    """
    Load a CSV (Comma-Separated Values) file into a data table.

    CSV files are simple text files where each line is a row of data,
    and columns are separated by commas. Example:
        name,age,city
        Alice,30,Boston
        Bob,25,Denver

    Args:
        file_path: Location of the CSV file.

    Returns:
        The loaded data as a pandas DataFrame (like a spreadsheet).

    Raises:
        DataLoaderError: If the file is empty or cannot be read.
    """
    _check_file_size_and_csv_rows(file_path, ".csv")
    try:
        # 'low_memory=False' helps with large files containing mixed data types
        return pd.read_csv(file_path, low_memory=False)
    except pd.errors.EmptyDataError as e:
        raise DataLoaderError(f"CSV file is empty: {file_path.name}") from e
    except pd.errors.ParserError as e:
        raise DataLoaderError(f"CSV parsing error: {e}") from e


def load_tsv(file_path: Path) -> pd.DataFrame:
    """
    Load a TSV (Tab-Separated Values) file into a data table.

    TSV files are like CSV files, but use tabs instead of commas
    to separate columns. Common in bioinformatics data.

    Args:
        file_path: Location of the TSV file.

    Returns:
        The loaded data as a pandas DataFrame.

    Raises:
        DataLoaderError: If the file is empty or cannot be read.
    """
    _check_file_size_and_csv_rows(file_path, ".tsv")
    try:
        # sep="\t" tells pandas to use tab characters as separators
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    except pd.errors.EmptyDataError as e:
        raise DataLoaderError(f"TSV file is empty: {file_path.name}") from e
    except pd.errors.ParserError as e:
        raise DataLoaderError(f"TSV parsing error: {e}") from e


def load_json(file_path: Path) -> pd.DataFrame:
    """
    Load a JSON file into a data table.

    JSON (JavaScript Object Notation) is a structured text format
    commonly used for data exchange. Only tabular JSON structures
    can be converted to a DataFrame.

    Args:
        file_path: Location of the JSON file.

    Returns:
        The loaded data as a pandas DataFrame.

    Raises:
        DataLoaderError: If the file cannot be read or isn't valid JSON.
    """
    try:
        size = file_path.stat().st_size
    except OSError as e:
        raise DataLoaderError(f"Cannot stat file: {e}") from e
    if size > MAX_FILE_SIZE_BYTES:
        raise DataLoaderError(
            f"File size {size} exceeds maximum allowed {MAX_FILE_SIZE_BYTES} bytes"
        )
    try:
        return pd.read_json(file_path)
    except ValueError as e:
        raise DataLoaderError(f"JSON parsing error: {e}") from e


def load_parquet(file_path: Path) -> pd.DataFrame:
    """
    Load a Parquet file into a data table.

    Parquet is an efficient, compressed file format designed for
    large datasets. It's commonly used in big data applications
    because it's faster to read than CSV for large files.

    Args:
        file_path: Location of the Parquet file.

    Returns:
        The loaded data as a pandas DataFrame.

    Raises:
        DataLoaderError: If the file cannot be read.
    """
    try:
        size = file_path.stat().st_size
    except OSError as e:
        raise DataLoaderError(f"Cannot stat file: {e}") from e
    if size > MAX_FILE_SIZE_BYTES:
        raise DataLoaderError(
            f"File size {size} exceeds maximum allowed {MAX_FILE_SIZE_BYTES} bytes"
        )
    try:
        return pd.read_parquet(file_path)
    except Exception as e:
        raise DataLoaderError(f"Parquet reading error: {e}") from e


def load_excel(file_path: Path) -> pd.DataFrame:
    """
    Load an Excel spreadsheet file into a data table.

    This reads Microsoft Excel files (.xlsx or .xls). Only the
    first sheet is loaded. For safety, macros and embedded code
    are not executed.

    Args:
        file_path: Location of the Excel file.

    Returns:
        The loaded data as a pandas DataFrame.

    Raises:
        DataLoaderError: If the file cannot be read.
    """
    try:
        size = file_path.stat().st_size
    except OSError as e:
        raise DataLoaderError(f"Cannot stat file: {e}") from e
    if size > MAX_FILE_SIZE_BYTES:
        raise DataLoaderError(
            f"File size {size} exceeds maximum allowed {MAX_FILE_SIZE_BYTES} bytes"
        )
    try:
        # engine="openpyxl" for .xlsx (legacy .xls not supported)
        return pd.read_excel(file_path, engine="openpyxl")
    except ValueError as e:
        raise DataLoaderError(f"Excel reading error: {e}") from e


def load_data(file_path: str) -> pd.DataFrame:
    """
    Load data from any supported file format into a data table.

    This is the main function to use when you want to load a data file.
    It automatically detects the file type based on the extension and
    uses the appropriate loader.

    Supported formats: CSV, TSV, JSON, Parquet, Excel (.xlsx only)

    Args:
        file_path: Location of the data file (e.g., "C:/data/myfile.csv").

    Returns:
        The loaded data as a pandas DataFrame.

    Raises:
        FileValidationError: If the file path is invalid or file doesn't exist.
        UnsupportedFormatError: If we don't know how to read this file type.
        DataLoaderError: If the file exists but cannot be read properly.
    """
    # STEP 1: Validate the file path and make sure the file exists
    validated_path = validate_file_path(file_path, must_exist=True)
    
    # STEP 2: Determine what type of file it is based on extension
    file_format = get_file_format(validated_path)

    # STEP 3: Map each file extension to its loader function
    # This dictionary tells us which function to use for each file type
    loaders = {
        ".csv": load_csv,       # CSV files
        ".tsv": load_tsv,       # TSV files (tab-separated)
        ".json": load_json,     # JSON files
        ".parquet": load_parquet,  # Parquet files
        ".xlsx": load_excel,    # Excel .xlsx (legacy .xls not supported)
    }

    # STEP 4: Get the appropriate loader function for this file type
    loader_func = loaders.get(file_format)
    if loader_func is None:
        # This shouldn't happen if get_file_format worked correctly
        raise UnsupportedFormatError(f"No loader implemented for {file_format}")

    # STEP 5: Use the loader function to read the file and return the data
    return loader_func(validated_path)


def load_fusion_data(file_path: str) -> pd.DataFrame:
    """
    Load and validate fusion dataset (the main data for analysis).

    This function is specifically for loading the primary fusion gene data.
    It not only loads the file but also checks that all required columns
    are present.

    Required columns:
        - fusion_id: Unique identifier for each fusion event
        - protein_length: Length of the resulting protein
        - recurrence_count: How many times this fusion was observed

    Args:
        file_path: Location of the fusion data file.

    Returns:
        The validated fusion data as a pandas DataFrame.

    Raises:
        FileValidationError: If the file path is invalid.
        UnsupportedFormatError: If the file type is not supported.
        DataLoaderError: If the file cannot be read.
        SchemaValidationError: If required columns are missing.
    """
    # STEP 1: Load the data from the file
    df = load_data(file_path)
    
    # STEP 2: Verify all required columns are present
    validate_schema(df, REQUIRED_FUSION_FIELDS)

    # STEP 3: Enforce numeric constraints (fusion_schema.yaml)
    _validate_fusion_numeric_constraints(df)

    # Return the validated data
    return df


def load_reference_data(file_path: Optional[str]) -> Optional[pd.DataFrame]:
    """
    Load optional reference data (like COSMIC database).

    Reference data is external data used for comparison or validation.
    For example, COSMIC is a database of known cancer mutations.
    
    Unlike fusion data, reference data doesn't need specific columns,
    so no schema validation is performed.

    If no file path is provided, this function simply returns None
    (indicating no reference data is available).

    Args:
        file_path: Location of the reference data file, or None.

    Returns:
        The reference data as a pandas DataFrame, or None if not provided.

    Raises:
        FileValidationError: If the file path is invalid.
        UnsupportedFormatError: If the file type is not supported.
        DataLoaderError: If the file cannot be read.
    """
    # If no path was provided, there's no reference data to load
    if file_path is None:
        return None

    # Load and return the reference data (no schema validation needed)
    return load_data(file_path)