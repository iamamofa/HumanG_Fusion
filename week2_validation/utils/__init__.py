"""
Utility modules for Week 2 validation pipeline.

This package provides core utilities for data loading, validation,
and secure file handling operations.

WHAT IS THIS FILE?
In Python, an __init__.py file marks a folder as a "package" (a collection
of related code files). This particular file makes it easy to import
commonly-used functions from the utils folder.

Instead of writing:
    from week2_validation.utils.data_loader import load_fusion_data

You can simply write:
    from week2_validation.utils import load_fusion_data

This is a convenience feature that makes the code cleaner and easier to use.
"""

# Import all the commonly-used items from the data_loader module
# This makes them available directly from 'week2_validation.utils'
from week2_validation.utils.data_loader import (
    DataLoaderError,           # Error when data can't be read
    FileValidationError,       # Error when file path is invalid
    SchemaValidationError,     # Error when required columns are missing
    UnsupportedFormatError,    # Error when file type isn't supported
    load_data,                 # Load any supported data file
    load_fusion_data,          # Load and validate fusion gene data
    load_reference_data,       # Load optional reference data (like COSMIC)
    validate_file_path,        # Check if a file path is valid and safe
    validate_output_directory, # Check if an output folder is valid
    validate_schema,           # Check if data has required columns
)

# __all__ defines what gets exported when someone writes 'from utils import *'
# This is a list of all the public items available from this package
__all__ = [
    "DataLoaderError",
    "FileValidationError",
    "SchemaValidationError",
    "UnsupportedFormatError",
    "load_data",
    "load_fusion_data",
    "load_reference_data",
    "validate_file_path",
    "validate_output_directory",
    "validate_schema",
]