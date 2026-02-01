"""
Data Integrity & Statistical Validation Pipeline — validates authenticity and reliability of fusion protein length data.

This package provides tools for validating cancer genomics fusion
datasets as part of a reproducible research pipeline.

WHAT IS THIS PROJECT?
This is a data validation pipeline for analyzing cancer genomics fusion data.
"Fusion genes" are abnormal genes formed when parts of two different genes
combine - they are often found in cancer cells.

This pipeline checks that fusion data is:
1. Properly formatted (has the right columns)
2. Statistically reasonable (passes Benford's Law checks, etc.)
3. Ready for further analysis

This pipeline is diagnostic only:
    - No hypothesis testing (not making scientific claims)
    - No modeling (not predicting anything)
    - No biological claims (just checking data quality)

WHAT'S IN THIS PACKAGE?
Modules:
    run_week2: The main program you run from the command line.
    utils: Helper functions for loading and validating files.
    benford: Tools for Benford's Law analysis (digit distribution checks).
"""

# Version number for this package
# Following semantic versioning: MAJOR.MINOR.PATCH
# 0.1.0 means early development, first minor version
__version__ = "0.1.0"