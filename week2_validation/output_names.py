"""
Data Integrity & Statistical Validation — output filename constants.

The only "week2" in this project is the folder week2_validation.
All output files use validation_* or data_integrity_* naming.
"""

# Status and diagnostic outputs
STATUS_FILENAME_PATTERN = "validation_status_{stem}.json"
DIAGNOSTIC_RESULTS_FILENAME_PATTERN = "validation_diagnostic_results_{stem}.json"
QUALITY_GATES_FILENAME_PATTERN = "validation_quality_gates_{stem}.json"
CERTIFICATION_FILENAME_PATTERN = "validation_dataset_certification_{stem}.json"
CLEANED_DATASET_FILENAME_PATTERN = "validation_cleaned_dataset_{stem}.csv"

# Intermediate / adapted outputs
ADAPTED_FUSION_FILENAME = "validation_adapted_fusion.csv"
FROZEN_FILENAME_PATTERN = "validation_frozen_{stem}.csv"

# Reports
DATA_INTEGRITY_REPORT_PDF = "data_integrity_validation_report.pdf"
DATA_INTEGRITY_REPORT_HTML = "data_integrity_validation_report.html"

# JSON envelope key (replaces week2_version)
VALIDATION_VERSION_KEY = "validation_version"

# Legacy: for backward compatibility when reading old outputs
LEGACY_STATUS_PREFIX = "week2_status_"
LEGACY_DIAGNOSTIC_PREFIX = "week2_diagnostic_results_"
LEGACY_CERT_PREFIX = "week2_dataset_certification_"
LEGACY_CLEANED_PREFIX = "week2_cleaned_dataset_"
