"""
Week 2: Data Integrity & Statistical Validation — Deterministic Exit Code Contract.

Defines machine-readable exit codes for pipeline integration.
No dependencies outside stdlib.
"""

from enum import IntEnum


class Week2ExitCode(IntEnum):
    """Deterministic exit codes for Data Integrity & Statistical Validation pipeline."""

    SUCCESS = 0
    INPUT_SCHEMA_ERROR = 10
    FREEZE_ERROR = 20
    QUALITY_GATE_REJECTED = 25  # One or more quality gates FAIL
    DIAGNOSTIC_RUNTIME_ERROR = 30
    CONFIG_ERROR = 40
    UNKNOWN_ERROR = 99


__all__ = ["Week2ExitCode"]
