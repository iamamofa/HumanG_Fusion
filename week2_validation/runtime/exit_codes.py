"""
Week 2 Deterministic Exit Code Contract.

Defines machine-readable exit codes for pipeline integration.
No dependencies outside stdlib.
"""

from enum import IntEnum


class Week2ExitCode(IntEnum):
    """Deterministic exit codes for Week 2 validation pipeline."""

    SUCCESS = 0
    INPUT_SCHEMA_ERROR = 10
    FREEZE_ERROR = 20
    DIAGNOSTIC_RUNTIME_ERROR = 30
    CONFIG_ERROR = 40
    UNKNOWN_ERROR = 99


__all__ = ["Week2ExitCode"]
