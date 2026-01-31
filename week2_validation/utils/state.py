"""
state.py - Dataset freeze state enforcement.

Enforces dataset freeze using:
  1. Flag file (authoritative)
  2. CLI argument (secondary)

Flag file takes precedence. Fails fast on invalid or indeterminate states.

WHAT IS "FREEZE STATE"?
In data validation pipelines, "freezing" a dataset means locking it down
so no changes can be made. This ensures that all analysis is performed
on the exact same data. Think of it like taking a snapshot.

WHY IS THIS IMPORTANT?
- Reproducibility: Everyone analyzes the same data
- Audit trail: We know exactly what data was validated
- Safety: Prevents accidental modifications during analysis

HOW DOES THIS WORK?
1. If a "flag file" exists (e.g., .frozen), the dataset is frozen
2. If no flag file exists, we check the CLI argument
3. If neither exists, we cannot determine the state (error)
"""

# 'dataclass' creates a simple class for holding data
from dataclasses import dataclass
# 'Path' helps work with file and folder locations
from pathlib import Path
# 'Optional' means a value might be present or might be None
from typing import Optional


# =============================================================================
# CUSTOM ERROR TYPE
# =============================================================================

class StateError(Exception):
    """
    Raised when state validation fails or cannot be determined.
    
    This error indicates that something is wrong with the freeze state:
    - Flag file cannot be accessed
    - State is ambiguous or undefined
    - Dataset is not frozen when it should be
    """


# =============================================================================
# FREEZE STATE DATA STRUCTURE
# =============================================================================

# The '@dataclass' decorator creates a simple class for holding data.
# 'frozen=True' means once created, the values cannot be changed.
@dataclass(frozen=True)
class FreezeState:
    """
    Immutable representation of dataset freeze state.
    
    This class holds two pieces of information:
    - is_frozen: Whether the dataset is currently frozen (True/False)
    - source: Where we got this information from ("flag_file" or "cli_argument")
    """

    is_frozen: bool   # True if dataset is frozen, False otherwise
    source: str       # Where this information came from


# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

def _check_flag_file(flag_path: Path) -> bool:
    """
    Check if freeze flag file exists.

    Args:
        flag_path: Path to the freeze flag file.

    Returns:
        True if flag file exists, False otherwise.

    Raises:
        StateError: If path is invalid, inaccessible, or not a regular file.
    """
    # CHECK 1: Make sure we received a valid path (not None)
    if flag_path is None:
        raise StateError("Flag path cannot be None")

    # CHECK 2: Make sure it's a Path object, not a string or something else
    if not isinstance(flag_path, Path):
        raise StateError("Flag path must be a Path object")

    # CHECK 3: Try to check if the file exists on the filesystem
    try:
        if flag_path.exists():
            # File exists - but make sure it's a regular file, not a directory
            if not flag_path.is_file():
                raise StateError("Flag path exists but is not a regular file")
            # Flag file exists, dataset is frozen
            return True
        # Flag file does not exist, dataset is NOT frozen (from flag file perspective)
        return False
    except PermissionError as exc:
        # Cannot access the file due to permissions
        raise StateError("Permission denied when accessing flag path") from exc
    except OSError as exc:
        # Some other filesystem error occurred
        raise StateError(f"OS error when checking flag path: {exc}") from exc


# =============================================================================
# MAIN PUBLIC FUNCTIONS
# =============================================================================

def resolve_freeze_state(
    flag_path: Path,
    cli_frozen: Optional[bool] = None,
) -> FreezeState:
    """
    Resolve the dataset freeze state from flag file and CLI argument.

    The flag file is authoritative. If the flag file exists, the dataset
    is considered frozen regardless of the CLI argument. If the flag file
    does not exist, the CLI argument is used. If neither is available,
    state cannot be determined and an error is raised.

    PRIORITY ORDER:
    1. Flag file (highest priority - if it exists, dataset is frozen)
    2. CLI argument (secondary - used if flag file doesn't exist)
    3. Error (if neither source provides information)

    Args:
        flag_path: Path to the freeze flag file.
        cli_frozen: Optional CLI argument indicating freeze state.

    Returns:
        FreezeState with is_frozen status and source indicator.

    Raises:
        StateError: If state cannot be determined.
    """
    # STEP 1: Check if the flag file exists
    # This is the authoritative source - if it exists, dataset is frozen
    flag_exists = _check_flag_file(flag_path)

    # PRIORITY 1: Flag file takes precedence over everything else
    if flag_exists:
        # Flag file exists, so dataset is definitely frozen
        return FreezeState(is_frozen=True, source="flag_file")

    # PRIORITY 2: Fall back to CLI argument if flag file doesn't exist
    if cli_frozen is not None:
        # Use the CLI argument's value (could be True or False)
        return FreezeState(is_frozen=cli_frozen, source="cli_argument")

    # PRIORITY 3: Neither source available - cannot determine state
    # This is an error condition - we MUST know the freeze state
    raise StateError(
        "Cannot determine freeze state: "
        "flag file does not exist and no CLI argument provided"
    )


def enforce_frozen(state: FreezeState) -> None:
    """
    Enforce that the dataset is frozen.

    This function acts as a "gate" - it only allows the pipeline to proceed
    if the dataset is frozen. Think of it as a security checkpoint.

    Args:
        state: The resolved freeze state.

    Raises:
        StateError: If the dataset is not frozen.
    """
    # CHECK 1: Make sure we received a valid FreezeState object
    if not isinstance(state, FreezeState):
        raise StateError("Invalid state object")

    # CHECK 2: Verify the dataset is actually frozen
    # If not frozen, raise an error to stop the pipeline
    if not state.is_frozen:
        raise StateError(
            f"Dataset is not frozen (source: {state.source}). "
            "Pipeline requires frozen dataset to proceed."
        )
    
    # If we get here, the dataset is frozen - pipeline can proceed