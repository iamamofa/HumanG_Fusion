"""
state.py - Dataset freeze state enforcement.

Enforces dataset freeze using:
  1. Flag file (authoritative)
  2. CLI argument (secondary)

Flag file takes precedence. Fails fast on invalid or indeterminate states.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


class StateError(Exception):
    """Raised when state validation fails or cannot be determined."""


@dataclass(frozen=True)
class FreezeState:
    """Immutable representation of dataset freeze state."""

    is_frozen: bool
    source: str


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
    if flag_path is None:
        raise StateError("Flag path cannot be None")

    if not isinstance(flag_path, Path):
        raise StateError("Flag path must be a Path object")

    try:
        if flag_path.exists():
            if not flag_path.is_file():
                raise StateError("Flag path exists but is not a regular file")
            return True
        return False
    except PermissionError as exc:
        raise StateError("Permission denied when accessing flag path") from exc
    except OSError as exc:
        raise StateError(f"OS error when checking flag path: {exc}") from exc


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

    Args:
        flag_path: Path to the freeze flag file.
        cli_frozen: Optional CLI argument indicating freeze state.

    Returns:
        FreezeState with is_frozen status and source indicator.

    Raises:
        StateError: If state cannot be determined.
    """
    flag_exists = _check_flag_file(flag_path)

    if flag_exists:
        return FreezeState(is_frozen=True, source="flag_file")

    if cli_frozen is not None:
        return FreezeState(is_frozen=cli_frozen, source="cli_argument")

    raise StateError(
        "Cannot determine freeze state: "
        "flag file does not exist and no CLI argument provided"
    )


def enforce_frozen(state: FreezeState) -> None:
    """
    Enforce that the dataset is frozen.

    Args:
        state: The resolved freeze state.

    Raises:
        StateError: If the dataset is not frozen.
    """
    if not isinstance(state, FreezeState):
        raise StateError("Invalid state object")

    if not state.is_frozen:
        raise StateError(
            f"Dataset is not frozen (source: {state.source}). "
            "Pipeline requires frozen dataset to proceed."
        )