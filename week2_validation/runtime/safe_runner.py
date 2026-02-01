"""
Week 2: Data Integrity & Statistical Validation — Atomic Failure Containment Wrapper.

Executes Data Integrity & Statistical Validation orchestration inside a hardened safety wrapper.
Guarantees: exception classification, no partial output, clean stderr, never masks errors.
"""

import os
import sys
import traceback

from week2_validation.runtime.exit_codes import Week2ExitCode

# Module-scope exception types for classification (avoids import-time fragility in _classify_exception)
try:
    from week2_validation.utils.freeze import FreezeError
    from week2_validation.run_week2 import ConfigurationError, PipelineError
except Exception:
    class _Unimported:
        pass
    FreezeError = _Unimported
    ConfigurationError = _Unimported
    PipelineError = _Unimported


def _classify_exception(exc: BaseException) -> Week2ExitCode:
    """Classify exception to deterministic exit code."""
    # Schema/input validation first (including in cause chain) for consistent exit 10
    try:
        from week2_validation.utils.data_loader import SchemaValidationError

        if isinstance(exc, SchemaValidationError):
            return Week2ExitCode.INPUT_SCHEMA_ERROR
    except ImportError:
        pass

    # Unwrap cause/context for wrapped exceptions
    cause = getattr(exc, "__cause__", None) or getattr(exc, "__context__", None)
    if cause is not None:
        code = _classify_exception(cause)
        if code != Week2ExitCode.UNKNOWN_ERROR:
            return code

    # Freeze failures (explicit type; no message-based detection)
    if isinstance(exc, FreezeError):
        return Week2ExitCode.FREEZE_ERROR

    # Pipeline / configuration errors (explicit types)
    if isinstance(exc, ConfigurationError):
        return Week2ExitCode.CONFIG_ERROR
    if isinstance(exc, PipelineError):
        return Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR

    # Schema / input validation (ValueError, e.g. from adapters)
    if isinstance(exc, ValueError):
        return Week2ExitCode.INPUT_SCHEMA_ERROR

    # Config / YAML errors (loader module)
    try:
        from week2_validation.config.loader import ConfigError

        if isinstance(exc, ConfigError):
            return Week2ExitCode.CONFIG_ERROR
    except ImportError:
        pass
    try:
        import yaml

        if isinstance(exc, yaml.YAMLError):
            return Week2ExitCode.CONFIG_ERROR
    except ImportError:
        pass

    # Other RuntimeError (e.g. diagnostic failures)
    if isinstance(exc, RuntimeError):
        return Week2ExitCode.DIAGNOSTIC_RUNTIME_ERROR

    return Week2ExitCode.UNKNOWN_ERROR


def run_week2_safely(main_callable, *args, **kwargs) -> tuple[int, str]:
    """
    Execute Data Integrity & Statistical Validation orchestration inside a hardened safety wrapper.

    Guarantees:
    - Exception classification → deterministic exit code
    - No partial machine output
    - Clean stderr messaging
    - Never masks freeze or schema errors

    Args:
        main_callable: The main Data Integrity & Statistical Validation entry point (e.g., run_pipeline).
        *args: Positional arguments passed to main_callable.
        **kwargs: Keyword arguments passed to main_callable.

    Returns:
        Tuple of (exit_code: int, status_message: str).
    """
    try:
        exit_code = main_callable(*args, **kwargs)
        if exit_code is None:
            exit_code = 0
        return (int(exit_code), "OK" if exit_code == 0 else f"Exit {exit_code}")
    except BaseException as exc:
        code = _classify_exception(exc)
        msg = str(exc)
        if not msg:
            msg = type(exc).__name__

        if os.environ.get("DEBUG", "").strip().lower() in ("1", "true", "yes"):
            traceback.print_exc(file=sys.stderr)

        print(f"Week2 error [{code.name}]: {msg}", file=sys.stderr)
        return (int(code), msg)
