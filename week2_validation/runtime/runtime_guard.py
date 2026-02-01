"""
Week 2: Data Integrity & Statistical Validation — Runtime Execution Guard.

Lightweight runtime monitoring: prevents runaway execution, unbounded memory.
Non-invasive. No thread spawning. No perf penalty.
"""

import time

MAX_RUNTIME_SECONDS = 3600
MAX_MEMORY_MB = 8192

_guard_start_time: float | None = None

try:
    import psutil

    _PSUTIL_AVAILABLE = True
except ImportError:
    _PSUTIL_AVAILABLE = False


def start_runtime_guard() -> None:
    """
    Start lightweight runtime monitoring.

    Records start time for duration checks. Call check_runtime_guard()
    periodically during execution to enforce limits.
    Non-invasive. No thread spawning. No perf penalty.
    """
    global _guard_start_time
    _guard_start_time = time.monotonic()


def check_runtime_guard() -> None:
    """
    Verify runtime limits. Call periodically during execution.

    Raises:
        RuntimeError: If MAX_RUNTIME_SECONDS or MAX_MEMORY_MB exceeded.
    """
    if _guard_start_time is None:
        return

    elapsed = time.monotonic() - _guard_start_time
    if elapsed > MAX_RUNTIME_SECONDS:
        raise RuntimeError("Runtime guard triggered")

    if _PSUTIL_AVAILABLE:
        try:
            process = psutil.Process()
            mem_mb = process.memory_info().rss / (1024 * 1024)
            if mem_mb > MAX_MEMORY_MB:
                raise RuntimeError("Runtime guard triggered")
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
