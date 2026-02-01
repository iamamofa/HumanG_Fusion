"""
Week 2: Data Integrity & Statistical Validation — Machine Status Envelope.

Extends (not replaces) JSON summary with runtime metadata.
Sanitizes notes to avoid leaking absolute paths, stack traces, or user home.
"""

import re
from datetime import datetime, timezone
from typing import Optional

try:
    from week2_validation import __version__ as WEEK2_VERSION
except ImportError:
    WEEK2_VERSION = "unknown"


def _sanitize_note(note: str) -> str:
    """
    Remove absolute paths, stack traces, environment paths, user home from a note.
    Allow error class name, short message, error category.
    Covers UNC paths, symlinks, long paths, and common path patterns.
    """
    if not note or not isinstance(note, str):
        return ""
    s = note.strip()
    # Remove lines that look like stack trace (File "...", line N)
    s = re.sub(r'\n\s*File "[^"]+", line \d+.*', "", s)
    s = re.sub(r'\n\s*File \'[^\']+\', line \d+.*', "", s)
    # Windows UNC (\\?\ or \\server\share)
    s = re.sub(r"\\\\\?\\[^\s]*", "<path>", s)
    s = re.sub(r"\\\\[^\s\\]+\\[^\s]*", "<path>", s)
    # Windows long path (\\?\C:\...)
    s = re.sub(r"\\\\\?\\[A-Za-z]:\\[^\s]*", "<path>", s, flags=re.IGNORECASE)
    # Absolute paths (Unix and Windows)
    s = re.sub(r"/[A-Za-z]:?[\w/\\\.\-]*", "<path>", s)
    s = re.sub(r"[A-Za-z]:\\[\w\\\.\-]*", "<path>", s)
    # User home patterns (/Users/name, /home/name, C:\Users\name)
    s = re.sub(r"(/Users|/home)/[^\s]*", "<path>", s)
    s = re.sub(r"[A-Za-z]:\\Users\\[^\s]*", "<path>", s, flags=re.IGNORECASE)
    # Symlink targets (path -> path)
    s = re.sub(r"(?:[/\\][\w/\\\.\-]+|[\w]:\\[\w\\\.\-]*)\s*->\s*[\w/\\\.\-]+", "<path>", s)
    # Multiple slash collapse (e.g. /// or \\\\) so no path-like run remains
    s = re.sub(r"/{2,}", "/", s)
    s = re.sub(r"\\\\+", r"\\", s)
    # Collapse repeated <path>
    s = re.sub(r"(<path>\s*)+", "<path> ", s)
    return s.strip()[:500]  # Cap length


def _sanitize_notes(notes: list[str]) -> list[str]:
    """Sanitize each note; drop empty after sanitization."""
    out = []
    for n in notes:
        if not isinstance(n, str):
            continue
        sanitized = _sanitize_note(n)
        if sanitized:
            out.append(sanitized)
    return out


def build_status_envelope(
    dataset_hash: str,
    diagnostics_run: list[str],
    exit_code: int,
    notes: list[str],
    cosmic_reference_loaded: Optional[bool] = None,
    cosmic_reference_requested: Optional[bool] = None,
    scipy_available: Optional[bool] = None,
    psutil_available: Optional[bool] = None,
    data_quality: Optional[dict] = None,
    diagnostics_skipped: Optional[bool] = None,
    skip_reason: Optional[str] = None,
) -> dict:
    """
    Build machine-readable status envelope with runtime metadata.

    Args:
        dataset_hash: Hash identifying the frozen dataset.
        diagnostics_run: List of diagnostic names executed.
        exit_code: Pipeline exit code (0 = success).
        notes: Optional notes or messages (sanitized: no paths/stack traces).
        cosmic_reference_loaded: If --cosmic-data was provided, True when load
            succeeded, False when load failed. Omitted when cosmic not requested.
        cosmic_reference_requested: True when --cosmic-data was provided.
        scipy_available: True if scipy import succeeded at runtime (observability only).
        psutil_available: True if psutil import succeeded at runtime (observability only).

    Returns:
        Dict with week2_version, run_timestamp_utc, dataset_hash,
        diagnostics_run, exit_code, status, notes.
        Optional keys added when provided. JSON-serializable.
    """
    if exit_code == 0:
        if data_quality and data_quality.get("high_risk_flag"):
            status = "SUCCESS_WITH_HIGH_DATA_RISK"
        elif data_quality and data_quality.get("warning_flag"):
            status = "SUCCESS_WITH_WARNINGS"
        else:
            status = "SUCCESS"
    elif exit_code < 50:
        status = "FAILED"
    else:
        status = "PARTIAL"

    envelope = {
        "week2_version": WEEK2_VERSION,
        "run_timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "dataset_hash": dataset_hash,
        "diagnostics_run": list(diagnostics_run),
        "exit_code": int(exit_code),
        "status": status,
        "notes": _sanitize_notes(list(notes)),
    }
    if cosmic_reference_loaded is not None:
        envelope["cosmic_reference_loaded"] = bool(cosmic_reference_loaded)
    if cosmic_reference_requested is not None:
        envelope["cosmic_reference_requested"] = bool(cosmic_reference_requested)
    if scipy_available is not None:
        envelope["scipy_available"] = bool(scipy_available)
    if psutil_available is not None:
        envelope["psutil_available"] = bool(psutil_available)
    if data_quality is not None:
        envelope["data_quality"] = dict(data_quality)
    if diagnostics_skipped:
        envelope["diagnostics"] = {
            "skipped": True,
            "skip_reason": skip_reason or "EMPTY_DATASET",
        }
    return envelope
