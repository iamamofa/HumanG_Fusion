"""
Week 2: Data Integrity & Statistical Validation — Logging configuration.

Provides run-scoped logging with optional file output for audit trail.
Supports both text and JSON logging formats.
"""
from __future__ import annotations

import json
import logging
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path


_RUN_ID: str | None = None


def get_run_id() -> str:
    global _RUN_ID
    if _RUN_ID is None:
        _RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8]
    return _RUN_ID


class JSONFormatter(logging.Formatter):
    """JSON formatter for structured logging."""
    
    def format(self, record: logging.LogRecord) -> str:
        log_entry = {
            "timestamp": datetime.fromtimestamp(record.created, tz=timezone.utc).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }
        if record.exc_info:
            log_entry["exception"] = self.formatException(record.exc_info)
        return json.dumps(log_entry)


def setup_logging(
    output_dir: Path | None = None,
    level: str = "INFO",
    json_format: bool = False,
) -> logging.Logger:
    """
    Setup logging with optional JSON format.
    
    Args:
        output_dir: Directory for log file (if None, no file handler).
        level: Logging level (INFO, DEBUG, WARNING, ERROR).
        json_format: If True, use JSON format; otherwise use text format.
    """
    logger = logging.getLogger("week2_validation")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    if json_format:
        fmt = JSONFormatter()
    else:
        fmt = logging.Formatter(
            fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S",
        )

    # Console handler (always text format for readability)
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    ))
    logger.addHandler(ch)

    # File handler (if output_dir provided)
    if output_dir is not None:
        log_filename = f"validation_{get_run_id()}.log"
        if json_format:
            log_filename = log_filename.replace(".log", ".jsonl")
        fh = logging.FileHandler(output_dir / log_filename, encoding="utf-8")
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger
