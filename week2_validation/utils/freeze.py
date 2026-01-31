"""
freeze.py - Week 2 Defensive Dataset Freeze Logic.

This module provides Week 2's own dataset freeze mechanism that ensures
data immutability BEFORE any diagnostics run. This is defensive freezing,
independent of any Week 1 behavior.

WHAT IS DEFENSIVE FREEZING?
Week 2 cannot assume that incoming data was frozen by Week 1. To ensure
data integrity and reproducibility, Week 2 enforces its own freeze logic:

1. Detect if incoming data is already frozen (has manifest + marker)
2. If not frozen, create a frozen snapshot owned by Week 2
3. All diagnostics operate on the frozen snapshot, never raw input

WHY IS THIS NEEDED?
- Week 1 behavior is unknown and cannot be assumed
- Week 2 must protect itself from data changes during analysis
- Audit trail requires explicit proof of data state at analysis time
- Reproducibility requires immutable data snapshots

FREEZE ARTIFACTS:
- frozen_inputs/<hash_prefix>/  -- frozen dataset directory
- frozen_inputs/<hash_prefix>/manifest.json  -- metadata + hash
- frozen_inputs/<hash_prefix>/.frozen  -- marker file
- frozen_inputs/<hash_prefix>/<original_filename>  -- copied data
"""

import hashlib
import json
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

# Import version from package
try:
    from week2_validation import __version__ as WEEK2_VERSION
except ImportError:
    WEEK2_VERSION = "unknown"


# =============================================================================
# CONSTANTS
# =============================================================================

# Names of freeze artifacts
MANIFEST_FILENAME = "manifest.json"
FROZEN_MARKER_FILENAME = ".frozen"

# Hash prefix length for directory naming (first 12 chars of SHA-256)
HASH_PREFIX_LENGTH = 12

# Buffer size for file hashing (64KB)
HASH_BUFFER_SIZE = 65536


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class FreezeManifest:
    """
    Immutable manifest for a frozen dataset.
    
    Contains all metadata needed to verify and audit the freeze.
    """
    source_path: str
    frozen_by: str
    frozen_at_utc: str
    file_hash_sha256: str
    week2_version: str
    original_filename: str


# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

def _compute_file_hash(file_path: Path) -> str:
    """
    Compute SHA-256 hash of a file.
    
    Args:
        file_path: Path to the file to hash.
    
    Returns:
        Hexadecimal SHA-256 hash string.
    
    Raises:
        RuntimeError: If file cannot be read or hashed.
    """
    try:
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(HASH_BUFFER_SIZE), b""):
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()
    except (IOError, OSError) as e:
        raise RuntimeError(f"Failed to compute hash for {file_path}: {e}") from e


def _load_manifest(manifest_path: Path) -> Optional[FreezeManifest]:
    """
    Load and parse a freeze manifest file.
    
    Args:
        manifest_path: Path to the manifest.json file.
    
    Returns:
        FreezeManifest if valid, None if file doesn't exist or is invalid.
    """
    if not manifest_path.exists():
        return None
    
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        return FreezeManifest(
            source_path=data["source_path"],
            frozen_by=data["frozen_by"],
            frozen_at_utc=data["frozen_at_utc"],
            file_hash_sha256=data["file_hash_sha256"],
            week2_version=data["week2_version"],
            original_filename=data.get("original_filename", ""),
        )
    except (json.JSONDecodeError, KeyError, IOError):
        return None


def _write_manifest(manifest: FreezeManifest, manifest_path: Path) -> None:
    """
    Write a freeze manifest to disk.
    
    Args:
        manifest: The FreezeManifest to write.
        manifest_path: Path where manifest.json will be written.
    
    Raises:
        RuntimeError: If manifest cannot be written.
    """
    manifest_dict = {
        "source_path": manifest.source_path,
        "frozen_by": manifest.frozen_by,
        "frozen_at_utc": manifest.frozen_at_utc,
        "file_hash_sha256": manifest.file_hash_sha256,
        "week2_version": manifest.week2_version,
        "original_filename": manifest.original_filename,
    }
    
    try:
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest_dict, f, indent=2)
    except (IOError, OSError) as e:
        raise RuntimeError(f"Failed to write manifest: {e}") from e


def _find_existing_frozen_dir(
    input_path: Path,
    frozen_root: Path,
) -> Optional[Path]:
    """
    Find an existing frozen directory for the given input file.
    
    Searches frozen_root for directories containing a valid manifest
    that matches the current input file's hash.
    
    Args:
        input_path: Path to the input file.
        frozen_root: Root directory for frozen datasets.
    
    Returns:
        Path to existing frozen directory if found and valid, None otherwise.
    """
    if not frozen_root.exists():
        return None
    
    # Compute hash of current input file
    try:
        current_hash = _compute_file_hash(input_path)
    except RuntimeError:
        return None
    
    # Check all subdirectories for matching manifest
    for subdir in frozen_root.iterdir():
        if not subdir.is_dir():
            continue
        
        manifest_path = subdir / MANIFEST_FILENAME
        marker_path = subdir / FROZEN_MARKER_FILENAME
        
        # Must have both manifest and marker
        if not manifest_path.exists() or not marker_path.exists():
            continue
        
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            continue
        
        # Check if hash matches
        if manifest.file_hash_sha256 == current_hash:
            # Verify the frozen file still exists and matches
            frozen_file_path = subdir / manifest.original_filename
            if frozen_file_path.exists():
                try:
                    frozen_hash = _compute_file_hash(frozen_file_path)
                    if frozen_hash == current_hash:
                        return subdir
                except RuntimeError:
                    continue
    
    return None


def _validate_frozen_snapshot(frozen_dir: Path) -> bool:
    """
    Validate that a frozen snapshot is complete and consistent.
    
    Args:
        frozen_dir: Path to the frozen directory.
    
    Returns:
        True if snapshot is valid, False otherwise.
    """
    manifest_path = frozen_dir / MANIFEST_FILENAME
    marker_path = frozen_dir / FROZEN_MARKER_FILENAME
    
    # Check required files exist
    if not manifest_path.exists():
        return False
    if not marker_path.exists():
        return False
    
    # Load and validate manifest
    manifest = _load_manifest(manifest_path)
    if manifest is None:
        return False
    
    # Check frozen file exists
    frozen_file_path = frozen_dir / manifest.original_filename
    if not frozen_file_path.exists():
        return False
    
    # Verify hash matches
    try:
        actual_hash = _compute_file_hash(frozen_file_path)
        return actual_hash == manifest.file_hash_sha256
    except RuntimeError:
        return False


# =============================================================================
# PUBLIC API
# =============================================================================

def ensure_frozen_input(
    input_path: Path,
    frozen_root: Path,
) -> Path:
    """
    Ensure that the input dataset is frozen before Week 2 diagnostics.

    If the dataset is already frozen (manifest + .frozen marker),
    validate and return frozen path.

    If not frozen, create a frozen snapshot owned by Week 2.

    This function is:
    - Deterministic: Same input always produces same frozen directory
    - Idempotent: Running twice does NOT create duplicate freezes
    - Auditable: All freeze events are logged with full metadata

    Args:
        input_path: Path to the input dataset file.
        frozen_root: Root directory for storing frozen snapshots.

    Returns:
        Path to frozen dataset directory to be used by diagnostics.

    Raises:
        RuntimeError: If freeze fails. Diagnostics must NOT proceed.
    """
    # =========================================================================
    # STEP 1: Validate input path
    # =========================================================================
    if not input_path.exists():
        raise RuntimeError(f"Input file does not exist: {input_path}")
    
    if not input_path.is_file():
        raise RuntimeError(f"Input path is not a file: {input_path}")
    
    # =========================================================================
    # STEP 2: Check for existing frozen snapshot
    # =========================================================================
    existing_frozen_dir = _find_existing_frozen_dir(input_path, frozen_root)
    
    if existing_frozen_dir is not None:
        # Validate the existing snapshot
        if _validate_frozen_snapshot(existing_frozen_dir):
            print("Dataset already frozen — using existing frozen snapshot")
            return existing_frozen_dir
        else:
            # Existing snapshot is corrupted, will create new one
            print("WARNING: Existing frozen snapshot is invalid, creating new freeze")
    
    # =========================================================================
    # STEP 3: Dataset not frozen — create frozen snapshot
    # =========================================================================
    print("Dataset not frozen — freezing input for Week 2 diagnostics")
    
    # Compute hash of input file
    try:
        file_hash = _compute_file_hash(input_path)
    except RuntimeError as e:
        raise RuntimeError(f"ERROR: Failed to freeze dataset — diagnostics aborted: {e}") from e
    
    # Create frozen directory using hash prefix for deterministic naming
    hash_prefix = file_hash[:HASH_PREFIX_LENGTH]
    frozen_dir = frozen_root / hash_prefix
    
    # Create frozen_root if it doesn't exist
    try:
        frozen_root.mkdir(parents=True, exist_ok=True)
    except (IOError, OSError) as e:
        raise RuntimeError(f"ERROR: Failed to freeze dataset — diagnostics aborted: {e}") from e
    
    # Check if directory already exists (idempotency check)
    if frozen_dir.exists():
        if _validate_frozen_snapshot(frozen_dir):
            print("Dataset already frozen — using existing frozen snapshot")
            return frozen_dir
        else:
            # Remove invalid frozen directory
            try:
                shutil.rmtree(frozen_dir)
            except (IOError, OSError) as e:
                raise RuntimeError(f"ERROR: Failed to freeze dataset — diagnostics aborted: {e}") from e
    
    # Create frozen directory
    try:
        frozen_dir.mkdir(parents=True, exist_ok=True)
    except (IOError, OSError) as e:
        raise RuntimeError(f"ERROR: Failed to freeze dataset — diagnostics aborted: {e}") from e
    
    # =========================================================================
    # STEP 4: Copy input file to frozen directory
    # =========================================================================
    original_filename = input_path.name
    frozen_file_path = frozen_dir / original_filename
    
    try:
        shutil.copy2(input_path, frozen_file_path)
    except (IOError, OSError) as e:
        # Clean up partial freeze
        try:
            shutil.rmtree(frozen_dir)
        except:
            pass
        raise RuntimeError(f"ERROR: Failed to freeze dataset — diagnostics aborted: {e}") from e
    
    # Verify copy integrity
    try:
        copied_hash = _compute_file_hash(frozen_file_path)
        if copied_hash != file_hash:
            shutil.rmtree(frozen_dir)
            raise RuntimeError("ERROR: Failed to freeze dataset — diagnostics aborted: Copy verification failed")
    except RuntimeError:
        try:
            shutil.rmtree(frozen_dir)
        except:
            pass
        raise
    
    # =========================================================================
    # STEP 5: Write manifest
    # =========================================================================
    manifest = FreezeManifest(
        source_path=str(input_path.resolve()),
        frozen_by="week2_validation",
        frozen_at_utc=datetime.now(timezone.utc).isoformat(),
        file_hash_sha256=file_hash,
        week2_version=WEEK2_VERSION,
        original_filename=original_filename,
    )
    
    manifest_path = frozen_dir / MANIFEST_FILENAME
    try:
        _write_manifest(manifest, manifest_path)
    except RuntimeError as e:
        # Clean up partial freeze
        try:
            shutil.rmtree(frozen_dir)
        except:
            pass
        raise RuntimeError(f"ERROR: Failed to freeze dataset — diagnostics aborted: {e}") from e
    
    # =========================================================================
    # STEP 6: Create .frozen marker file
    # =========================================================================
    marker_path = frozen_dir / FROZEN_MARKER_FILENAME
    try:
        marker_path.touch()
    except (IOError, OSError) as e:
        # Clean up partial freeze
        try:
            shutil.rmtree(frozen_dir)
        except:
            pass
        raise RuntimeError(f"ERROR: Failed to freeze dataset — diagnostics aborted: {e}") from e
    
    # =========================================================================
    # STEP 7: Final validation and success
    # =========================================================================
    if not _validate_frozen_snapshot(frozen_dir):
        try:
            shutil.rmtree(frozen_dir)
        except:
            pass
        raise RuntimeError("ERROR: Failed to freeze dataset — diagnostics aborted: Final validation failed")
    
    print("Freeze complete — diagnostics may now proceed")
    
    return frozen_dir


def get_frozen_data_path(frozen_dir: Path) -> Path:
    """
    Get the path to the frozen data file within a frozen directory.
    
    Args:
        frozen_dir: Path to the frozen directory.
    
    Returns:
        Path to the frozen data file.
    
    Raises:
        RuntimeError: If frozen directory is invalid or data file not found.
    """
    manifest_path = frozen_dir / MANIFEST_FILENAME
    manifest = _load_manifest(manifest_path)
    
    if manifest is None:
        raise RuntimeError(f"Cannot read manifest from frozen directory: {frozen_dir}")
    
    frozen_file_path = frozen_dir / manifest.original_filename
    
    if not frozen_file_path.exists():
        raise RuntimeError(f"Frozen data file not found: {frozen_file_path}")
    
    return frozen_file_path
