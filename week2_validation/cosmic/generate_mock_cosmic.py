"""
Mock COSMIC Census Generator
Week 2: Data Integrity & Statistical Validation

Generates a synthetic COSMIC fusion census dataset for fallback use when
real COSMIC data is unavailable. This ensures the pipeline remains operational
without requiring external COSMIC downloads.

REQUIREMENTS:
- Minimum 50 fusion pairs
- Includes real known fusions (BCR::ABL1, EML4::ALK, etc.)
- Recurrence counts follow power-law/heavy-tail distribution (Zipf or log-normal tail)

DESIGN PRINCIPLE:
Mock COSMIC exists ONLY to guarantee pipeline operability.
Real COSMIC always overrides mock when provided.
"""

import logging
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

_logger = logging.getLogger(__name__)

# Known real fusion pairs that must be included
KNOWN_FUSIONS: List[Tuple[str, str]] = [
    ("BCR", "ABL1"),
    ("EML4", "ALK"),
    ("TMPRSS2", "ERG"),
    ("EWSR1", "FLI1"),
    ("PML", "RARA"),
    ("KIF5B", "RET"),
    ("FGFR3", "TACC3"),
    ("SLC34A2", "ROS1"),
    ("TPM3", "NTRK1"),
]

# Additional synthetic fusion pairs to reach minimum 50
# These are realistic gene names but synthetic pairs
SYNTHETIC_FUSIONS: List[Tuple[str, str]] = [
    ("MYC", "IGH"),
    ("IGH", "BCL2"),
    ("IGH", "BCL6"),
    ("MYCN", "ALK"),
    ("RET", "NCOA4"),
    ("NTRK3", "ETV6"),
    ("ALK", "NPM1"),
    ("ROS1", "SLC34A2"),
    ("BRAF", "KIAA1549"),
    ("MET", "TPR"),
    ("FGFR1", "TACC1"),
    ("FGFR2", "TACC2"),
    ("NTRK2", "AFAP1"),
    ("NTRK3", "ETV6"),
    ("RET", "CCDC6"),
    ("ALK", "EML4"),
    ("ROS1", "CD74"),
    ("BRAF", "AKAP9"),
    ("MET", "CAPZA2"),
    ("FGFR1", "ZNF703"),
    ("FGFR2", "BICC1"),
    ("NTRK1", "LMNA"),
    ("NTRK2", "QKI"),
    ("NTRK3", "SQSTM1"),
    ("RET", "HOOK3"),
    ("ALK", "STRN"),
    ("ROS1", "EZR"),
    ("BRAF", "FAM131B"),
    ("MET", "KIF5B"),
    ("FGFR1", "OPHN1"),
    ("FGFR2", "OFD1"),
    ("NTRK1", "TPR"),
    ("NTRK2", "AFAP1"),
    ("NTRK3", "SQSTM1"),
    ("RET", "NCOA4"),
    ("ALK", "HIP1"),
    ("ROS1", "SLC34A2"),
    ("BRAF", "KIAA1549"),
    ("MET", "TPR"),
    ("FGFR1", "TACC1"),
    ("FGFR2", "TACC2"),
]

# Minimum number of fusion pairs required
MIN_FUSION_PAIRS = 50


def _generate_power_law_recurrence_counts(
    n: int, 
    seed: int = 42,
    min_count: int = 1,
    max_count: int = 10000,
    add_noise: bool = True,
) -> np.ndarray:
    """
    Generate recurrence counts following a power-law/heavy-tail distribution.
    
    Uses a Zipf-like distribution (power-law) where:
    - Most fusions have low recurrence
    - Few fusions have very high recurrence
    - Follows: P(k) ~ k^(-alpha) where alpha ~ 1.5-2.0
    
    With realism upgrades:
    - Adds sampling noise
    - Small random dropout
    - Recurrence jitter
    
    Args:
        n: Number of fusion pairs
        seed: Random seed for reproducibility
        min_count: Minimum recurrence count
        max_count: Maximum recurrence count
        add_noise: If True, add realistic noise and jitter
    
    Returns:
        Array of recurrence counts following power-law distribution
    """
    rng = np.random.default_rng(seed)
    
    # Generate Zipf-distributed ranks (alpha=1.5 gives realistic heavy tail)
    # Higher alpha = steeper tail (fewer high-recurrence fusions)
    alpha = 1.5
    ranks = rng.zipf(alpha, size=n)
    
    # Normalize ranks to desired range
    # Use log-space to preserve power-law structure
    min_rank = ranks.min()
    max_rank = ranks.max()
    
    # Map to log-space counts
    log_min = np.log(min_count)
    log_max = np.log(max_count)
    
    # Normalize ranks to [0, 1] then map to log-space
    normalized = (ranks - min_rank) / (max_rank - min_rank + 1e-10)
    log_counts = log_min + normalized * (log_max - log_min)
    
    # Convert back to linear space and round to integers
    counts = np.exp(log_counts)
    counts = np.round(counts).astype(int)
    counts = np.clip(counts, min_count, max_count)
    
    # Add realism: sampling noise, dropout, and jitter
    if add_noise:
        # Add small random jitter (5% coefficient of variation)
        noise_factor = 1.0 + rng.normal(0, 0.05, size=n)
        counts = (counts * noise_factor).astype(int)
        counts = np.clip(counts, min_count, max_count)
        
        # Small random dropout: 1-2% of counts become 0 (then reset to min_count)
        dropout_mask = rng.random(n) < 0.015
        counts[dropout_mask] = min_count
        
        # Recurrence jitter: add small random variations
        jitter = rng.integers(-2, 3, size=n)
        counts = counts + jitter
        counts = np.clip(counts, min_count, max_count)
    
    return counts


def generate_mock_cosmic_census(
    output_path: Path,
    min_pairs: int = MIN_FUSION_PAIRS,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Generate a mock COSMIC fusion census CSV file.
    
    Creates a synthetic dataset with:
    - Known real fusion pairs (BCR::ABL1, EML4::ALK, etc.)
    - Additional synthetic pairs to reach minimum count
    - Power-law distributed recurrence counts
    
    Args:
        output_path: Path where mock_cosmic_census.csv will be written
        min_pairs: Minimum number of fusion pairs to generate
        seed: Random seed for reproducibility
    
    Returns:
        DataFrame containing the generated mock COSMIC data
    
    Raises:
        ValueError: If min_pairs is less than number of known fusions
    """
    if min_pairs < len(KNOWN_FUSIONS):
        raise ValueError(
            f"min_pairs ({min_pairs}) must be >= {len(KNOWN_FUSIONS)} "
            f"(number of known fusions)"
        )
    
    # Collect all fusion pairs
    all_pairs: List[Tuple[str, str]] = []
    
    # Add known real fusions first (these get highest recurrence counts)
    all_pairs.extend(KNOWN_FUSIONS)
    
    # Add synthetic pairs to reach minimum
    n_needed = min_pairs - len(KNOWN_FUSIONS)
    synthetic_to_use = SYNTHETIC_FUSIONS[:n_needed]
    all_pairs.extend(synthetic_to_use)
    
    # If we still need more, generate additional synthetic pairs
    if len(all_pairs) < min_pairs:
        rng = np.random.default_rng(seed)
        # Generate additional synthetic gene names
        gene_prefixes = ["GENE", "FUS", "ONC", "TUM"]
        gene_suffixes = list(range(1, 1000))
        
        while len(all_pairs) < min_pairs:
            g1 = f"{rng.choice(gene_prefixes)}{rng.choice(gene_suffixes)}"
            g2 = f"{rng.choice(gene_prefixes)}{rng.choice(gene_suffixes)}"
            pair = (g1, g2)
            if pair not in all_pairs:
                all_pairs.append(pair)
    
    # Generate power-law recurrence counts with noise for realism
    n_pairs = len(all_pairs)
    recurrence_counts = _generate_power_law_recurrence_counts(
        n_pairs, 
        seed=seed,
        min_count=1,
        max_count=10000,
        add_noise=True,  # Enable noise for more realistic distribution
    )
    
    # Sort by recurrence count (descending) so known fusions get higher counts
    # This ensures known fusions appear in top ranks
    sorted_indices = np.argsort(recurrence_counts)[::-1]
    
    # Reorder pairs and counts
    sorted_pairs = [all_pairs[i] for i in sorted_indices]
    sorted_counts = recurrence_counts[sorted_indices]
    
    # Ensure known fusions get top recurrence counts
    known_indices = []
    other_indices = []
    for i, pair in enumerate(sorted_pairs):
        if pair in KNOWN_FUSIONS:
            known_indices.append(i)
        else:
            other_indices.append(i)
    
    # Assign highest counts to known fusions
    if known_indices:
        known_counts = sorted_counts[:len(known_indices)]
        other_counts = sorted_counts[len(known_indices):]
        
        # Rebuild with known fusions first
        final_pairs = []
        final_counts = []
        
        # Add known fusions with highest counts
        for idx in known_indices:
            final_pairs.append(sorted_pairs[idx])
        final_counts.extend(known_counts)
        
        # Add other fusions
        for idx in other_indices:
            final_pairs.append(sorted_pairs[idx])
        final_counts.extend(other_counts)
        
        sorted_pairs = final_pairs
        sorted_counts = np.array(final_counts)
    
    # Create DataFrame
    data = {
        "gene_1": [pair[0] for pair in sorted_pairs],
        "gene_2": [pair[1] for pair in sorted_pairs],
        "recurrence_count": sorted_counts.tolist(),
    }
    
    df = pd.DataFrame(data)
    
    # Write to CSV
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    
    _logger.info(
        f"Generated mock COSMIC census with {len(df)} fusion pairs "
        f"at {output_path}"
    )
    
    return df


def ensure_mock_cosmic_exists(output_dir: Path) -> Path:
    """
    Ensure mock COSMIC census file exists, creating it if necessary.
    
    Args:
        output_dir: Directory where mock_cosmic_census.csv should exist
    
    Returns:
        Path to the mock COSMIC census file
    """
    mock_path = output_dir / "mock_cosmic_census.csv"
    
    if not mock_path.exists():
        _logger.info(f"Mock COSMIC census not found at {mock_path}, generating...")
        generate_mock_cosmic_census(mock_path)
    
    return mock_path


if __name__ == "__main__":
    # Standalone execution: generate mock COSMIC census
    import sys
    
    if len(sys.argv) > 1:
        output_path = Path(sys.argv[1])
    else:
        # Default: generate in week2_validation/cosmic/
        script_dir = Path(__file__).parent
        output_path = script_dir / "mock_cosmic_census.csv"
    
    print(f"Generating mock COSMIC census at {output_path}...")
    df = generate_mock_cosmic_census(output_path)
    print(f"Generated {len(df)} fusion pairs")
    print(f"Top 10 fusions by recurrence:")
    print(df.head(10)[["gene_1", "gene_2", "recurrence_count"]].to_string(index=False))
