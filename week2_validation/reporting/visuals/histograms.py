"""
Week 2: Data Integrity & Statistical Validation — Histogram Generation.

Generates scientific publication-style histograms for protein length distributions.
This module is optional and does not affect pipeline execution if unavailable.
"""

import sys
from pathlib import Path
from typing import Optional

try:
    import matplotlib
    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def generate_protein_length_histogram(
    df: pd.DataFrame,
    output_path: Path,
    column: str = "protein_length",
) -> Optional[Path]:
    """
    Generate a scientific publication-style histogram of protein length distribution.
    
    Args:
        df: DataFrame containing protein length data.
        output_path: Path where PNG histogram should be saved.
        column: Name of column containing protein length values (default: "protein_length").
    
    Returns:
        Path to generated histogram file, or None if generation failed.
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Warning: Matplotlib not available. Histogram generation skipped.", file=sys.stderr)
        return None
    
    if column not in df.columns:
        print(f"Warning: Column '{column}' not found in DataFrame. Histogram generation skipped.", file=sys.stderr)
        return None
    
    try:
        # Extract valid protein length values
        values = df[column].dropna()
        values = values[values > 0]  # Only positive values
        
        if len(values) == 0:
            print(f"Warning: No valid protein length values found. Histogram generation skipped.", file=sys.stderr)
            return None
        
        # Create figure with white background
        fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
        ax.set_facecolor("white")
        
        # Create histogram
        n, bins, patches = ax.hist(
            values,
            bins=50,
            alpha=0.8,
            edgecolor="black",
            facecolor="#4a90e2",
            linewidth=0.5,
        )
        
        # Set labels and title
        ax.set_xlabel("Protein Length (aa)", fontsize=12, color="black")
        ax.set_ylabel("Frequency", fontsize=12, color="black")
        ax.set_title("Protein Length Distribution", fontsize=14, fontweight="bold", color="black")
        
        # Style axes
        ax.spines["top"].set_color("black")
        ax.spines["bottom"].set_color("black")
        ax.spines["left"].set_color("black")
        ax.spines["right"].set_color("black")
        ax.tick_params(colors="black")
        
        # Grid for readability
        ax.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)
        
        # Tight layout
        plt.tight_layout()
        
        # Save at 300 DPI
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            output_path,
            dpi=300,
            format="png",
            facecolor="white",
            edgecolor="none",
            bbox_inches="tight",
        )
        
        plt.close(fig)
        
        return output_path
        
    except Exception as e:
        print(f"Warning: Histogram generation failed: {e}", file=sys.stderr)
        return None


def generate_log_protein_length_histogram(
    df: pd.DataFrame,
    output_path: Path,
    column: str = "protein_length",
) -> Optional[Path]:
    """
    Generate a log10-transformed histogram of protein length distribution.
    
    Args:
        df: DataFrame containing protein length data.
        output_path: Path where PNG histogram should be saved.
        column: Name of column containing protein length values (default: "protein_length").
    
    Returns:
        Path to generated histogram file, or None if generation failed.
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Warning: Matplotlib not available. Log histogram generation skipped.", file=sys.stderr)
        return None
    
    if column not in df.columns:
        print(f"Warning: Column '{column}' not found in DataFrame. Log histogram generation skipped.", file=sys.stderr)
        return None
    
    try:
        # Extract valid positive protein length values
        values = df[column].dropna()
        values = values[values > 0]  # Only positive values for log transform
        
        if len(values) == 0:
            print(f"Warning: No valid positive protein length values found. Log histogram generation skipped.", file=sys.stderr)
            return None
        
        # Apply log10 transform
        log_values = np.log10(values)
        
        # Create figure with white background
        fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
        ax.set_facecolor("white")
        
        # Create histogram
        n, bins, patches = ax.hist(
            log_values,
            bins=50,
            alpha=0.8,
            edgecolor="black",
            facecolor="#4a90e2",
            linewidth=0.5,
        )
        
        # Set labels and title
        ax.set_xlabel("Log10(Protein Length)", fontsize=12, color="black")
        ax.set_ylabel("Frequency", fontsize=12, color="black")
        ax.set_title("Log10 Protein Length Distribution", fontsize=14, fontweight="bold", color="black")
        
        # Style axes
        ax.spines["top"].set_color("black")
        ax.spines["bottom"].set_color("black")
        ax.spines["left"].set_color("black")
        ax.spines["right"].set_color("black")
        ax.tick_params(colors="black")
        
        # Grid for readability
        ax.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)
        
        # Tight layout
        plt.tight_layout()
        
        # Save at 300 DPI
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            output_path,
            dpi=300,
            format="png",
            facecolor="white",
            edgecolor="none",
            bbox_inches="tight",
        )
        
        plt.close(fig)
        
        return output_path
        
    except Exception as e:
        print(f"Warning: Log histogram generation failed: {e}", file=sys.stderr)
        return None
