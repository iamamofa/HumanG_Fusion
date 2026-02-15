"""
Week 2: Data Integrity & Statistical Validation — Skewness Visualization.

Generates scientific publication-style skewness visualization (histogram + boxplot).
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


def generate_skewness_visualization(
    df: pd.DataFrame,
    output_path: Path,
    column: str = "protein_length",
) -> Optional[Path]:
    """
    Generate skewness visualization with histogram and boxplot side by side.
    
    Args:
        df: DataFrame containing protein length data.
        output_path: Path where PNG visualization should be saved.
        column: Name of column containing protein length values (default: "protein_length").
    
    Returns:
        Path to generated visualization file, or None if generation failed.
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Warning: Matplotlib not available. Skewness visualization skipped.", file=sys.stderr)
        return None
    
    if column not in df.columns:
        print(f"Warning: Column '{column}' not found in DataFrame. Skewness visualization skipped.", file=sys.stderr)
        return None
    
    try:
        # Extract valid protein length values
        values = df[column].dropna()
        values = values[values > 0]  # Only positive values
        
        if len(values) == 0:
            print(f"Warning: No valid protein length values found. Skewness visualization skipped.", file=sys.stderr)
            return None
        
        # Create figure with two subplots side by side
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), facecolor="white")
        
        # Left panel: Histogram
        ax1.set_facecolor("white")
        n, bins, patches = ax1.hist(
            values,
            bins=50,
            alpha=0.8,
            edgecolor="black",
            facecolor="#4a90e2",
            linewidth=0.5,
        )
        ax1.set_xlabel("Protein Length (aa)", fontsize=11, color="black")
        ax1.set_ylabel("Frequency", fontsize=11, color="black")
        ax1.set_title("Histogram", fontsize=12, fontweight="bold", color="black")
        ax1.spines["top"].set_color("black")
        ax1.spines["bottom"].set_color("black")
        ax1.spines["left"].set_color("black")
        ax1.spines["right"].set_color("black")
        ax1.tick_params(colors="black")
        ax1.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)
        
        # Right panel: Boxplot
        ax2.set_facecolor("white")
        bp = ax2.boxplot(
            [values],
            vert=True,
            patch_artist=True,
            widths=0.6,
            showmeans=True,
            meanline=True,
        )
        # Style boxplot
        for patch in bp["boxes"]:
            patch.set_facecolor("#4a90e2")
            patch.set_edgecolor("black")
            patch.set_alpha(0.8)
        for element in ["whiskers", "fliers", "medians", "caps"]:
            for item in bp[element]:
                item.set_color("black")
                item.set_linewidth(1)
        ax2.set_ylabel("Protein Length (aa)", fontsize=11, color="black")
        ax2.set_title("Boxplot", fontsize=12, fontweight="bold", color="black")
        ax2.spines["top"].set_color("black")
        ax2.spines["bottom"].set_color("black")
        ax2.spines["left"].set_color("black")
        ax2.spines["right"].set_color("black")
        ax2.tick_params(colors="black")
        ax2.grid(True, alpha=0.3, linestyle="--", linewidth=0.5, axis="y")
        ax2.set_xticklabels(["Protein Length"])
        
        # Overall title
        fig.suptitle(
            "Distribution Shape Visualization (Skewness Assessment)",
            fontsize=14,
            fontweight="bold",
            color="black",
            y=1.02,
        )
        
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
        print(f"Warning: Skewness visualization generation failed: {e}", file=sys.stderr)
        return None
