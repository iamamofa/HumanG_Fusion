"""
Week 2: Data Integrity & Statistical Validation — Bootstrap Visualization.

Generates bootstrap correlation distribution plots.
Uses ONLY matplotlib (no color specification).
"""

import sys
from pathlib import Path
from typing import List, Optional

try:
    import matplotlib
    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def generate_bootstrap_rho_plot(
    bootstrap_rho_values: List[float],
    output_path: Path,
) -> Optional[Path]:
    """
    Generate bootstrap correlation distribution plot.
    
    Args:
        bootstrap_rho_values: List of bootstrap rho values.
        output_path: Path where PNG plot should be saved.
    
    Returns:
        Path to generated plot file, or None if generation failed.
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Warning: Matplotlib not available. Bootstrap plot generation skipped.", file=sys.stderr)
        return None
    
    if not bootstrap_rho_values or len(bootstrap_rho_values) < 10:
        print("Warning: Insufficient bootstrap samples for plot generation.", file=sys.stderr)
        return None
    
    try:
        # Create figure with white background
        fig, ax = plt.subplots(figsize=(8, 5), facecolor="white")
        ax.set_facecolor("white")
        
        # Create histogram
        n, bins, patches = ax.hist(
            bootstrap_rho_values,
            bins=30,
            alpha=0.8,
            edgecolor="black",
            facecolor="#4a90e2",
            linewidth=0.5,
        )
        
        # Add vertical line for mean
        mean_rho = np.mean(bootstrap_rho_values)
        ax.axvline(mean_rho, color="black", linestyle="--", linewidth=1.5, label=f"Mean ρ = {mean_rho:.3f}")
        
        # Add CI lines
        ci_lower = np.percentile(bootstrap_rho_values, 2.5)
        ci_upper = np.percentile(bootstrap_rho_values, 97.5)
        ax.axvline(ci_lower, color="black", linestyle=":", linewidth=1, label=f"95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
        ax.axvline(ci_upper, color="black", linestyle=":", linewidth=1)
        
        # Set labels and title
        ax.set_xlabel("Bootstrap Spearman ρ", fontsize=11, color="black")
        ax.set_ylabel("Frequency", fontsize=11, color="black")
        ax.set_title("Bootstrap Correlation Distribution", fontsize=12, fontweight="bold", color="black")
        
        # Style axes
        ax.spines["top"].set_color("black")
        ax.spines["bottom"].set_color("black")
        ax.spines["left"].set_color("black")
        ax.spines["right"].set_color("black")
        ax.tick_params(colors="black")
        
        # Grid for readability
        ax.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)
        
        # Legend
        ax.legend(fontsize=9, framealpha=0.9)
        
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
        print(f"Warning: Bootstrap plot generation failed: {e}", file=sys.stderr)
        return None
