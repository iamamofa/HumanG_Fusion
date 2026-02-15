"""
Week 2: Data Integrity & Statistical Validation — COSMIC Bias Quantification.

This module quantifies sampling bias in COSMIC data by computing:
- Gene concentration index (top 5 gene dominance)
- Long tail coverage ratio
- Gini coefficient for gene recurrence inequality

SAFETY: This module does NOT modify COSMIC data or core formulas.
It only analyzes and reports bias metrics.
"""

import sys
from typing import Dict, Optional

import numpy as np
import pandas as pd


def compute_gini_coefficient(values: np.ndarray) -> float:
    """
    Compute Gini coefficient for inequality measurement.
    
    Gini coefficient ranges from 0 (perfect equality) to 1 (perfect inequality).
    
    Args:
        values: Array of non-negative values
    
    Returns:
        Gini coefficient (0.0 to 1.0)
    """
    if len(values) == 0:
        return 0.0
    
    # Remove zeros and sort
    values = values[values > 0]
    if len(values) == 0:
        return 0.0
    
    values = np.sort(values)
    n = len(values)
    
    # Compute Gini coefficient using formula:
    # G = (2 * sum(i * x_i)) / (n * sum(x_i)) - (n + 1) / n
    # where x_i are sorted values
    cumsum = np.cumsum(values)
    numerator = np.sum((np.arange(1, n + 1)) * values)
    denominator = n * np.sum(values)
    
    if denominator == 0:
        return 0.0
    
    gini = (2 * numerator) / denominator - (n + 1) / n
    
    return float(np.clip(gini, 0.0, 1.0))


def quantify_cosmic_sampling_bias(cosmic_df: pd.DataFrame) -> Dict[str, Optional[float]]:
    """
    Quantify COSMIC sampling bias by computing concentration and inequality metrics.
    
    METRICS:
    - top5_gene_fraction: Fraction of total recurrence accounted for by top 5 genes
    - gini_gene_recurrence: Gini coefficient for gene recurrence inequality
    - tail_coverage_ratio: Ratio of unique genes in long tail vs top genes
    - bias_severity_classification: Classification of bias severity
    
    Args:
        cosmic_df: DataFrame with columns: gene_1, gene_2, recurrence_count
    
    Returns:
        Dictionary with:
        - top5_gene_fraction: Fraction of recurrence from top 5 genes (0.0 to 1.0)
        - gini_gene_recurrence: Gini coefficient for recurrence inequality
        - tail_coverage_ratio: Ratio of tail genes to top genes
        - bias_severity_classification: "LOW", "MODERATE", "HIGH", or "N/A"
        - n_unique_genes: Number of unique genes in dataset
        - n_fusion_pairs: Number of unique fusion pairs
    """
    base_result = {
        "top5_gene_fraction": None,
        "gini_gene_recurrence": None,
        "tail_coverage_ratio": None,
        "bias_severity_classification": "N/A",
        "n_unique_genes": 0,
        "n_fusion_pairs": 0,
    }
    
    try:
        # Build gene recurrence dictionary (count total recurrence per gene)
        gene_recurrence = {}
        
        if "gene_1" in cosmic_df.columns and "gene_2" in cosmic_df.columns:
            for _, row in cosmic_df.iterrows():
                gene_1 = str(row["gene_1"]).strip().upper()
                gene_2 = str(row["gene_2"]).strip().upper()
                recurrence = float(row.get("recurrence_count", 0))
                
                # Count recurrence for each gene
                gene_recurrence[gene_1] = gene_recurrence.get(gene_1, 0) + recurrence
                gene_recurrence[gene_2] = gene_recurrence.get(gene_2, 0) + recurrence
        
        if len(gene_recurrence) == 0:
            return base_result
        
        # Count unique fusion pairs
        fusion_pairs = set()
        if "gene_1" in cosmic_df.columns and "gene_2" in cosmic_df.columns:
            for _, row in cosmic_df.iterrows():
                gene_1 = str(row["gene_1"]).strip().upper()
                gene_2 = str(row["gene_2"]).strip().upper()
                pair = tuple(sorted([gene_1, gene_2]))
                fusion_pairs.add(pair)
        
        base_result["n_unique_genes"] = len(gene_recurrence)
        base_result["n_fusion_pairs"] = len(fusion_pairs)
        
        # Sort genes by recurrence (descending)
        gene_recurrence_sorted = sorted(gene_recurrence.items(), key=lambda x: x[1], reverse=True)
        recurrence_values = np.array([count for _, count in gene_recurrence_sorted])
        
        total_recurrence = np.sum(recurrence_values)
        if total_recurrence == 0:
            return base_result
        
        # Compute top 5 gene fraction
        n_top = min(5, len(gene_recurrence_sorted))
        top5_recurrence = np.sum(recurrence_values[:n_top])
        top5_fraction = float(top5_recurrence / total_recurrence)
        base_result["top5_gene_fraction"] = top5_fraction
        
        # Compute Gini coefficient for gene recurrence inequality
        gini = compute_gini_coefficient(recurrence_values)
        base_result["gini_gene_recurrence"] = gini
        
        # Compute tail coverage ratio
        # Define "top" as top 10% of genes, "tail" as bottom 50% of genes
        n_genes = len(gene_recurrence_sorted)
        if n_genes >= 10:
            n_top_genes = max(1, int(n_genes * 0.10))
            n_tail_genes = max(1, int(n_genes * 0.50))
            
            # Top genes: highest recurrence
            top_genes = set(gene for gene, _ in gene_recurrence_sorted[:n_top_genes])
            
            # Tail genes: lowest recurrence (but exclude zeros)
            tail_start_idx = max(0, n_genes - n_tail_genes)
            tail_genes = set(gene for gene, _ in gene_recurrence_sorted[tail_start_idx:] if gene_recurrence[gene] > 0)
            
            if len(top_genes) > 0:
                tail_coverage_ratio = float(len(tail_genes) / len(top_genes))
                base_result["tail_coverage_ratio"] = tail_coverage_ratio
            else:
                base_result["tail_coverage_ratio"] = None
        else:
            base_result["tail_coverage_ratio"] = None
        
        # Classify bias severity
        # LOW: top5_fraction < 0.3, gini < 0.5
        # MODERATE: top5_fraction 0.3-0.6, gini 0.5-0.7
        # HIGH: top5_fraction > 0.6, gini > 0.7
        if top5_fraction < 0.3 and gini < 0.5:
            classification = "LOW"
        elif top5_fraction < 0.6 and gini < 0.7:
            classification = "MODERATE"
        elif top5_fraction >= 0.6 or gini >= 0.7:
            classification = "HIGH"
        else:
            classification = "MODERATE"  # Default
        
        base_result["bias_severity_classification"] = classification
        
        return base_result
        
    except Exception as e:
        print(f"Warning: COSMIC bias quantification failed: {e}", file=sys.stderr)
        return base_result
