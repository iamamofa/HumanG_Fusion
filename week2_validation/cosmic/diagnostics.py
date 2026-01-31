"""
COSMIC Recurrence Diagnostic - Descriptive Rank-Order Comparison.

This module provides a DESCRIPTIVE ONLY diagnostic that compares
rank ordering between fusion recurrence data and COSMIC recurrence data.

WHAT DOES THIS MODULE DO?
This module compares the rank order of fusion pairs between two datasets:
1. Your fusion recurrence data
2. COSMIC fusion recurrence data

IT DOES NOT:
- Compute correlations (no Spearman, no Kendall)
- Compute p-values
- Compare distributions statistically
- Conclude agreement or disagreement
- Validate or invalidate data

IT ONLY:
- Compares rank ordering
- Reports counts and mismatches
- Prints descriptive summaries
"""


def run_cosmic_recurrence_diagnostic(
    *,
    fusion_df,
    cosmic_df,
    top_n: int = 10,
) -> dict:
    """
    Run a descriptive rank-order comparison between fusion and COSMIC data.
    
    This function performs a DESCRIPTIVE ONLY analysis:
    - NO statistical tests
    - NO correlations
    - NO inference
    - NO validation conclusions
    
    Args:
        fusion_df: DataFrame with columns: gene_1, gene_2, recurrence_count.
        cosmic_df: DataFrame with columns: gene_1, gene_2, recurrence_count.
        top_n: Number of top rank discrepancies to report (default: 10).
    
    Returns:
        A dictionary containing descriptive information only:
        - total_fusions_ours: Count of fusion pairs in our data
        - total_fusions_cosmic: Count of fusion pairs in COSMIC data
        - overlap_count: Count of pairs present in both datasets
        - only_in_ours_count: Count of pairs only in our data
        - only_in_cosmic_count: Count of pairs only in COSMIC data
        - top_rank_discrepancies: List of dicts with top rank mismatches
    """
    # -------------------------------------------------------------------------
    # STEP 1: Validate required columns
    # -------------------------------------------------------------------------
    required_columns = {"gene_1", "gene_2", "recurrence_count"}
    
    fusion_columns = set(fusion_df.columns) if fusion_df is not None else set()
    cosmic_columns = set(cosmic_df.columns) if cosmic_df is not None else set()
    
    # Check fusion_df columns
    if fusion_df is None:
        return {
            "total_fusions_ours": 0,
            "total_fusions_cosmic": 0,
            "overlap_count": 0,
            "only_in_ours_count": 0,
            "only_in_cosmic_count": 0,
            "top_rank_discrepancies": [],
            "message": "Fusion DataFrame is None.",
        }
    
    missing_fusion_cols = required_columns - fusion_columns
    if missing_fusion_cols:
        return {
            "total_fusions_ours": 0,
            "total_fusions_cosmic": 0,
            "overlap_count": 0,
            "only_in_ours_count": 0,
            "only_in_cosmic_count": 0,
            "top_rank_discrepancies": [],
            "message": f"Fusion data missing required columns: {missing_fusion_cols}",
        }
    
    # Check cosmic_df columns
    if cosmic_df is None:
        return {
            "total_fusions_ours": len(fusion_df),
            "total_fusions_cosmic": 0,
            "overlap_count": 0,
            "only_in_ours_count": len(fusion_df),
            "only_in_cosmic_count": 0,
            "top_rank_discrepancies": [],
            "message": "COSMIC DataFrame is None.",
        }
    
    missing_cosmic_cols = required_columns - cosmic_columns
    if missing_cosmic_cols:
        return {
            "total_fusions_ours": len(fusion_df),
            "total_fusions_cosmic": 0,
            "overlap_count": 0,
            "only_in_ours_count": len(fusion_df),
            "only_in_cosmic_count": 0,
            "top_rank_discrepancies": [],
            "message": f"COSMIC data missing required columns: {missing_cosmic_cols}",
        }
    
    # -------------------------------------------------------------------------
    # STEP 2: Normalize fusion pairs (uppercase, strip, sort alphabetically)
    # -------------------------------------------------------------------------
    def normalize_pair(gene_1, gene_2):
        """Normalize a gene pair: uppercase, strip, sort alphabetically."""
        g1 = str(gene_1).strip().upper()
        g2 = str(gene_2).strip().upper()
        # Return sorted tuple for consistent ordering
        return tuple(sorted([g1, g2]))
    
    # Build normalized lookup for fusion data
    fusion_pairs = {}
    for idx, row in fusion_df.iterrows():
        pair = normalize_pair(row["gene_1"], row["gene_2"])
        count = row["recurrence_count"]
        # If duplicate pairs, keep the one with higher count
        if pair not in fusion_pairs or count > fusion_pairs[pair]:
            fusion_pairs[pair] = count
    
    # Build normalized lookup for COSMIC data
    cosmic_pairs = {}
    for idx, row in cosmic_df.iterrows():
        pair = normalize_pair(row["gene_1"], row["gene_2"])
        count = row["recurrence_count"]
        # If duplicate pairs, keep the one with higher count
        if pair not in cosmic_pairs or count > cosmic_pairs[pair]:
            cosmic_pairs[pair] = count
    
    # -------------------------------------------------------------------------
    # STEP 3: Compute counts
    # -------------------------------------------------------------------------
    fusion_set = set(fusion_pairs.keys())
    cosmic_set = set(cosmic_pairs.keys())
    
    overlap = fusion_set & cosmic_set
    only_in_ours = fusion_set - cosmic_set
    only_in_cosmic = cosmic_set - fusion_set
    
    total_fusions_ours = len(fusion_set)
    total_fusions_cosmic = len(cosmic_set)
    overlap_count = len(overlap)
    only_in_ours_count = len(only_in_ours)
    only_in_cosmic_count = len(only_in_cosmic)
    
    # -------------------------------------------------------------------------
    # STEP 4: Compute ranks for overlapping pairs
    # -------------------------------------------------------------------------
    if overlap_count == 0:
        # No overlap - nothing to compare
        return {
            "total_fusions_ours": total_fusions_ours,
            "total_fusions_cosmic": total_fusions_cosmic,
            "overlap_count": overlap_count,
            "only_in_ours_count": only_in_ours_count,
            "only_in_cosmic_count": only_in_cosmic_count,
            "top_rank_discrepancies": [],
            "message": "No overlapping fusion pairs found between datasets.",
        }
    
    # Rank fusion pairs by recurrence_count (descending)
    # Rank 1 = highest recurrence count
    fusion_sorted = sorted(fusion_pairs.items(), key=lambda x: x[1], reverse=True)
    fusion_ranks = {pair: rank + 1 for rank, (pair, _) in enumerate(fusion_sorted)}
    
    cosmic_sorted = sorted(cosmic_pairs.items(), key=lambda x: x[1], reverse=True)
    cosmic_ranks = {pair: rank + 1 for rank, (pair, _) in enumerate(cosmic_sorted)}
    
    # -------------------------------------------------------------------------
    # STEP 5: Compute rank discrepancies for overlapping pairs
    # -------------------------------------------------------------------------
    discrepancies = []
    for pair in overlap:
        our_rank = fusion_ranks[pair]
        cosmic_rank = cosmic_ranks[pair]
        abs_diff = abs(our_rank - cosmic_rank)
        discrepancies.append({
            "gene_1": pair[0],
            "gene_2": pair[1],
            "our_rank": our_rank,
            "cosmic_rank": cosmic_rank,
            "absolute_rank_difference": abs_diff,
            "our_recurrence_count": fusion_pairs[pair],
            "cosmic_recurrence_count": cosmic_pairs[pair],
        })
    
    # Sort by absolute rank difference (descending) to get top mismatches
    discrepancies.sort(key=lambda x: x["absolute_rank_difference"], reverse=True)
    
    # Take top N
    top_discrepancies = discrepancies[:top_n]
    
    # -------------------------------------------------------------------------
    # STEP 6: Return descriptive dictionary
    # -------------------------------------------------------------------------
    return {
        "total_fusions_ours": total_fusions_ours,
        "total_fusions_cosmic": total_fusions_cosmic,
        "overlap_count": overlap_count,
        "only_in_ours_count": only_in_ours_count,
        "only_in_cosmic_count": only_in_cosmic_count,
        "top_rank_discrepancies": top_discrepancies,
    }
