"""
Week 1: Pipeline Execution & Data Generation — Compatibility Adapter.

Converts Week 1 (Pipeline Execution & Data Generation) output format to
Week 2: Data Integrity & Statistical Validation schema.
Never modifies original DataFrame; always returns a new copy.
"""

import pandas as pd

# Clamp recurrence_frequency before * 1000 to avoid float/int overflow
# int32 max ~2.1e9; 1000 * freq <= 2e9 => freq <= 2e6
SAFE_FREQ_LIMIT = 2.0e6


def adapt_week1_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adapt DataFrame to Week 2: Data Integrity & Statistical Validation schema.
    Accepts Week 1 (Pipeline Execution & Data Generation) style columns.

    If input already matches the schema (fusion_id, gene_1, gene_2,
    protein_length, recurrence_count), returns unchanged copy.

    Otherwise maps Week 1 (Pipeline Execution & Data Generation) columns:
        geneA -> gene_1
        geneB -> gene_2
        recurrence_frequency OR samples_detected -> recurrence_count
        fusion_id = gene_1 + "::" + gene_2 if missing
        protein_length = NaN if missing (diagnostics that require it may skip)

    Safety: Never modifies original DataFrame. Always returns new copy.

    Args:
        df: Input DataFrame (Week 1: Pipeline Execution & Data Generation or Week 2: Data Integrity & Statistical Validation format).

    Returns:
        New DataFrame conforming to Week 2: Data Integrity & Statistical Validation schema.

    Raises:
        ValueError: If adaptation is not possible.
    """
    if df is None or not isinstance(df, pd.DataFrame):
        raise ValueError("Input must be a non-null pandas DataFrame")

    out = df.copy()
    cols = set(out.columns)
    adaptation_applied = []

    # Week 2: Data Integrity & Statistical Validation required columns
    w2_fusion_id = "fusion_id" in cols
    w2_gene_1 = "gene_1" in cols
    w2_gene_2 = "gene_2" in cols
    w2_protein_length = "protein_length" in cols
    w2_recurrence_count = "recurrence_count" in cols

    if w2_fusion_id and w2_gene_1 and w2_gene_2 and w2_recurrence_count and w2_protein_length:
        return out

    # Map geneA -> gene_1
    if not w2_gene_1 and "geneA" in cols:
        out["gene_1"] = out["geneA"].astype(str)
        adaptation_applied.append("geneA -> gene_1")
    elif not w2_gene_1:
        raise ValueError("Cannot adapt: missing gene_1 and geneA")

    # Map geneB -> gene_2
    if not w2_gene_2 and "geneB" in cols:
        out["gene_2"] = out["geneB"].astype(str)
        adaptation_applied.append("geneB -> gene_2")
    elif not w2_gene_2:
        raise ValueError("Cannot adapt: missing gene_2 and geneB")

    # Map recurrence
    if not w2_recurrence_count:
        if "samples_detected" in cols:
            out["recurrence_count"] = pd.to_numeric(out["samples_detected"], errors="coerce").fillna(0).astype(int)
            adaptation_applied.append("samples_detected -> recurrence_count")
        elif "recurrence_frequency" in cols:
            freq = pd.to_numeric(out["recurrence_frequency"], errors="coerce").fillna(0)
            freq = freq.clip(upper=SAFE_FREQ_LIMIT)
            scaled = (freq * 1000).round(0).clip(lower=0)
            out["recurrence_count"] = scaled.astype("int64")
            adaptation_applied.append("recurrence_frequency -> recurrence_count (scaled)")
        else:
            raise ValueError("Cannot adapt: missing recurrence_count, samples_detected, recurrence_frequency")

    # Generate fusion_id if missing
    if not w2_fusion_id:
        g1 = out["gene_1"].astype(str)
        g2 = out["gene_2"].astype(str)
        out["fusion_id"] = g1 + "::" + g2
        adaptation_applied.append("fusion_id = gene_1::gene_2")

    # protein_length: if missing, set NaN
    if not w2_protein_length:
        out["protein_length"] = float("nan")
        adaptation_applied.append("protein_length = NaN (not available)")

    if adaptation_applied:
        print("  Week 1 (Pipeline Execution & Data Generation) adapter applied:", ", ".join(adaptation_applied))

    return out
