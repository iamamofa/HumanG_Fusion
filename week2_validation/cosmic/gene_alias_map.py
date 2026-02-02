"""
HGNC-Style Gene Alias Normalization Map
Week 2: Data Integrity & Statistical Validation

This module provides gene alias mapping for canonical gene name normalization.
Common gene aliases are mapped to their canonical HGNC-approved symbols to ensure
consistent fusion pair matching across datasets.

WHAT DOES THIS MODULE DO?
This module maps common gene aliases to their canonical names, ensuring that
fusion pairs like "P53::ABL1" and "TP53::ABL1" are recognized as the same pair.

DESIGN PRINCIPLE:
- Maps aliases to canonical symbols (not bidirectional)
- Applied BEFORE pair canonicalization
- Preserves biological accuracy while improving matching
"""

from typing import Dict

# HGNC-style alias mapping: alias -> canonical_symbol
# This dictionary maps common gene aliases to their canonical HGNC-approved symbols
GENE_ALIAS_MAP: Dict[str, str] = {
    # Tumor suppressor genes
    "P53": "TP53",
    "P73": "TP73",
    "P63": "TP63",
    
    # Histone methyltransferases
    "MLL": "KMT2A",
    "MLL1": "KMT2A",
    "MLL2": "KMT2D",
    "MLL3": "KMT2C",
    "MLL4": "KMT2B",
    
    # Receptor tyrosine kinases
    "HER2": "ERBB2",
    "HER1": "EGFR",
    "HER3": "ERBB3",
    "HER4": "ERBB4",
    
    # Transcription factors
    "EWS": "EWSR1",
    "EWS1": "EWSR1",
    
    # Additional common aliases
    "ALK1": "ACVRL1",  # Activin A receptor like type 1
    "ALK2": "ACVR1",   # Activin A receptor type 1
    "ALK5": "TGFBR1",  # Transforming growth factor beta receptor 1
    
    # Note: ALK (anaplastic lymphoma kinase) is canonical, no alias needed
}

def normalize_gene_alias(gene_name: str) -> str:
    """
    Normalize a gene name by mapping aliases to canonical symbols.
    
    If the gene name is an alias, returns the canonical symbol.
    If the gene name is already canonical or not in the map, returns it unchanged.
    
    Args:
        gene_name: Gene name (may be alias or canonical)
    
    Returns:
        Canonical gene symbol (or original if not in alias map)
    """
    if not gene_name:
        return gene_name
    
    # Convert to uppercase for case-insensitive matching
    gene_upper = str(gene_name).strip().upper()
    
    # Check if it's an alias
    canonical = GENE_ALIAS_MAP.get(gene_upper)
    
    if canonical:
        return canonical
    
    # Not an alias, return original (normalized to uppercase)
    return gene_upper
