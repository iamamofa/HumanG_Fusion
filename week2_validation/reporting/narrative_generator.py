"""
Week 2: Data Integrity & Statistical Validation — Narrative Report Generator.

Converts week2_status.json, week2_dataset_certification.json, and
week2_diagnostic_results.json into a researcher-friendly Markdown report.

The JSON files are for machines; this report is for humans.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional


STATUS_FILENAME_PATTERN = "week2_status_{stem}.json"
CERTIFICATION_FILENAME_PATTERN = "week2_dataset_certification_{stem}.json"
DIAGNOSTIC_RESULTS_FILENAME_PATTERN = "week2_diagnostic_results_{stem}.json"
REPORT_FILENAME_PATTERN = "Statistical_Integrity_Report_{stem}.md"
FAILURE_REPORT_FILENAME_PATTERN = "Failure_Report_{stem}.md"


def _load_json(path: Path) -> Optional[Dict[str, Any]]:
    """Load JSON file; return None if missing or invalid."""
    if not path.is_file():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


def _interpret_skewness(skewness: float) -> str:
    """Provide a brief interpretive sentence for skewness."""
    if abs(skewness) < 0.5:
        return "indicates a roughly symmetrical distribution."
    if skewness > 0:
        return "indicates right-skewed (positive skew) distribution."
    return "indicates left-skewed (negative skew) distribution."


def _section_1_executive_summary(
    status: Dict[str, Any],
    certification: Optional[Dict[str, Any]],
) -> str:
    """Build Section 1: Executive Summary."""
    dataset_hash = status.get("dataset_hash", "N/A")
    dq = status.get("data_quality") or {}
    total_records = dq.get("total_rows_original", dq.get("rows_used_for_analysis", "N/A"))
    if isinstance(total_records, int):
        total_records = str(total_records)

    approved = False
    if certification:
        approved = certification.get("approved_for_powerlaw_modeling", False)

    status_word = "**PASS**" if approved else "Not approved"
    return f"""## 1. Executive Summary

| Field | Value |
|-------|-------|
| **Dataset Hash** | `{dataset_hash}` |
| **Total Records** | {total_records} |
| **Approval for Week 3 (Power-Law Modeling)** | {status_word} |

The dataset has completed the Data Integrity and Statistical Validation layer. 
{'The dataset is approved for downstream power-law modeling (Week 3).' if approved else 'The dataset has not been approved for downstream modeling.'}
"""


def _section_2_distribution(
    status: Dict[str, Any],
    diag: Optional[Dict[str, Any]],
    output_dir: Optional[Path] = None,
) -> str:
    """Build Section 2: Distribution Analysis."""
    dq = status.get("data_quality") or {}
    dist = (diag or {}).get("distribution") or {}

    min_val = dist.get("minimum", dq.get("min"))
    max_val = dist.get("maximum", dq.get("max"))
    mean_val = dist.get("mean", dq.get("mean"))
    median_val = dist.get("median", dq.get("median"))
    skewness = dq.get("skewness", dist.get("skewness"))

    def fmt(v):
        if v is None:
            return "N/A"
        try:
            return f"{float(v):.2f}"
        except (TypeError, ValueError):
            return str(v)

    skew_str = f"{float(skewness):.4f}" if skewness is not None and isinstance(skewness, (int, float)) else fmt(skewness)

    # Determine range label for narrative
    if skewness is not None and isinstance(skewness, (int, float)):
        sk = float(skewness)
        if abs(sk) < 0.5:
            range_label = "Symmetrical"
        elif sk > 0:
            range_label = "Positive Skew"
        else:
            range_label = "Negative Skew"
        interp = f"\n\nThe calculated skewness of {skew_str} is visualized in the Live Skewness Indicator above, placing the dataset in the **{range_label}** range."
    else:
        interp = ""

    # Build section: table, then actual histogram (if exists), then skewness diagnostic (if exists)
    fig_num = 1
    img_blocks = []
    if output_dir:
        if (output_dir / "protein_distribution.png").is_file():
            img_blocks.append(
                "![Protein Length Distribution](protein_distribution.png)\n\n"
                "*Figure 1: Actual protein length distribution (histogram) from the dataset.*"
            )
            fig_num = 2
        if (output_dir / "skewness_diagnostic.png").is_file():
            img_blocks.append(
                "![Live Skewness Indicator](skewness_diagnostic.png)\n\n"
                f"*Figure {fig_num}: Skewness gauge and theoretical density curve matching the dataset skewness.*"
            )
    else:
        img_blocks.append(
            "![Protein Length Distribution](protein_distribution.png)\n\n"
            "*Figure 1: Actual protein length distribution (histogram) from the dataset.*\n\n"
            "![Live Skewness Indicator](skewness_diagnostic.png)\n\n"
            "*Figure 2: Skewness gauge and theoretical density curve matching the dataset skewness.*"
        )
    img_section = "\n\n".join(img_blocks) if img_blocks else ""

    return f"""## 2. Distribution Analysis

| Statistic | Value |
|-----------|-------|
| Min | {fmt(min_val)} |
| Max | {fmt(max_val)} |
| Mean | {fmt(mean_val)} |
| Median | {fmt(median_val)} |
| Skewness | {skew_str} |

{img_section}
{interp}
"""


def _section_3_statistical_diagnostics(diag: Optional[Dict[str, Any]]) -> str:
    """Build Section 3: Statistical Diagnostics (KS, Anderson-Darling)."""
    ln = (diag or {}).get("log_normality") or {}
    ks = ln.get("ks_statistic")
    ad = ln.get("ad_statistic")

    def fmt(v):
        if v is None:
            return "N/A"
        try:
            return f"{float(v):.6f}"
        except (TypeError, ValueError):
            return str(v)

    return f"""## 3. Statistical Diagnostics

| Test | Value |
|------|-------|
| **Kolmogorov-Smirnov statistic** | {fmt(ks)} |
| **Anderson-Darling statistic** | {fmt(ad)} |

*These are descriptive statistics only. No inference or pass/fail determination is made.*
"""


def _section_4_benford(diag: Optional[Dict[str, Any]], notes: list) -> str:
    """Build Section 4: Benford's Law."""
    ben = (diag or {}).get("benford") or {}
    observed = ben.get("observed_frequencies") or {}
    expected = ben.get("expected_frequencies") or {}
    scale_span = ben.get("scale_span_orders_of_magnitude")
    applicability = ben.get("applicability")
    reason = ben.get("reason_if_not_applicable", "")

    scale_warning = ""
    for n in notes or []:
        if "scale" in n.lower() or "benford" in n.lower() or "power-law" in n.lower():
            scale_warning = f"\n\n**{n}**"
            break
    if scale_span is not None and float(scale_span) < 2.0 and not scale_warning:
        scale_warning = (
            f"\n\n**WARNING:** Scale span ({float(scale_span):.2f} orders of magnitude) "
            "is below the minimum threshold (2.0) typically associated with "
            "Benford-distributed data. Statistical range insufficient for reliable "
            "Benford and Power-Law modeling."
        )

    rows = []
    for d in range(1, 10):
        obs = observed.get(d, observed.get(str(d), 0))
        exp = expected.get(d, expected.get(str(d), 0)) if expected else 0
        if isinstance(obs, (int, float)):
            obs_f = f"{float(obs):.4f}"
        else:
            obs_f = str(obs)
        if isinstance(exp, (int, float)):
            exp_f = f"{float(exp):.4f}"
        else:
            exp_f = str(exp)
        rows.append(f"| {d} | {obs_f} | {exp_f} |")

    table_body = "\n".join(rows) if rows else "| *No data* | - | - |"

    appl_note = ""
    if applicability is not None:
        appl_note = f"\n\n**Applicability (heuristic):** {'Applicable' if applicability else 'Not applicable'}"
    if reason:
        appl_note += f"\n\n**Reason:** {reason}"

    return f"""## 4. Benford's Law — First Significant Digit Analysis

| Digit | Observed | Benford Expected |
|-------|----------|------------------|
{table_body}
{appl_note}
{scale_warning}

*Diagnostic only. No inference drawn.*
"""


def generate_failure_report(
    output_dir: Path,
    dataset_stem: str,
    failure_reason: str,
    failure_details: str = "",
) -> Optional[Path]:
    """
    Generate a Markdown Failure Report when validation fails.

    Provides a paper trail for every file that fails validation.

    Args:
        output_dir: Directory for the failure report.
        dataset_stem: Stem from input filename.
        failure_reason: Human-readable reason (e.g. "Missing Columns", "Empty File").
        failure_details: Additional context or error message.

    Returns:
        Path to the written report, or None on write failure.
    """
    output_dir = Path(output_dir).resolve()
    if not output_dir.is_dir():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            return None

    path = output_dir / FAILURE_REPORT_FILENAME_PATTERN.format(stem=dataset_stem)
    content = f"""# Data Integrity Validation — Failure Report

*Dataset: {dataset_stem}*

---

## Validation Status: **FAILED**

## Failure Reason

**{failure_reason}**

"""
    if failure_details:
        content += f"""
## Details

```
{failure_details}
```

"""

    content += """
---

*This report was generated automatically when the dataset failed validation.*
*No statistical analysis or certification was performed.*
"""
    try:
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        return path
    except OSError:
        return None


def _build_method_justification_section() -> str:
    """Build the statistical method justification section."""
    return """

### Statistical Method Justification

The statistical methods used in this COSMIC cross-validation are specifically chosen for 
their appropriateness to biological fusion recurrence data:

**Spearman Rank Correlation**: Non-parametric method robust to the heavy-tailed, non-normal 
distributions characteristic of recurrence count data. Measures rank-order agreement rather 
than exact value correspondence, which is appropriate for comparing recurrence rankings 
across different cohorts.

**Bootstrap Confidence Intervals**: Provide uncertainty quantification without parametric 
assumptions. Particularly important when overlap sample sizes are limited.

**Hypergeometric Enrichment Test**: Tests whether top-ranked fusions overlap more than 
expected by chance. Appropriate discrete model for sampling without replacement.

**Negative Control Validation**: Distinguishes genuine biological signal from statistical 
artifacts by comparing real correlation to shuffled null distribution.

*For detailed methodology documentation, see `cosmic/statistical_methodology.md`.*
"""


def _build_biological_bias_disclosure() -> str:
    """Build the biological bias disclosure section."""
    return """

### Biological Bias Disclosure

**Important limitations of COSMIC reference data:**

COSMIC (Catalogue of Somatic Mutations in Cancer) reflects inherent biases that affect 
comparisons with any dataset:

1. **Tumor Sampling Bias**: COSMIC over-represents commonly studied cancer types (e.g., 
   breast, lung, colorectal) and under-represents rare cancers. Fusion frequencies reflect 
   research attention, not true biological prevalence.

2. **Detection Technology Bias**: Historical data reflects older sequencing technologies 
   with different sensitivity profiles. Modern datasets may detect fusions missed in 
   earlier studies.

3. **Cohort Representation Bias**: COSMIC cohorts are not population-representative. 
   Certain demographics, geographic regions, and healthcare systems are over-represented.

4. **Publication Bias**: COSMIC aggregates published data, which skews toward positive 
   findings and known driver fusions.

**Interpretation guidance**: Statistical comparisons evaluate *plausibility and consistency* 
with COSMIC patterns, not *exact biological equivalence*. Strong agreement suggests the 
dataset contains biologically realistic fusion patterns; weak agreement does not necessarily 
indicate data quality problems.
"""


def _build_mock_cosmic_compatibility_statement() -> str:
    """Build the mock vs real COSMIC compatibility statement."""
    return """

### Mock COSMIC Compatibility Statement

**Note:** This validation used a synthetic mock COSMIC dataset for reference.

**Mock COSMIC is designed to preserve:**
- Heavy tail recurrence distribution (power-law/Zipf-like)
- Driver gene enrichment structure (BCR-ABL1, EML4-ALK, etc.)
- Non-uniform recurrence behavior
- Realistic gene pair counts (50+ fusion pairs)

**Mock COSMIC does NOT replicate:**
- Real COSMIC cohort structure
- Actual variant frequencies from patient samples
- Sample population demographics
- Temporal sampling patterns
- Complete fusion catalog coverage

**Appropriate use:** Pipeline testing, development, and methodology validation.

**Inappropriate use:** Drawing scientific conclusions about real fusion biology.

*When real COSMIC data is available, re-run validation with user-provided COSMIC reference 
for scientifically meaningful results.*
"""


def _section_5_cosmic(diag: Optional[Dict[str, Any]]) -> str:
    """Build Section 5: COSMIC Cross-Validation."""
    cosmic = (diag or {}).get("cosmic") or {}
    ours = cosmic.get("total_fusions_ours", "N/A")
    cosmic_count = cosmic.get("total_fusions_cosmic", "N/A")
    overlap = cosmic.get("overlap_count", "N/A")
    only_ours = cosmic.get("only_in_ours_count", "N/A")
    only_cosmic = cosmic.get("only_in_cosmic_count", "N/A")
    
    # Statistical metrics
    spearman_rho = cosmic.get("spearman_rho")
    spearman_p = cosmic.get("spearman_p_value")
    rho_ci_lower = cosmic.get("rho_ci_lower")
    rho_ci_upper = cosmic.get("rho_ci_upper")
    bootstrap_iterations = cosmic.get("bootstrap_iterations")
    top_overlap = cosmic.get("top_fusion_overlap")
    enrichment_ratio = cosmic.get("top_fusion_enrichment_ratio")
    enrichment_p = cosmic.get("enrichment_p_value")
    expected_overlap = cosmic.get("expected_overlap_random")
    negative_control_rho = cosmic.get("negative_control_rho")
    negative_control_p = cosmic.get("negative_control_p_value")
    validation_score = cosmic.get("cosmic_validation_score")
    validation_classification = cosmic.get("cosmic_validation_classification", "N/A")
    reference_source = cosmic.get("cosmic_reference_source", "N/A")
    reference_version = cosmic.get("cosmic_reference_version", "N/A")
    reference_hash = cosmic.get("cosmic_reference_file_hash", "N/A")
    reference_timestamp = cosmic.get("cosmic_reference_load_timestamp", "N/A")
    
    # Score component breakdown
    score_breakdown = cosmic.get("score_component_breakdown", {})
    reproducibility_lock = cosmic.get("reproducibility_lock", {})
    
    def fmt_stat(v):
        if v is None:
            return "N/A"
        try:
            if isinstance(v, float):
                return f"{v:.4f}"
            return str(v)
        except (TypeError, ValueError):
            return str(v)
    
    def fmt_p_value(v):
        if v is None:
            return "N/A"
        try:
            if isinstance(v, float):
                return f"{v:.6f}"
            return str(v)
        except (TypeError, ValueError):
            return str(v)
    
    # Scientific interpretation text
    interpretation_parts = []
    
    if spearman_rho is not None and isinstance(spearman_rho, (int, float)):
        rho_val = float(spearman_rho)
        abs_rho = abs(rho_val)
        if abs_rho < 0.3:
            strength = "weak"
        elif abs_rho < 0.7:
            strength = "moderate"
        else:
            strength = "strong"
        
        interpretation_parts.append(
            f"A Spearman correlation of {rho_val:.4f} indicates "
            f"{strength} rank-order agreement with COSMIC recurrence patterns."
        )
    
    if enrichment_p is not None and isinstance(enrichment_p, (int, float)):
        enrich_p_val = float(enrichment_p)
        if enrich_p_val < 0.05:
            interpretation_parts.append(
                f"Hypergeometric enrichment test (p={enrich_p_val:.6f}) demonstrates "
                f"statistically significant overlap of top-ranked fusions with COSMIC, "
                f"suggesting biological relevance beyond chance expectation."
            )
    
    if negative_control_rho is not None and spearman_rho is not None:
        real_abs = abs(float(spearman_rho))
        shuffled_abs = abs(float(negative_control_rho))
        if real_abs > shuffled_abs + 0.1:
            interpretation_parts.append(
                f"Negative control validation confirms that observed agreement "
                f"(ρ={real_abs:.4f}) significantly exceeds random expectation "
                f"(shuffled ρ={shuffled_abs:.4f}), supporting biological validity."
            )
    
    interpretation = ""
    if interpretation_parts:
        interpretation = "\n\n### Biological Interpretation\n\n" + "\n\n".join(interpretation_parts)
    
    # Build comprehensive metrics table with bootstrap CI
    ci_str = "N/A"
    if rho_ci_lower is not None and rho_ci_upper is not None:
        ci_str = f"[{rho_ci_lower:.4f}, {rho_ci_upper:.4f}]"
    
    table_rows = [
        f"| **Spearman rho** | {fmt_stat(spearman_rho)} |",
        f"| **Spearman 95% CI** | {ci_str} |",
        f"| **Spearman p-value** | {fmt_p_value(spearman_p)} |",
        f"| **Enrichment p-value** | {fmt_p_value(enrichment_p)} |",
        f"| **Expected Random Overlap** | {fmt_stat(expected_overlap)} |",
        f"| **Negative Control rho** | {fmt_stat(negative_control_rho)} |",
        f"| **COSMIC Validation Score** | {fmt_stat(validation_score)} |",
    ]
    
    table_body = "\n".join(table_rows)
    
    # Add classification section
    classification_section = ""
    if validation_classification != "N/A":
        classification_section = f"\n\n### COSMIC Validation Classification\n\n**{validation_classification}**\n"

    # Provenance section
    provenance_section = ""
    if reference_source != "N/A" or reference_version != "N/A":
        provenance_section = "\n\n### COSMIC Reference Provenance\n\n"
        provenance_rows = []
        if reference_source != "N/A":
            provenance_rows.append(f"| **Source** | {reference_source} |")
        if reference_version != "N/A":
            provenance_rows.append(f"| **Version** | {reference_version} |")
        if reference_hash != "N/A" and reference_hash:
            hash_short = str(reference_hash)[:16] + "..." if len(str(reference_hash)) > 16 else str(reference_hash)
            provenance_rows.append(f"| **File Hash** | `{hash_short}` |")
        if reference_timestamp != "N/A":
            provenance_rows.append(f"| **Load Timestamp** | {reference_timestamp} |")
        if provenance_rows:
            provenance_section += "| Field | Value |\n|-------|-------|\n" + "\n".join(provenance_rows)
    
    # Score component breakdown section
    breakdown_section = ""
    if score_breakdown:
        corr_comp = score_breakdown.get("correlation_component", "N/A")
        enrich_comp = score_breakdown.get("enrichment_component", "N/A")
        neg_ctrl_comp = score_breakdown.get("negative_control_component", "N/A")
        overlap_comp = score_breakdown.get("overlap_component", "N/A")
        
        def fmt_comp(v):
            if v is None or v == "N/A":
                return "N/A"
            try:
                return f"{float(v):.4f}"
            except (TypeError, ValueError):
                return str(v)
        
        breakdown_section = f"""

### Score Component Breakdown

| Component | Score |
|-----------|-------|
| **Correlation (rho + significance)** | {fmt_comp(corr_comp)} |
| **Enrichment** | {fmt_comp(enrich_comp)} |
| **Negative Control** | {fmt_comp(neg_ctrl_comp)} |
| **Overlap** | {fmt_comp(overlap_comp)} |

*Components sum to produce the final COSMIC Validation Score.*
"""
    
    # Method justification section
    method_justification = _build_method_justification_section()
    
    # Biological bias disclosure section
    bias_disclosure = _build_biological_bias_disclosure()
    
    # Mock vs Real COSMIC compatibility statement
    compatibility_statement = ""
    if reference_source in ("mock", "mock_fallback", "mock_v1"):
        compatibility_statement = _build_mock_cosmic_compatibility_statement()
    
    # Reproducibility lock section
    reproducibility_section = ""
    if reproducibility_lock:
        repro_rows = []
        if reproducibility_lock.get("python_version"):
            repro_rows.append(f"| **Python Version** | {reproducibility_lock.get('python_version')} |")
        if reproducibility_lock.get("numpy_version"):
            repro_rows.append(f"| **NumPy Version** | {reproducibility_lock.get('numpy_version')} |")
        if reproducibility_lock.get("scipy_version"):
            repro_rows.append(f"| **SciPy Version** | {reproducibility_lock.get('scipy_version')} |")
        if reproducibility_lock.get("pandas_version"):
            repro_rows.append(f"| **Pandas Version** | {reproducibility_lock.get('pandas_version')} |")
        if reproducibility_lock.get("cosmic_validation_code_version"):
            repro_rows.append(f"| **Validation Code Version** | {reproducibility_lock.get('cosmic_validation_code_version')} |")
        if reproducibility_lock.get("random_seed_mock_generation") is not None:
            repro_rows.append(f"| **Mock Generation Seed** | {reproducibility_lock.get('random_seed_mock_generation')} |")
        
        if repro_rows:
            reproducibility_section = "\n\n### Reproducibility Lock\n\n| Field | Value |\n|-------|-------|\n" + "\n".join(repro_rows)
    
    return f"""## 5. COSMIC Cross-Validation

### Statistical Metrics

| Metric | Value |
|--------|-------|
{table_body}
{classification_section}
{breakdown_section}
### Summary Statistics

| Statistic | Value |
|-----------|-------|
| **Total fusion pairs (ours)** | {ours} |
| **Total fusion pairs (COSMIC reference)** | {cosmic_count} |
| **Overlapping pairs** | {overlap} |
| **Only in our dataset** | {only_ours} |
| **Only in COSMIC** | {only_cosmic} |
{interpretation}
{method_justification}
{bias_disclosure}
{compatibility_statement}
{provenance_section}
{reproducibility_section}

*Statistical metrics and biological interpretation provided for COSMIC cross-validation assessment.*
"""


def generate_narrative_report(
    output_dir: Path,
    dataset_stem: str = "week2",
) -> Optional[Path]:
    """
    Generate a researcher-friendly Markdown report from Week 2 output files.

    Loads status, certification, and diagnostic JSON from output_dir and produces
    Statistical_Integrity_Report_{stem}.md.

    Args:
        output_dir: Directory containing the JSON output files.
        dataset_stem: Stem from input filename (e.g. "demo_fusion" for demo_fusion.parquet).

    Returns:
        Path to the generated report file, or None if generation failed.
    """
    output_dir = Path(output_dir).resolve()
    if not output_dir.is_dir():
        return None

    status_path = output_dir / STATUS_FILENAME_PATTERN.format(stem=dataset_stem)
    cert_path = output_dir / CERTIFICATION_FILENAME_PATTERN.format(stem=dataset_stem)
    diag_path = output_dir / DIAGNOSTIC_RESULTS_FILENAME_PATTERN.format(stem=dataset_stem)

    status = _load_json(status_path)
    if not status:
        return None

    certification = _load_json(cert_path)
    diag = _load_json(diag_path)
    notes = status.get("notes") or []

    sections = []
    sections.append("# Forensic Data Integrity & Statistical Validation Report\n")
    sections.append("*Data Integrity & Statistical Validation Layer — Week 2*\n")
    sections.append("---\n")
    sections.append(_section_1_executive_summary(status, certification))
    sections.append(_section_2_distribution(status, diag, output_dir))
    sections.append(_section_3_statistical_diagnostics(diag))
    sections.append(_section_4_benford(diag, notes))
    sections.append(_section_5_cosmic(diag))

    report_path = output_dir / REPORT_FILENAME_PATTERN.format(stem=dataset_stem)
    try:
        with open(report_path, "w", encoding="utf-8") as f:
            f.write("\n".join(sections))
        return report_path
    except OSError:
        return None
