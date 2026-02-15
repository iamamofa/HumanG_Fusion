"""
Data Integrity & Statistical Validation — MultiQC Integration.

Writes validation data as multiple MultiQC custom content modules (*_mqc.yaml).
MultiQC automatically discovers all *_mqc.yaml files in the output directory.
Each module has its own section with plots and/or tables.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# MultiQC custom content filenames (must end with _mqc.yaml for discovery)
MQC_DISTRIBUTION = "validation_distribution_mqc.yaml"
MQC_BENFORD = "validation_benford_mqc.yaml"
MQC_LOGNORMALITY = "validation_lognormality_mqc.yaml"
MQC_COSMIC = "validation_cosmic_mqc.yaml"
MQC_GATES = "validation_gates_mqc.yaml"


def _safe_float(v: Any) -> Optional[float]:
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _distribution_bins_from_csv(output_dir: Path, dataset_stem: str, n_bins: int = 10) -> Optional[Dict[str, int]]:
    """Load cleaned dataset, bin protein_length into n_bins, return bin_name -> count."""
    from week2_validation.output_names import CLEANED_DATASET_FILENAME_PATTERN
    cleaned_path = output_dir / CLEANED_DATASET_FILENAME_PATTERN.format(stem=dataset_stem)
    if not cleaned_path.is_file():
        return None
    try:
        import csv
        lengths: List[float] = []
        with open(cleaned_path, "r", encoding="utf-8") as f:
            r = csv.DictReader(f)
            if "protein_length" not in (r.fieldnames or []):
                return None
            for row in r:
                v = _safe_float(row.get("protein_length"))
                if v is not None and v > 0:
                    lengths.append(v)
        if not lengths:
            return None
        lo, hi = min(lengths), max(lengths)
        if hi <= lo:
            hi = lo + 1
        bin_edges = [lo + (hi - lo) * i / n_bins for i in range(n_bins + 1)]
        counts = [0] * n_bins
        for x in lengths:
            for i in range(n_bins):
                if i < n_bins - 1:
                    if bin_edges[i] <= x < bin_edges[i + 1]:
                        counts[i] += 1
                        break
                else:
                    if bin_edges[i] <= x <= bin_edges[i + 1]:
                        counts[i] += 1
                        break
        return {f"Bin_{i+1}": counts[i] for i in range(n_bins)}
    except (OSError, KeyError):
        return None


def _write_distribution_mqc(
    output_dir: Path,
    diagnostic_data: Optional[Dict[str, Any]],
    dataset_stem: str,
    sample_id: str,
) -> Optional[Path]:
    """Write validation_distribution_mqc.yaml: bargraph (10 bins) + table (N, min, max, mean, median, skewness)."""
    dist = (diagnostic_data or {}).get("distribution") or {}
    dq = (diagnostic_data or {}).get("data_quality") or {}
    if not dist and not dq:
        return None

    # Table: N, min, max, mean, median, skewness (always written)
    n_val = dq.get("total_rows_original") or dq.get("rows_used_for_analysis")
    min_val = dist.get("minimum")
    max_val = dist.get("maximum")
    mean_val = dist.get("mean")
    median_val = dist.get("median")
    skew_val = dq.get("skewness") or dist.get("skewness")
    table_row = {
        "N": n_val if n_val is not None else "N/A",
        "min": round(float(min_val), 2) if _safe_float(min_val) is not None else "N/A",
        "max": round(float(max_val), 2) if _safe_float(max_val) is not None else "N/A",
        "mean": round(float(mean_val), 2) if _safe_float(mean_val) is not None else "N/A",
        "median": round(float(median_val), 2) if _safe_float(median_val) is not None else "N/A",
        "skewness": round(float(skew_val), 4) if _safe_float(skew_val) is not None else "N/A",
    }
    table_content = {
        "id": "validation_distribution_table",
        "section_name": "Distribution Analysis",
        "description": "Summary statistics for protein length (N, min, max, mean, median, skewness).",
        "plot_type": "table",
        "pconfig": {"id": "validation_distribution_table", "namespace": "Validation"},
        "headers": {
            "N": {"title": "N", "description": "Sample size"},
            "min": {"title": "Min", "description": "Minimum protein length"},
            "max": {"title": "Max", "description": "Maximum protein length"},
            "mean": {"title": "Mean", "description": "Mean protein length"},
            "median": {"title": "Median", "description": "Median protein length"},
            "skewness": {"title": "Skewness", "description": "Distribution skewness"},
        },
        "data": {sample_id: table_row},
    }

    # Bargraph: binned protein length distribution (separate file so MultiQC gets one plot per file)
    bin_data = _distribution_bins_from_csv(output_dir, dataset_stem, n_bins=10)
    if bin_data:
        bargraph_content = {
            "id": "validation_distribution_bargraph",
            "section_name": "Distribution Analysis",
            "description": "Protein length distribution binned into 10 bins.",
            "plot_type": "bargraph",
            "pconfig": {
                "id": "validation_distribution_bargraph",
                "title": "Validation: Protein length distribution",
                "ylab": "Count",
                "xlab": "Bin",
            },
            "data": {sample_id: bin_data},
        }
        path_bargraph = output_dir / MQC_DISTRIBUTION
        path_bargraph.write_text(_dump_mqc_yaml(bargraph_content), encoding="utf-8")
        path_table = output_dir / "validation_distribution_table_mqc.yaml"
        path_table.write_text(_dump_mqc_yaml(table_content), encoding="utf-8")
        return path_bargraph  # caller can also check path_table.exists()
    # No bin data: single file with table only
    path_table = output_dir / MQC_DISTRIBUTION
    path_table.write_text(_dump_mqc_yaml(table_content), encoding="utf-8")
    return path_table


def _dump_mqc_yaml(obj: Dict[str, Any]) -> str:
    """Dump dict to YAML string (MultiQC custom content format)."""
    import yaml
    return yaml.dump(obj, default_flow_style=False, allow_unicode=True, sort_keys=False)


def _write_benford_mqc(
    output_dir: Path,
    diagnostic_data: Optional[Dict[str, Any]],
    sample_id: str,
) -> Optional[Path]:
    """Write validation_benford_mqc.yaml: bargraph (observed vs expected) + table."""
    benford = (diagnostic_data or {}).get("benford") or {}
    if not benford:
        return None

    obs = benford.get("observed_frequencies") or {}
    exp = benford.get("expected_frequencies") or {}
    digits = [str(d) for d in range(1, 10)]
    # Keys in JSON may be int or str
    def _get_freq(dct: Dict, key: str) -> float:
        v = dct.get(key) or dct.get(int(key))
        return round(float(v), 4) if v is not None else 0.0
    observed_vals = {d: _get_freq(obs, d) for d in digits}
    expected_vals = {d: _get_freq(exp, d) for d in digits}

    # Bargraph: two "samples" = Observed and Expected
    bargraph_content = {
        "id": "validation_benford",
        "section_name": "Benford's Law Analysis",
        "description": "First-digit observed vs expected frequencies (Benford's Law).",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "validation_benford_bargraph",
            "title": "Validation: Benford first-digit frequencies",
            "ylab": "Frequency",
            "xlab": "Digit",
            "stacking": "group",
        },
        "data": {
            "Observed": observed_vals,
            "Expected": expected_vals,
        },
    }

    chi2 = benford.get("chi_squared_statistic")
    pval = benford.get("p_value")
    scale_span = benford.get("scale_span_orders_of_magnitude")
    applicability = benford.get("applicability")
    table_row = {
        "chi_squared": round(float(chi2), 4) if _safe_float(chi2) is not None else "N/A",
        "p_value": round(float(pval), 6) if _safe_float(pval) is not None else "N/A",
        "scale_span": round(float(scale_span), 2) if _safe_float(scale_span) is not None else "N/A",
        "applicability": "Yes" if applicability is True else "No" if applicability is False else "N/A",
    }
    table_content = {
        "id": "validation_benford_table",
        "section_name": "Benford's Law Summary",
        "description": "Chi-squared test and scale span.",
        "plot_type": "table",
        "pconfig": {"id": "validation_benford_table", "namespace": "Validation"},
        "headers": {
            "chi_squared": {"title": "Chi-squared", "description": "Goodness-of-fit statistic"},
            "p_value": {"title": "p-value", "description": "Chi-squared p-value"},
            "scale_span": {"title": "Scale span", "description": "Orders of magnitude"},
            "applicability": {"title": "Applicable", "description": "Benford applicable (scale span ≥ 2)"},
        },
        "data": {sample_id: table_row},
    }

    path = output_dir / MQC_BENFORD
    # MultiQC: one plot per file for custom content; use first file as section with plot, second as table-only section
    path.write_text(_dump_mqc_yaml(bargraph_content), encoding="utf-8")
    path_table = output_dir / "validation_benford_table_mqc.yaml"
    path_table.write_text(_dump_mqc_yaml(table_content), encoding="utf-8")
    return path


def _write_lognormality_mqc(
    output_dir: Path,
    diagnostic_data: Optional[Dict[str, Any]],
    sample_id: str,
) -> Optional[Path]:
    """Write validation_lognormality_mqc.yaml: table (KS, AD, scipy_available)."""
    logn = (diagnostic_data or {}).get("log_normality") or {}
    if not logn:
        return None

    ks = logn.get("ks_statistic")
    ad = logn.get("ad_statistic")
    scipy_avail = logn.get("scipy_available")
    table_row = {
        "ks_statistic": round(float(ks), 6) if _safe_float(ks) is not None else "N/A",
        "ad_statistic": round(float(ad), 6) if _safe_float(ad) is not None else "N/A",
        "scipy_available": "Yes" if scipy_avail else "No",
    }
    content = {
        "id": "validation_lognormality",
        "section_name": "Log-Normality Tests",
        "description": "Kolmogorov-Smirnov and Anderson-Darling statistics for log-normal fit.",
        "plot_type": "table",
        "pconfig": {"id": "validation_lognormality_table", "namespace": "Validation"},
        "headers": {
            "ks_statistic": {"title": "KS statistic", "description": "Kolmogorov-Smirnov"},
            "ad_statistic": {"title": "AD statistic", "description": "Anderson-Darling"},
            "scipy_available": {"title": "SciPy", "description": "SciPy available for AD"},
        },
        "data": {sample_id: table_row},
    }
    path = output_dir / MQC_LOGNORMALITY
    path.write_text(_dump_mqc_yaml(content), encoding="utf-8")
    return path


def _write_cosmic_mqc(
    output_dir: Path,
    diagnostic_data: Optional[Dict[str, Any]],
    sample_id: str,
) -> Optional[Path]:
    """Write validation_cosmic_mqc.yaml: table (spearman_rho, p-value, overlap_count, score, classification). Scatter if paired data available."""
    cosmic = (diagnostic_data or {}).get("cosmic") or {}
    if not cosmic:
        return None

    rho = cosmic.get("spearman_rho")
    pval = cosmic.get("spearman_p_value")
    overlap = cosmic.get("overlap_count")
    score = cosmic.get("cosmic_validation_score")
    classification = cosmic.get("cosmic_validation_classification") or "N/A"
    table_row = {
        "spearman_rho": round(float(rho), 4) if _safe_float(rho) is not None else "N/A",
        "spearman_p": round(float(pval), 6) if _safe_float(pval) is not None else "N/A",
        "overlap_count": int(overlap) if overlap is not None else "N/A",
        "validation_score": round(float(score), 4) if _safe_float(score) is not None else "N/A",
        "classification": str(classification),
    }
    content = {
        "id": "validation_cosmic",
        "section_name": "COSMIC Cross-Validation",
        "description": "Spearman rank correlation and validation classification vs COSMIC reference.",
        "plot_type": "table",
        "pconfig": {"id": "validation_cosmic_table", "namespace": "Validation"},
        "headers": {
            "spearman_rho": {"title": "Spearman rho", "description": "Rank correlation"},
            "spearman_p": {"title": "p-value", "description": "Spearman p-value"},
            "overlap_count": {"title": "Overlap", "description": "Overlapping fusion pairs"},
            "validation_score": {"title": "Score", "description": "Validation score"},
            "classification": {"title": "Classification", "description": "COSMIC agreement"},
        },
        "data": {sample_id: table_row},
    }
    path = output_dir / MQC_COSMIC
    path.write_text(_dump_mqc_yaml(content), encoding="utf-8")
    return path


def _write_gates_mqc(
    output_dir: Path,
    status_data: Dict[str, Any],
    sample_id: str,
) -> Optional[Path]:
    """Write validation_gates_mqc.yaml: table with color-coded PASS/WARN/FAIL per module."""
    quality_gates = status_data.get("quality_gates") or {}
    gates = quality_gates.get("gates", [])
    if not gates:
        return None

    # One row per module: module_name, status, metric_name, metric_value, reason
    rows: Dict[str, Dict[str, Any]] = {}
    for g in gates:
        mod = g.get("module_name", "unknown")
        status = g.get("status", "N/A")
        metric_name = g.get("metric_name", "")
        metric_value = g.get("metric_value")
        reason = (g.get("reason") or "")[:80]
        rows[mod] = {
            "Module": mod,
            "Status": status,
            "Metric": metric_name,
            "Value": round(float(metric_value), 4) if _safe_float(metric_value) is not None else "—",
            "Reason": reason,
        }

    content = {
        "id": "validation_gates",
        "section_name": "Quality Gate Summary",
        "description": "Per-module validation gate status (PASS/WARN/FAIL).",
        "plot_type": "table",
        "pconfig": {"id": "validation_gates_table", "namespace": "Validation"},
        "headers": {
            "Module": {"title": "Module", "description": "Validation module"},
            "Status": {
                "title": "Status",
                "description": "Gate result",
                "bgcols": {"PASS": "#d1e7dd", "WARN": "#fff3cd", "FAIL": "#f8d7da", "SKIPPED": "#e2e3e5", "NOT_APPLICABLE": "#e2e3e5"},
            },
            "Metric": {"title": "Metric", "description": "Primary metric"},
            "Value": {"title": "Value", "description": "Metric value"},
            "Reason": {"title": "Reason", "description": "Gate reason"},
        },
        "data": rows,
    }
    path = output_dir / MQC_GATES
    path.write_text(_dump_mqc_yaml(content), encoding="utf-8")
    return path


def _write_validation_mqc_yaml(
    output_dir: Path,
    status_data: Dict[str, Any],
    diagnostic_data: Optional[Dict[str, Any]] = None,
    dataset_stem: str = "validation",
) -> List[Path]:
    """
    Write multiple MultiQC custom content YAML files. MultiQC will automatically
    find all *_mqc.yaml files in the output directory.

    Writes:
      - validation_distribution_mqc.yaml (and optional _bargraph if bins available)
      - validation_benford_mqc.yaml + validation_benford_table_mqc.yaml
      - validation_lognormality_mqc.yaml
      - validation_cosmic_mqc.yaml
      - validation_gates_mqc.yaml

    Returns list of paths written.
    """
    output_dir = Path(output_dir)
    sample_id = dataset_stem or "validation"
    written: List[Path] = []

    # 1. Distribution (may write validation_distribution_mqc.yaml + validation_distribution_table_mqc.yaml)
    p = _write_distribution_mqc(output_dir, diagnostic_data, dataset_stem, sample_id)
    if p:
        written.append(p)
        dist_table = output_dir / "validation_distribution_table_mqc.yaml"
        if dist_table.exists():
            written.append(dist_table)

    # 2. Benford
    p = _write_benford_mqc(output_dir, diagnostic_data, sample_id)
    if p:
        written.append(p)
        table_path = output_dir / "validation_benford_table_mqc.yaml"
        if table_path.exists():
            written.append(table_path)

    # 3. Log-normality
    p = _write_lognormality_mqc(output_dir, diagnostic_data, sample_id)
    if p:
        written.append(p)

    # 4. COSMIC
    p = _write_cosmic_mqc(output_dir, diagnostic_data, sample_id)
    if p:
        written.append(p)

    # 5. Quality gates
    p = _write_gates_mqc(output_dir, status_data, sample_id)
    if p:
        written.append(p)

    return written


def run_multiqc(output_dir: Path, report_name: str = "multiqc_report") -> Tuple[Optional[Path], Optional[str]]:
    """
    Run MultiQC on the output directory.

    MultiQC automatically discovers all *_mqc.yaml (and *_mqc.yml, etc.) files
    in the directory, so no need to pass file names. Running `multiqc {output_dir}`
    picks up all validation modules and produces a report with multiple sections
    and interactive plots.

    Args:
        output_dir: Directory containing validation outputs and *_mqc.yaml files
        report_name: Base name for the report (default: multiqc_report)

    Returns:
        (Path to multiqc_report.html, error_message). On success: (path, None).
    """
    try:
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "multiqc",
                str(output_dir),
                "-o",
                str(output_dir),
                "-n",
                report_name,
                "--force",
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            err = (result.stderr or result.stdout or "").strip() or f"Exit code {result.returncode}"
            return (None, err)
        report_path = output_dir / f"{report_name}.html"
        if report_path.exists():
            return (report_path, None)
        return (None, "Report file not created")
    except FileNotFoundError:
        return (None, "MultiQC not found. Install with: pip install multiqc")
    except subprocess.TimeoutExpired:
        return (None, "MultiQC timed out")
    except OSError as e:
        return (None, str(e))


def generate_multiqc_report(
    output_dir: Path,
    status_path: Path,
    diagnostic_path: Optional[Path] = None,
    dataset_stem: str = "validation",
) -> Tuple[Optional[Path], Optional[str]]:
    """
    Write all MultiQC custom content YAML files and run MultiQC.

    Args:
        output_dir: Output directory (e.g. week2_validation/results/{folder}/)
        status_path: Path to validation_status_{stem}.json
        diagnostic_path: Path to validation_diagnostic_results_{stem}.json (optional)
        dataset_stem: Dataset stem for sample naming

    Returns:
        (Path to multiqc_report.html, error_message). On success: (path, None).
    """
    if not status_path.exists():
        return (None, "Status file not found")

    try:
        status_data = json.loads(status_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return (None, "Could not read status JSON")

    diagnostic_data = None
    if diagnostic_path and diagnostic_path.exists():
        try:
            diagnostic_data = json.loads(diagnostic_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            pass

    _write_validation_mqc_yaml(output_dir, status_data, diagnostic_data, dataset_stem)
    return run_multiqc(output_dir)
