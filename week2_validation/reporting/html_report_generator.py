"""
Data Integrity & Statistical Validation — HTML Report Generator.

Generates an HTML report with the same content as the PDF report.
Uses all information generated for the PDF. Does not modify the PDF or any other components.
"""

import json
import sys
from base64 import b64encode
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional


def _safe_format(value: Any, precision: int = 4) -> str:
    """Safely format a value for display."""
    if value is None:
        return "N/A"
    try:
        if isinstance(value, (int, float)):
            if abs(value) < 1e-10:
                return "0.0"
            return f"{float(value):.{precision}f}"
        return str(value)
    except (TypeError, ValueError):
        return str(value)


def _safe_format_correlation(value: Any) -> str:
    """Format correlation values with 3 decimal places."""
    if value is None:
        return "N/A"
    try:
        corr_val = float(value)
        if abs(corr_val) < 1e-10:
            return "0.000"
        return f"{corr_val:.3f}"
    except (TypeError, ValueError):
        return str(value)


def _safe_format_p_value(value: Any) -> str:
    """Format p-value with 4 decimal places or scientific notation."""
    if value is None:
        return "N/A"
    try:
        p_val = float(value)
        if p_val < 0.0001:
            return f"{p_val:.4e}"
        return f"{p_val:.4f}"
    except (TypeError, ValueError):
        return str(value)


def _load_json_file(path: Path) -> Optional[Dict[str, Any]]:
    """Load JSON file, return None if missing or invalid."""
    if not path.is_file():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


def _get_stability_badge(ci_width: float) -> tuple[str, str]:
    """Get stability badge based on CI width."""
    if ci_width < 0.3:
        return ("&#x1F7E2;", "Stable")
    elif ci_width < 0.6:
        return ("&#x1F7E1;", "Moderate")
    else:
        return ("&#x1F534;", "Unstable")


def _status_css_class(status: str) -> str:
    """Get CSS class for status (MultiQC-style)."""
    s = str(status).upper()
    if "PASS" in s or "STRONG" in s or "HIGH" in s:
        return "status-pass"
    if "MODERATE" in s or "WARNING" in s:
        return "status-warning"
    if "FAIL" in s or "LIMITED" in s or "WEAK" in s or "LOW" in s:
        return "status-fail"
    return "status-neutral"


def _img_to_data_uri(path: Path, max_size_mb: float = 2.0) -> Optional[str]:
    """Embed image as data URI if file exists and is small enough."""
    if not path or not path.exists():
        return None
    try:
        size_mb = path.stat().st_size / (1024 * 1024)
        if size_mb > max_size_mb:
            return None
        data = path.read_bytes()
        b64 = b64encode(data).decode("ascii")
        suffix = path.suffix.lower()
        mime = "image/png" if suffix == ".png" else "image/jpeg" if suffix in (".jpg", ".jpeg") else "image/png"
        return f"data:{mime};base64,{b64}"
    except Exception:
        return None


def _html_escape(s: str) -> str:
    """Escape HTML special characters."""
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


# MultiQC-style CSS (orange #F18046, dark #160F26)
MULTIQC_CSS = """
:root {
  --mqc-orange: #F18046;
  --mqc-dark: #160F26;
  --mqc-gray: #666666;
  --mqc-light: #f5f5f5;
  --mqc-border: #e0e0e0;
  --mqc-pass: #e8f5e9;
  --mqc-warning: #fff3e0;
  --mqc-fail: #ffebee;
}
* { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  font-size: 14px;
  line-height: 1.5;
  color: var(--mqc-dark);
  margin: 0;
  padding: 0;
  background: #fff;
}
.mqc-header {
  background: linear-gradient(135deg, var(--mqc-dark) 0%, #2a1f3d 100%);
  color: #fff;
  padding: 24px 32px;
  margin-bottom: 0;
}
.mqc-header h1 {
  margin: 0;
  font-size: 24px;
  font-weight: 600;
}
.mqc-header .subtitle {
  margin-top: 8px;
  font-size: 13px;
  opacity: 0.9;
}
.mqc-meta {
  background: var(--mqc-light);
  border: 1px solid var(--mqc-border);
  padding: 12px 20px;
  margin: 20px 32px;
  border-radius: 4px;
  font-size: 13px;
}
.mqc-section {
  margin: 24px 32px 32px;
  border: 1px solid var(--mqc-border);
  border-radius: 6px;
  overflow: hidden;
}
.mqc-section-header {
  background: var(--mqc-light);
  padding: 12px 20px;
  cursor: pointer;
  font-weight: 600;
  font-size: 16px;
  border-bottom: 1px solid var(--mqc-border);
  display: flex;
  align-items: center;
  justify-content: space-between;
}
.mqc-section-header:hover { background: #eee; }
.mqc-section-header::after {
  content: "\\25BC";
  font-size: 12px;
  transition: transform 0.2s;
}
.mqc-section.collapsed .mqc-section-header::after { transform: rotate(-90deg); }
.mqc-section-content {
  padding: 20px;
  background: #fff;
}
.mqc-section.collapsed .mqc-section-content { display: none; }
.mqc-table {
  width: 100%;
  border-collapse: collapse;
  margin: 12px 0;
  font-size: 13px;
}
.mqc-table th, .mqc-table td {
  padding: 10px 14px;
  text-align: left;
  border: 1px solid var(--mqc-border);
}
.mqc-table th {
  background: var(--mqc-light);
  font-weight: 600;
}
.mqc-table tr:nth-child(even) td { background: #fafafa; }
.mqc-table .num { text-align: right; }
.mqc-status-block {
  padding: 14px 18px;
  margin: 12px 0;
  border-radius: 4px;
  border-left: 4px solid var(--mqc-orange);
}
.mqc-status-block.status-pass { background: var(--mqc-pass); border-left-color: #4caf50; }
.mqc-status-block.status-warning { background: var(--mqc-warning); border-left-color: var(--mqc-orange); }
.mqc-status-block.status-fail { background: var(--mqc-fail); border-left-color: #f44336; }
.mqc-status-block.status-neutral { background: var(--mqc-light); }
.mqc-highlight {
  background: var(--mqc-light);
  border: 1px solid var(--mqc-border);
  padding: 14px 18px;
  margin: 12px 0;
  border-radius: 4px;
}
.mqc-warning-box {
  background: #fff3cd;
  border: 1px solid #ffc107;
  padding: 14px 18px;
  margin: 12px 0;
  border-radius: 4px;
}
.mqc-img {
  max-width: 100%;
  height: auto;
  border: 1px solid var(--mqc-border);
  border-radius: 4px;
  margin: 12px 0;
}
.mqc-footer {
  margin-top: 32px;
  padding: 20px 32px;
  font-size: 12px;
  color: var(--mqc-gray);
  border-top: 1px solid var(--mqc-border);
}
.mqc-h2 { font-size: 18px; margin: 20px 0 12px; font-weight: 600; }
.mqc-h3 { font-size: 15px; margin: 16px 0 8px; font-weight: 600; }
/* Pipeline status bar (Data Provenance) */
.pipeline-bar { display: flex; flex-wrap: wrap; align-items: center; gap: 4px 0; margin: 16px 0; padding: 12px 0; overflow-x: auto; }
.pipeline-step { display: flex; align-items: center; flex-shrink: 0; }
.pipeline-step-dot { width: 14px; height: 14px; border-radius: 50%; margin-right: 6px; }
.pipeline-step-dot.status-pass { background: #4caf50; }
.pipeline-step-dot.status-approved { background: #4caf50; }
.pipeline-step-dot.status-warn, .pipeline-step-dot.status-conditional { background: #ff9800; }
.pipeline-step-dot.status-fail, .pipeline-step-dot.status-rejected { background: #f44336; }
.pipeline-step-dot.status-skipped, .pipeline-step-dot.status-neutral { background: #9e9e9e; }
.pipeline-step-label { font-size: 12px; white-space: nowrap; }
.pipeline-arrow { color: #999; margin: 0 6px; font-size: 10px; }
"""


def _build_html(
    status_data: dict,
    diagnostic_data: Optional[dict],
    cert_data: Optional[dict],
    output_dir: Path,
    dataset_stem: str,
    histogram_paths: Optional[Dict[str, Optional[Path]]],
    interpretation_text: Optional[str],
    skewness_path: Optional[Path],
    scientific_narrative: Optional[str],
    confidence_statement: Optional[str],
    limitations_text: Optional[str],
    robustness_data: Optional[dict],
    null_model_data: Optional[dict],
    external_validity_data: Optional[dict],
    bias_analysis_data: Optional[dict],
    claim_strength_data: Optional[dict],
    cosmic_version: str,
) -> str:
    """Build full HTML report content."""
    run_timestamp = status_data.get("run_timestamp_utc", "N/A")
    dataset_hash = status_data.get("dataset_hash", "N/A")
    exit_code = status_data.get("exit_code", -1)
    dq = status_data.get("data_quality", {})

    # Dashboard statuses (same logic as PDF)
    data_integrity_status = "PASS" if exit_code == 0 else "FAIL"
    distribution_status = "PASS"
    if diagnostic_data:
        dist = diagnostic_data.get("distribution", {})
        if not dist or len(dist) == 0:
            distribution_status = "LIMITED"

    benford_status = "N/A"
    if diagnostic_data:
        benford = diagnostic_data.get("benford", {})
        applicability = benford.get("applicability")
        benford_status = "APPLICABLE" if applicability is True else "NOT APPLICABLE"

    cosmic_status = "N/A"
    cosmic_confidence = "N/A"
    if diagnostic_data:
        cosmic = diagnostic_data.get("cosmic", {})
        validation_class = cosmic.get("cosmic_validation_classification", "N/A")
        if validation_class != "N/A":
            cosmic_status = validation_class
            if "STRONG" in str(validation_class).upper():
                cosmic_confidence = "HIGH"
            elif "MODERATE" in str(validation_class).upper():
                cosmic_confidence = "MODERATE"
            else:
                cosmic_confidence = "LOW"

    external_validity_status = "N/A"
    if external_validity_data:
        ci_lower = external_validity_data.get("stability_ci_lower")
        ci_upper = external_validity_data.get("stability_ci_upper")
        if ci_lower is not None and ci_upper is not None:
            ci_width = ci_upper - ci_lower
            external_validity_status = "STABLE" if ci_width < 0.3 else "MODERATE" if ci_width < 0.6 else "LIMITED"

    overall_status = "N/A"
    overall_confidence = "N/A"
    if claim_strength_data:
        overall_status = claim_strength_data.get("classification", "N/A")
        overall_confidence = claim_strength_data.get("confidence_level", "N/A")

    # Resolve image paths relative to HTML output (assume HTML is in output_dir root)
    def rel_img(p: Optional[Path]) -> Optional[str]:
        if not p or not p.exists():
            return None
        try:
            return str(p.relative_to(output_dir)).replace("\\", "/")
        except ValueError:
            return str(p).replace("\\", "/")

    skewness_src = None
    if skewness_path and skewness_path.exists():
        data_uri = _img_to_data_uri(skewness_path)
        skewness_src = data_uri or rel_img(skewness_path)

    linear_hist_src = None
    log_hist_src = None
    if histogram_paths:
        linear_path = histogram_paths.get("linear")
        log_path = histogram_paths.get("log")
        if linear_path and linear_path.exists():
            linear_hist_src = _img_to_data_uri(linear_path) or rel_img(linear_path)
        if log_path and log_path.exists():
            log_hist_src = _img_to_data_uri(log_path) or rel_img(log_path)

    bootstrap_plot_src = None
    if robustness_data:
        bp = robustness_data.get("bootstrap_plot_path")
        if bp:
            bp_path = Path(bp) if isinstance(bp, str) else bp
            if bp_path.exists():
                bootstrap_plot_src = _img_to_data_uri(bp_path) or rel_img(bp_path)

    benford_img_path = output_dir / dataset_stem / "benford_analysis.png"
    if not benford_img_path.exists():
        benford_img_path = output_dir / "benford_analysis.png"
    benford_src = _img_to_data_uri(benford_img_path) if benford_img_path.exists() else None
    if not benford_src and benford_img_path.exists():
        benford_src = rel_img(benford_img_path)

    sections = []

    # Data Provenance: interactive pipeline status bar (horizontal steps with colored dots)
    provenance_chain = (diagnostic_data or {}).get("provenance") or []
    def _dot_class(st: str) -> str:
        u = (st or "").upper()
        if u in ("PASS", "APPROVED"):
            return "status-pass"
        if u in ("WARN", "CONDITIONAL"):
            return "status-warn"
        if u in ("FAIL", "REJECTED"):
            return "status-fail"
        return "status-neutral"
    pipeline_items = []
    for i, s in enumerate(provenance_chain):
        step_name = _html_escape(s.get("step", ""))
        status = (s.get("status") or "N/A").upper()
        details = _html_escape((s.get("details") or "")[:80])
        dot_class = _dot_class(status)
        title = f"{step_name}: {status}" + (f" — {details}" if details else "")
        pipeline_items.append(
            f'<span class="pipeline-step" title="{title}">'
            f'<span class="pipeline-step-dot {dot_class}" aria-label="{status}"></span>'
            f'<span class="pipeline-step-label">{step_name}</span></span>'
        )
    pipeline_html = ""
    if pipeline_items:
        pipeline_html = "<div class=\"pipeline-bar\">" + " <span class=\"pipeline-arrow\">&#8594;</span> ".join(pipeline_items) + "</div>"
    else:
        pipeline_html = "<p>Provenance chain not available (run diagnostics to populate).</p>"
    sections.append(f"""
<div class="mqc-section" id="sec-provenance">
  <div class="mqc-section-header">Data Provenance</div>
  <div class="mqc-section-content">
    <p>Chain of custody from raw input to validated output. Hover over a step for details.</p>
    {pipeline_html}
  </div>
</div>""")

    # Section 1: General Statistics (MultiQC-style)
    sections.append(f"""
<div class="mqc-section" id="sec-general">
  <div class="mqc-section-header">General Statistics</div>
  <div class="mqc-section-content">
    <table class="mqc-table">
      <thead><tr><th>QC Layer</th><th>Status</th><th>Confidence</th></tr></thead>
      <tbody>
        <tr><td>Data Integrity</td><td class="{_status_css_class(data_integrity_status)}">{data_integrity_status}</td><td>{"HIGH" if data_integrity_status == "PASS" else "LOW"}</td></tr>
        <tr><td>Distribution Validity</td><td class="{_status_css_class(distribution_status)}">{distribution_status}</td><td>{"HIGH" if distribution_status == "PASS" else "MODERATE"}</td></tr>
        <tr><td>Benford Applicability</td><td>{_html_escape(benford_status)}</td><td>N/A</td></tr>
        <tr><td>COSMIC Biological Agreement</td><td class="{_status_css_class(cosmic_status)}">{_html_escape(cosmic_status)}</td><td>{_html_escape(cosmic_confidence)}</td></tr>
        <tr><td>External Validity</td><td class="{_status_css_class(external_validity_status)}">{external_validity_status}</td><td>{"HIGH" if external_validity_status == "STABLE" else "MODERATE" if external_validity_status == "MODERATE" else "LOW"}</td></tr>
        <tr><td>Overall Scientific Claim Strength</td><td class="{_status_css_class(overall_status)}">{_html_escape(overall_status)}</td><td>{_html_escape(overall_confidence)}</td></tr>
      </tbody>
    </table>
    <div class="mqc-highlight">
      <strong>Interpretation Rules</strong><br/>
      &bull; Effect size evaluated before p-value<br/>
      &bull; Low overlap reduces statistical power<br/>
      &bull; Wide CI = unstable biological inference
    </div>
  </div>
</div>""")

    # Quality Gates (formal PASS/WARN/FAIL and overall approval)
    quality_gates = status_data.get("quality_gates") or {}
    gates_list = quality_gates.get("gates", [])
    overall_approval = quality_gates.get("overall_approval", "N/A")
    def _gate_badge(s: str) -> str:
        u = (s or "").upper()
        if u == "PASS":
            return '<span class="status-pass">&#x1F7E2; PASS</span>'
        if u == "WARN":
            return '<span class="status-warning">&#x1F7E1; WARN</span>'
        if u == "FAIL":
            return '<span class="status-fail">&#x1F534; FAIL</span>'
        if u in ("SKIPPED", "NOT_APPLICABLE"):
            return '<span class="status-neutral">&#x26AA; ' + _html_escape(s or "SKIPPED") + "</span>"
        return _html_escape(str(s))
    qg_rows = ""
    for g in gates_list:
        mod = _html_escape(str(g.get("module_name", "")))
        status = _gate_badge(g.get("status", ""))
        metric = _html_escape(str(g.get("metric_name", "")))
        val = g.get("metric_value")
        thresh = g.get("threshold")
        _r = (g.get("reason") or "")
        reason = _html_escape(_r[:100] + ("..." if len(_r) > 100 else ""))
        val_str = f"{float(val):.4f}" if val is not None else "—"
        thresh_str = f"{float(thresh):.4f}" if thresh is not None else "—"
        qg_rows += f"<tr><td>{mod}</td><td>{status}</td><td>{metric}</td><td class=\"num\">{val_str}</td><td class=\"num\">{thresh_str}</td><td>{reason}</td></tr>"
    overall_class = _status_css_class(overall_approval) if overall_approval != "N/A" else "status-neutral"
    qg_section = f"""
<div class="mqc-section" id="sec-quality-gates">
  <div class="mqc-section-header">Quality Gates</div>
  <div class="mqc-section-content">
    <p>Per-module verdicts and overall dataset approval. <strong>Overall:</strong> <span class="{overall_class}">{_html_escape(overall_approval)}</span></p>
    <table class="mqc-table">
      <thead><tr><th>Module</th><th>Status</th><th>Metric</th><th>Value</th><th>Threshold</th><th>Reason</th></tr></thead>
      <tbody>{qg_rows if qg_rows else '<tr><td colspan="6">No quality gate results available.</td></tr>'}</tbody>
    </table>
    <div class="mqc-highlight">
      <strong>Rules:</strong> APPROVED = no FAIL, &ge;2 PASS; CONDITIONAL = 1+ WARN, no FAIL; REJECTED = any FAIL.
    </div>
  </div>
</div>"""
    sections.append(qg_section)

    # Section 2: Data Integrity
    schema_status = "PASSED" if exit_code == 0 else "FAILED"
    excluded_fraction = dq.get("excluded_fraction", 0.0)
    warning_flag = dq.get("warning_flag", False)
    high_risk_flag = dq.get("high_risk_flag", False)
    total_rows = dq.get("total_rows_original", "N/A")
    rows_used = dq.get("rows_used_for_analysis", "N/A")

    sections.append(f"""
<div class="mqc-section" id="sec-integrity">
  <div class="mqc-section-header">2. Data Integrity Validation</div>
  <div class="mqc-section-content">
    <div class="mqc-status-block {_status_css_class(data_integrity_status)}">
      <strong>DATA INTEGRITY STATUS</strong><br/>
      Schema: {schema_status}. Rows: {total_rows} total, {rows_used} used.
    </div>
    <table class="mqc-table">
      <thead><tr><th>Metric</th><th>Value</th><th>Interpretation</th></tr></thead>
      <tbody>
        <tr><td>Schema Validation</td><td>{schema_status}</td><td>{"PASS" if exit_code == 0 else "FAIL"}</td></tr>
        <tr><td>Excluded Fraction</td><td class="num">{_safe_format(excluded_fraction * 100, 2)}%</td><td>{"Low exclusion" if excluded_fraction < 0.1 else "Moderate exclusion" if excluded_fraction < 0.3 else "High exclusion"}</td></tr>
        <tr><td>Warning Flag</td><td>{"True" if warning_flag else "False"}</td><td>{"Data quality concerns" if warning_flag else "No warnings"}</td></tr>
        <tr><td>High Risk Flag</td><td>{"True" if high_risk_flag else "False"}</td><td>{"Critical issues detected" if high_risk_flag else "No critical issues"}</td></tr>
      </tbody>
    </table>
  </div>
</div>""")

    # Section 3: Distribution & Statistical Testing
    dist_html = ""
    if diagnostic_data:
        dist = diagnostic_data.get("distribution", {})
        log_norm = diagnostic_data.get("log_normality", {})
        dist_html = f"""
    <table class="mqc-table">
      <thead><tr><th>Metric</th><th>Value</th><th>Interpretation</th></tr></thead>
      <tbody>
        <tr><td>Mean</td><td class="num">{_safe_format(dist.get("mean"))}</td><td>Central tendency</td></tr>
        <tr><td>Median</td><td class="num">{_safe_format(dist.get("median"))}</td><td>Robust central value</td></tr>
        <tr><td>Std Dev</td><td class="num">{_safe_format(dist.get("std"))}</td><td>Dispersion measure</td></tr>
        <tr><td>Skewness</td><td class="num">{_safe_format(dist.get("skewness"))}</td><td>Distribution asymmetry</td></tr>
        <tr><td>Minimum</td><td class="num">{_safe_format(dist.get("minimum"))}</td><td>Lower bound</td></tr>
        <tr><td>Maximum</td><td class="num">{_safe_format(dist.get("maximum"))}</td><td>Upper bound</td></tr>
      </tbody>
    </table>
    <table class="mqc-table">
      <thead><tr><th>Test</th><th>Statistic</th><th>Interpretation</th></tr></thead>
      <tbody>
        <tr><td>KS Test</td><td class="num">{_safe_format(log_norm.get("ks_statistic"), 6) if log_norm.get("ks_statistic") is not None else "N/A"}</td><td>Smaller = better log-normal fit</td></tr>
        <tr><td>AD Test</td><td class="num">{_safe_format(log_norm.get("ad_statistic"), 6) if log_norm.get("ad_statistic") is not None else "N/A"}</td><td>Weighted tail assessment</td></tr>
      </tbody>
    </table>"""
        if interpretation_text:
            dist_html += f'<p><strong>Interpretation:</strong> {_html_escape(interpretation_text[:500])}{"..." if len(interpretation_text) > 500 else ""}</p>'

    skewness_img_html = ""
    if skewness_src:
        skewness_img_html = f'<p><strong>Skewness Diagnostic</strong></p><img class="mqc-img" src="{skewness_src}" alt="Skewness" style="max-width:800px;">'

    hist_html = ""
    if linear_hist_src:
        hist_html += f'<p><strong>Protein Length Distribution Histogram</strong></p><img class="mqc-img" src="{linear_hist_src}" alt="Histogram" style="max-width:700px;">'
    if log_hist_src:
        hist_html += f'<p><strong>Protein Length Distribution (Log Scale)</strong></p><img class="mqc-img" src="{log_hist_src}" alt="Log Histogram" style="max-width:700px;">'
    if not hist_html:
        hist_html = "<p>Histogram visualizations unavailable.</p>"

    sections.append(f"""
<div class="mqc-section" id="sec-distribution">
  <div class="mqc-section-header">3. Distribution &amp; Statistical Testing</div>
  <div class="mqc-section-content">
    <div class="mqc-h3">Distribution Diagnostics</div>
    {dist_html or "<p>No statistical validation data available.</p>"}
    {skewness_img_html}
    {hist_html}
  </div>
</div>""")

    # Section 4: Benford
    benford_html = "<p>No Benford analysis data available.</p>"
    if diagnostic_data:
        benford = diagnostic_data.get("benford", {})
        if benford:
            applicability = benford.get("applicability")
            reason = benford.get("reason_if_not_applicable", "")
            scale_span = benford.get("scale_span_orders_of_magnitude")
            chi_sq = benford.get("chi_squared_statistic")
            df_val = benford.get("degrees_of_freedom")
            p_val = benford.get("p_value")
            p_note = ""
            if p_val is not None and isinstance(p_val, (int, float)):
                try:
                    p_note = " p &gt; 0.05 indicates consistency with Benford distribution." if float(p_val) > 0.05 else " p &lt; 0.05 indicates deviation from Benford distribution (diagnostic only, not inferential)."
                except (TypeError, ValueError):
                    pass
            benford_html = f"""
    <p><strong>Benford Applicable:</strong> {applicability}</p>
    {f'<p><strong>Reason:</strong> {_html_escape(reason)}</p>' if reason else ''}
    {f'<p><strong>Scale Span:</strong> {_safe_format(scale_span)} orders of magnitude</p>' if scale_span is not None else ''}
    <table class="mqc-table">
      <thead><tr><th>Metric</th><th class="num">Value</th></tr></thead>
      <tbody>
        {f'<tr><td>Chi-squared statistic</td><td class="num">{_safe_format(chi_sq, 4)}</td></tr>' if chi_sq is not None else ''}
        {f'<tr><td>Degrees of freedom</td><td class="num">{df_val}</td></tr>' if df_val is not None else ''}
        {f'<tr><td>p-value</td><td class="num">{_safe_format_p_value(p_val)}</td></tr>' if p_val is not None else ''}
      </tbody>
    </table>
    {f'<p><em>Diagnostic note (descriptive only):</em>{_html_escape(p_note)}</p>' if p_note else ''}
    <table class="mqc-table">
      <thead><tr><th>Digit</th><th class="num">Observed</th><th class="num">Expected</th></tr></thead>
      <tbody>"""
            obs_freq = benford.get("observed_frequencies", {})
            exp_freq = benford.get("expected_frequencies", {})
            for d in range(1, 10):
                obs = obs_freq.get(d, obs_freq.get(str(d), 0))
                exp = exp_freq.get(d, exp_freq.get(str(d), 0))
                benford_html += f"<tr><td>{d}</td><td class='num'>{_safe_format(obs, 4)}</td><td class='num'>{_safe_format(exp, 4)}</td></tr>"
            benford_html += "</tbody></table>"
            if benford_src:
                benford_html += f'<img class="mqc-img" src="{benford_src}" alt="Benford" style="max-width:600px;">'

    sections.append(f"""
<div class="mqc-section" id="sec-benford">
  <div class="mqc-section-header">4. Benford Analysis</div>
  <div class="mqc-section-content">
    <p>In many real-world datasets, numbers starting with 1 appear more frequently than those starting with 9. This test only applies to datasets spanning multiple orders of magnitude.</p>
    {benford_html}
  </div>
</div>""")

    # Section 5: COSMIC
    cosmic_html = "<p>No COSMIC validation data available.</p>"
    if diagnostic_data:
        cosmic = diagnostic_data.get("cosmic", {})
        if cosmic:
            overlap_count = cosmic.get("overlap_count", 0)
            cosmic_html = f"""
    <table class="mqc-table">
      <thead><tr><th>Statistic</th><th>Value</th><th>Interpretation</th></tr></thead>
      <tbody>
        <tr><td>Total Fusions (Ours)</td><td>{cosmic.get("total_fusions_ours", "N/A")}</td><td>Dataset size</td></tr>
        <tr><td>Total Fusions (COSMIC)</td><td>{cosmic.get("total_fusions_cosmic", "N/A")}</td><td>Reference database size</td></tr>
        <tr><td>Overlap Count</td><td>{cosmic.get("overlap_count", "N/A")}</td><td>Common fusion pairs</td></tr>
        <tr><td>Spearman Rho</td><td class="num">{_safe_format_correlation(cosmic.get("spearman_rho"))}</td><td>Rank agreement</td></tr>
        <tr><td>Spearman P Value</td><td class="num">{_safe_format_p_value(cosmic.get("spearman_p_value"))}</td><td>Significance</td></tr>
        <tr><td>Enrichment P Value</td><td class="num">{_safe_format_p_value(cosmic.get("enrichment_p_value"))}</td><td>Top fusion overlap</td></tr>
      </tbody>
    </table>"""
            if isinstance(overlap_count, (int, float)) and overlap_count < 10:
                cosmic_html = f'<div class="mqc-warning-box">COSMIC overlap &lt; 10 fusion pairs (n={int(overlap_count)}). Statistical uncertainty expected.</div>{cosmic_html}'
            validation_class = cosmic.get("cosmic_validation_classification", "N/A")
            if validation_class != "N/A":
                cosmic_html += f'<div class="mqc-status-block {_status_css_class(validation_class)}"><strong>COSMIC VALIDATION:</strong> {_html_escape(validation_class)}</div>'

    sections.append(f"""
<div class="mqc-section" id="sec-cosmic">
  <div class="mqc-section-header">5. COSMIC Cross-Validation</div>
  <div class="mqc-section-content">
    {cosmic_html}
  </div>
</div>""")

    # Section 6: Robustness
    robustness_html = "<p>Robustness analysis unavailable due to insufficient COSMIC overlap.</p>"
    if robustness_data:
        stability_df = robustness_data.get("stability_table")
        robustness_html = ""
        if stability_df is not None and len(stability_df) > 0:
            robustness_html = "<table class='mqc-table'><thead><tr><th>Overlap Size</th><th class='num'>Mean &rho;</th><th class='num'>Std Dev</th><th class='num'>CI Width</th></tr></thead><tbody>"
            for _, row in stability_df.iterrows():
                robustness_html += f"<tr><td>{int(row['overlap_size'])}</td><td class='num'>{_safe_format(row.get('mean_rho'))}</td><td class='num'>{_safe_format(row.get('rho_std'))}</td><td class='num'>{_safe_format(row.get('ci_width'))}</td></tr>"
            robustness_html += "</tbody></table>"
        if bootstrap_plot_src:
            robustness_html += f'<p><strong>Bootstrap Correlation Distribution</strong></p><img class="mqc-img" src="{bootstrap_plot_src}" alt="Bootstrap" style="max-width:600px;">'
        jackknife_df = robustness_data.get("jackknife_table")
        if jackknife_df is not None and len(jackknife_df) > 0:
            jk_sorted = jackknife_df.copy()
            jk_sorted["abs_delta"] = jk_sorted["delta_from_full_rho"].abs()
            jk_sorted = jk_sorted.sort_values("abs_delta", ascending=False).head(10)
            robustness_html += "<table class='mqc-table'><thead><tr><th>Fusion Pair</th><th class='num'>ρ After Removal</th><th class='num'>Δ from Full</th></tr></thead><tbody>"
            for _, row in jk_sorted.iterrows():
                robustness_html += f"<tr><td>{_html_escape(str(row['fusion_pair_removed']))}</td><td class='num'>{_safe_format(row.get('rho_after_removal'))}</td><td class='num'>{_safe_format(row.get('delta_from_full_rho'))}</td></tr>"
            robustness_html += "</tbody></table>"
        effect_interp = robustness_data.get("effect_interpretation")
        if effect_interp:
            robustness_html += f"<p><strong>Biological Interpretation:</strong> {_html_escape(effect_interp[:400])}{'...' if len(effect_interp) > 400 else ''}</p>"

    sections.append(f"""
<div class="mqc-section" id="sec-robustness">
  <div class="mqc-section-header">6. Robustness &amp; Sensitivity Analysis</div>
  <div class="mqc-section-content">
    {robustness_html}
  </div>
</div>""")

    # Section 7: External Validity
    ext_html = "<p>External validity tests unavailable.</p>"
    if null_model_data:
        null_p = null_model_data.get("empirical_p_value")
        n_sim = null_model_data.get("n_simulations", 1000)
        ext_html = f"""
    <p><strong>7.1 Null Model Falsification Test</strong></p>
    <p>We randomly shuffle gene pairs {n_sim} times. Lower p-value (&lt; 0.05) = overlap unlikely due to chance.</p>
    <table class="mqc-table">
      <thead><tr><th>Metric</th><th>Value</th><th>Interpretation</th></tr></thead>
      <tbody>
        <tr><td>Observed Overlap</td><td class="num">{_safe_format(null_model_data.get("observed_overlap"))}</td><td>Actual overlap</td></tr>
        <tr><td>Null Mean Overlap</td><td class="num">{_safe_format(null_model_data.get("null_mean_overlap"))}</td><td>Random expectation</td></tr>
        <tr><td>Empirical P-Value</td><td class="num">{_safe_format_p_value(null_p)}</td><td>{"Significant" if null_p is not None and null_p < 0.05 else "Not significant"}</td></tr>
      </tbody>
    </table>"""

    if external_validity_data:
        ci_lower = external_validity_data.get("stability_ci_lower")
        ci_upper = external_validity_data.get("stability_ci_upper")
        n_boot = external_validity_data.get("bootstrap_iterations", 100)
        if ci_lower is not None and ci_upper is not None:
            ci_width = ci_upper - ci_lower
            badge_sym, badge_txt = _get_stability_badge(ci_width)
            ext_html += f"""
    <p><strong>7.2 External Validity Stability</strong></p>
    <p>Bootstrap resampling ({n_boot} iterations). Narrow CI = stable results.</p>
    <div class="mqc-status-block status-pass"><strong>STABILITY: {badge_sym} {badge_txt}</strong></div>
    <table class="mqc-table">
      <thead><tr><th>Metric</th><th>Value</th></tr></thead>
      <tbody>
        <tr><td>95% CI Lower</td><td class="num">{_safe_format_correlation(ci_lower)}</td></tr>
        <tr><td>95% CI Upper</td><td class="num">{_safe_format_correlation(ci_upper)}</td></tr>
        <tr><td>CI Width</td><td class="num">{_safe_format_correlation(ci_width)}</td></tr>
      </tbody>
    </table>"""

    if bias_analysis_data:
        mw_p = bias_analysis_data.get("mannwhitney_p_value")
        bias_detected = bias_analysis_data.get("bias_detected")
        ext_html += f"""
    <p><strong>7.3 COSMIC Sampling Bias</strong></p>
    <table class="mqc-table">
      <thead><tr><th>Metric</th><th>Value</th></tr></thead>
      <tbody>
        <tr><td>Mann-Whitney P-Value</td><td class="num">{_safe_format_p_value(mw_p)}</td></tr>
        <tr><td>Bias Detected</td><td>{"Yes" if bias_detected else "No" if bias_detected is not None else "N/A"}</td></tr>
      </tbody>
    </table>"""

    sections.append(f"""
<div class="mqc-section" id="sec-external">
  <div class="mqc-section-header">7. External Validity Assessment</div>
  <div class="mqc-section-content">
    {ext_html}
  </div>
</div>""")

    # Section 8: Claim Strength
    claim_html = "<p>Scientific claim strength classification unavailable.</p>"
    if claim_strength_data:
        classification = claim_strength_data.get("classification", "N/A")
        confidence = claim_strength_data.get("confidence_level", "N/A")
        interpretation = claim_strength_data.get("interpretation", "")
        claim_html = f"""
    <div class="mqc-highlight">
      <strong>Classification:</strong> {_html_escape(classification)}<br/>
      <strong>Confidence Level:</strong> {_html_escape(confidence)}
    </div>
    <p>{_html_escape(interpretation) if interpretation else "Scientific claim strength assessment completed."}</p>"""

    sections.append(f"""
<div class="mqc-section" id="sec-claim">
  <div class="mqc-section-header">8. Final Evidence Synthesis &amp; Claim Strength</div>
  <div class="mqc-section-content">
    {claim_html}
  </div>
</div>""")

    # Section 9: Limitations
    lim_html = limitations_text if limitations_text else "Limitations section unavailable."
    sections.append(f"""
<div class="mqc-section" id="sec-limitations">
  <div class="mqc-section-header">9. Limitations</div>
  <div class="mqc-section-content">
    <p>{_html_escape(lim_html)}</p>
  </div>
</div>""")

    # Section 10: Methods
    methods_html = ""
    if scientific_narrative:
        methods_html += f"<p>{_html_escape(scientific_narrative)}</p>"
    if confidence_statement:
        methods_html += f"<p>{_html_escape(confidence_statement)}</p>"
    if not methods_html:
        methods_html = "<p>Methods transparency section unavailable.</p>"

    sections.append(f"""
<div class="mqc-section" id="sec-methods">
  <div class="mqc-section-header">10. Methods Transparency</div>
  <div class="mqc-section-content">
    {methods_html}
  </div>
</div>""")

    # Section 11: Reproducibility
    try:
        from week2_validation import __version__ as pipeline_version
    except ImportError:
        pipeline_version = "unknown"
    numpy_ver = "N/A"
    pandas_ver = "N/A"
    scipy_ver = "N/A"
    try:
        import numpy
        numpy_ver = numpy.__version__
    except ImportError:
        pass
    try:
        import pandas
        pandas_ver = pandas.__version__
    except ImportError:
        pass
    try:
        import scipy
        scipy_ver = scipy.__version__
    except ImportError:
        pass
    python_ver = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

    sections.append(f"""
<div class="mqc-section" id="sec-repro">
  <div class="mqc-section-header">11. Reproducibility Metadata</div>
  <div class="mqc-section-content">
    <table class="mqc-table">
      <thead><tr><th>Field</th><th>Value</th></tr></thead>
      <tbody>
        <tr><td>COSMIC Source</td><td>{_html_escape(diagnostic_data.get("cosmic", {}).get("cosmic_reference_source", "N/A") if diagnostic_data else "N/A")}</td></tr>
        <tr><td>COSMIC Version</td><td>{_html_escape(cosmic_version)}</td></tr>
        <tr><td>Python Version</td><td>{python_ver}</td></tr>
        <tr><td>NumPy Version</td><td>{numpy_ver}</td></tr>
        <tr><td>Pandas Version</td><td>{pandas_ver}</td></tr>
        <tr><td>SciPy Version</td><td>{scipy_ver}</td></tr>
        <tr><td>Pipeline Version</td><td>{pipeline_version}</td></tr>
      </tbody>
    </table>
  </div>
</div>""")

    # Assemble full HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Data Integrity &amp; Statistical Validation Report</title>
  <style>{MULTIQC_CSS}</style>
</head>
<body>
  <div class="mqc-header">
    <h1>Data Integrity &amp; Statistical Validation Report</h1>
    <div class="subtitle">Data Integrity &amp; Statistical Validation</div>
  </div>
  <div class="mqc-meta">
    <strong>Run:</strong> {_html_escape(run_timestamp)} &nbsp;|&nbsp; <strong>Dataset:</strong> {_html_escape(dataset_hash)} &nbsp;|&nbsp; <strong>COSMIC:</strong> {_html_escape(cosmic_version)}
  </div>
  {"".join(sections)}
  <div class="mqc-footer">
    Generated on {datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")} UTC &nbsp;|&nbsp;
    Data Integrity &amp; Statistical Validation layer
  </div>
  <script>
    document.querySelectorAll('.mqc-section-header').forEach(function(h) {{
      h.addEventListener('click', function() {{
        this.parentElement.classList.toggle('collapsed');
      }});
    }});
  </script>
</body>
</html>"""
    return html


def generate_validation_html_report(
    results_json_path: Path,
    output_html_path: Path,
    run_metadata: dict,
    histogram_paths: Optional[Dict[str, Optional[Path]]] = None,
    interpretation_text: Optional[str] = None,
    skewness_image_path: Optional[Path] = None,
    scientific_narrative: Optional[str] = None,
    confidence_statement: Optional[str] = None,
    limitations_text: Optional[str] = None,
    robustness_data: Optional[Dict[str, Any]] = None,
    null_model_data: Optional[Dict[str, Any]] = None,
    external_validity_data: Optional[Dict[str, Any]] = None,
    bias_analysis_data: Optional[Dict[str, Any]] = None,
    claim_strength_data: Optional[Dict[str, Any]] = None,
) -> Optional[Path]:
    """
    Generate an HTML report from Data Integrity & Statistical Validation results.
    Uses the same data and structure as the PDF report.

    Args:
        results_json_path: Path to validation_status_{stem}.json (or legacy week2_status_{stem}.json).
        output_html_path: Path where HTML should be written.
        run_metadata: Runtime metadata (optional).
        histogram_paths, interpretation_text, skewness_image_path, etc.: Same as PDF generator.

    Returns:
        Path to generated HTML file, or None if generation failed.
    """
    status_data = _load_json_file(results_json_path)
    if not status_data:
        print(f"Warning: Could not load status JSON from {results_json_path}", file=sys.stderr)
        return None

    output_dir = output_html_path.parent
    filename = results_json_path.name
    # Support both validation_status_* and legacy week2_status_*
    if filename.startswith("validation_status_") and filename.endswith(".json"):
        dataset_stem = filename[18:-5]
    elif filename.startswith("week2_status_") and filename.endswith(".json"):
        dataset_stem = filename[13:-5]
    else:
        dataset_stem = status_data.get("dataset_hash", "unknown")

    from week2_validation.output_names import (
        DIAGNOSTIC_RESULTS_FILENAME_PATTERN,
        CERTIFICATION_FILENAME_PATTERN,
    )
    diagnostic_path = results_json_path.parent / DIAGNOSTIC_RESULTS_FILENAME_PATTERN.format(stem=dataset_stem)
    diagnostic_data = _load_json_file(diagnostic_path)
    # Fallback to legacy naming if new file not found
    if diagnostic_data is None:
        legacy_diag = results_json_path.parent / f"week2_diagnostic_results_{dataset_stem}.json"
        diagnostic_data = _load_json_file(legacy_diag)

    cert_path = results_json_path.parent / CERTIFICATION_FILENAME_PATTERN.format(stem=dataset_stem)
    cert_data = _load_json_file(cert_path)
    if cert_data is None:
        legacy_cert = results_json_path.parent / f"week2_dataset_certification_{dataset_stem}.json"
        cert_data = _load_json_file(legacy_cert)

    # Load advanced inference (same logic as PDF)
    advanced_inference = None
    if diagnostic_data and "advanced_inference" in diagnostic_data:
        advanced_inference = diagnostic_data["advanced_inference"]
    elif diagnostic_data and "diagnostic_results" in diagnostic_data:
        adv = diagnostic_data["diagnostic_results"].get("advanced_inference")
        if adv:
            advanced_inference = adv

    if null_model_data is None and advanced_inference:
        null_model_data = advanced_inference.get("null_model")
    if external_validity_data is None and advanced_inference:
        stab = advanced_inference.get("stability")
        if stab and isinstance(stab, dict):
            external_validity_data = {
                "stability_ci_lower": stab.get("stability_ci_lower"),
                "stability_ci_upper": stab.get("stability_ci_upper"),
                "bootstrap_iterations": stab.get("bootstrap_iterations"),
            }
    if bias_analysis_data is None and advanced_inference:
        bias_analysis_data = advanced_inference.get("bias")
    if claim_strength_data is None and advanced_inference:
        claim_strength_data = advanced_inference.get("claim_strength")

    cosmic_version = "N/A"
    if diagnostic_data:
        cosmic = diagnostic_data.get("cosmic", {})
        cosmic_version = cosmic.get("cosmic_reference_version", "N/A")

    # Resolve skewness path (check output_dir and dataset subfolder)
    skewness_path = skewness_image_path
    if not skewness_path or not skewness_path.exists():
        skewness_path = results_json_path.parent / "skewness_diagnostic.png"
    if not skewness_path or not skewness_path.exists():
        skewness_path = output_dir / "skewness_diagnostic.png"

    try:
        html_content = _build_html(
            status_data=status_data,
            diagnostic_data=diagnostic_data,
            cert_data=cert_data,
            output_dir=output_dir,
            dataset_stem=dataset_stem,
            histogram_paths=histogram_paths,
            interpretation_text=interpretation_text,
            skewness_path=skewness_path,
            scientific_narrative=scientific_narrative,
            confidence_statement=confidence_statement,
            limitations_text=limitations_text,
            robustness_data=robustness_data,
            null_model_data=null_model_data,
            external_validity_data=external_validity_data,
            bias_analysis_data=bias_analysis_data,
            claim_strength_data=claim_strength_data,
            cosmic_version=cosmic_version,
        )
        output_html_path.write_text(html_content, encoding="utf-8")
        return output_html_path
    except Exception as e:
        print(f"Warning: HTML report generation failed: {e}", file=sys.stderr)
        return None


# Backward compatibility alias
generate_week2_html_report = generate_validation_html_report
