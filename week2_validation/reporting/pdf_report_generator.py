"""
Week 2: Data Integrity & Statistical Validation — PDF Report Generator.

Generates a comprehensive PDF report summarizing validation results using ReportLab.
Professional layout: title page, TOC, headers/footers, page numbers, certification page.
"""

import hashlib
import json
import sys
import tempfile
import copy
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib.utils import ImageReader
    from reportlab.platypus import (
        CondPageBreak,
        Image,
        KeepTogether,
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False
    ImageReader = None
    CondPageBreak = None

# Spacing standardization (MultiQC-style layout)
TITLE_SPACE = 16
SECTION_SPACE = 14
SUBSECTION_SPACE = 10
PARAGRAPH_SPACE = 6
TABLE_SPACE = 10
FIGURE_SPACE = 12

# Safe minimum heights for flow-based page breaks (points)
MIN_SECTION_HEADER = 120
MIN_METRIC_TABLE = 180
MIN_FIGURE_BLOCK = 260
MIN_PAGE_FILL_PERCENT = 35  # Never leave page < 35% filled

# Legacy aliases for compatibility
SPACE_BEFORE_SECTION = SECTION_SPACE - 6
SPACE_AFTER_SECTION = SUBSECTION_SPACE
SPACE_BETWEEN_PARAGRAPHS = PARAGRAPH_SPACE
SPACE_BEFORE_TABLE = TABLE_SPACE
SPACE_AFTER_TABLE = TABLE_SPACE

# Debug: set True to print page fill percentages
DEBUG_PAGE_UTILIZATION = False


def _make_doc_template(output_path: str, add_footer_callback):
    """Create DocTemplate, with optional page utilization logging when DEBUG_PAGE_UTILIZATION."""
    if not DEBUG_PAGE_UTILIZATION:
        doc = SimpleDocTemplate(
            output_path,
            pagesize=letter,
            rightMargin=0.75 * inch,
            leftMargin=0.75 * inch,
            topMargin=0.75 * inch,
            bottomMargin=0.75 * inch,
        )
        doc.onFirstPage = add_footer_callback
        doc.onLaterPages = add_footer_callback
        return doc

    # Custom template with page utilization logger (debug only)
    from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame

    class DocTemplateWithPageLogger(BaseDocTemplate):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._page_count = 0

        def afterPage(self):
            self._page_count += 1
            try:
                frame = self.frame
                if frame and hasattr(frame, "_y") and hasattr(frame, "_y1") and hasattr(frame, "_y2"):
                    usable = frame._y2 - frame._y1
                    used = frame._y2 - frame._y if usable > 0 else 0
                    fill_pct = min(100, max(0, (used / usable) * 100)) if usable > 0 else 0
                    print(f"Page {self._page_count} fill: {fill_pct:.0f}%", file=sys.stderr)
            except Exception:
                pass

    frame = Frame(
        0.75 * inch, 0.75 * inch,
        letter[0] - 1.5 * inch, letter[1] - 1.5 * inch,
        id="normal",
    )
    doc = DocTemplateWithPageLogger(
        output_path,
        pagesize=letter,
        rightMargin=0.75 * inch,
        leftMargin=0.75 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
    )
    doc.addPageTemplates([PageTemplate(id="first", frames=[frame])])
    doc.onFirstPage = add_footer_callback
    doc.onLaterPages = add_footer_callback
    return doc


def _add_cond_break(story: list, height_pt: float) -> None:
    """Add conditional page break if less than height_pt remains (orphan protection, flow-based)."""
    if REPORTLAB_AVAILABLE and CondPageBreak is not None:
        story.append(CondPageBreak(height_pt))


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


def _safe_format_ci(value: Any) -> str:
    """Format confidence interval bounds with 3 decimal places."""
    if value is None:
        return "N/A"
    try:
        ci_val = float(value)
        if abs(ci_val) < 1e-10:
            return "0.000"
        return f"{ci_val:.3f}"
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


def _build_scientific_styles():
    """
    Create professional validation report styles.
    Typography: Times New Roman (Times-Roman) body, Arial-like (Helvetica-Bold) headings, Courier for code/hashes.
    """
    styles = getSampleStyleSheet()
    MULTIQC_ORANGE = colors.HexColor("#F18046")
    MULTIQC_DARK = colors.HexColor("#160F26")
    MULTIQC_GRAY = colors.HexColor("#666666")
    MULTIQC_LIGHT_GRAY = colors.HexColor("#f5f5f5")

    title_style = ParagraphStyle(
        "CustomTitle",
        parent=styles["Title"],
        fontSize=26,
        textColor=MULTIQC_DARK,
        spaceAfter=TITLE_SPACE,
        alignment=1,
        fontName="Helvetica-Bold",
    )
    section_style = ParagraphStyle(
        "CustomHeading1",
        parent=styles["Heading1"],
        fontSize=18,
        textColor=MULTIQC_DARK,
        spaceBefore=SECTION_SPACE,
        spaceAfter=SUBSECTION_SPACE,
        fontName="Helvetica-Bold",
    )
    subsection_style = ParagraphStyle(
        "CustomHeading2",
        parent=styles["Heading2"],
        fontSize=14,
        textColor=MULTIQC_DARK,
        spaceBefore=SUBSECTION_SPACE,
        spaceAfter=PARAGRAPH_SPACE,
        fontName="Helvetica-Bold",
    )
    body_style = ParagraphStyle(
        "CustomNormal",
        parent=styles["Normal"],
        fontSize=11,
        leading=16,
        textColor=MULTIQC_DARK,
        fontName="Times-Roman",
    )
    mono_style = ParagraphStyle(
        "CustomMono",
        parent=styles["Normal"],
        fontSize=10,
        leading=14,
        textColor=MULTIQC_DARK,
        fontName="Courier",
    )
    return {
        "title": title_style,
        "section": section_style,
        "subsection": subsection_style,
        "body": body_style,
        "mono": mono_style,
        "multiqc_orange": MULTIQC_ORANGE,
        "multiqc_dark": MULTIQC_DARK,
        "multiqc_gray": MULTIQC_GRAY,
        "multiqc_light_gray": MULTIQC_LIGHT_GRAY,
    }


def _build_highlight_box(text_lines: list[str], styles: dict) -> Table:
    """
    Create a highlight box for executive summary items.
    
    Args:
        text_lines: List of text lines to display in the box.
        styles: Style dictionary from _build_scientific_styles().
    
    Returns:
        Table element styled as a highlight box.
    """
    # Create table with single cell for highlight effect
    highlight_data = [[Paragraph("<br/>".join(text_lines), styles["body"])]]
    
    highlight_table = Table(highlight_data, colWidths=[6.5 * inch])
    # Use MultiQC light gray background
    multiqc_light_gray = colors.HexColor("#f5f5f5")
    multiqc_border = colors.HexColor("#e0e0e0")
    highlight_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), multiqc_light_gray),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 14),
        ("RIGHTPADDING", (0, 0), (-1, -1), 14),
        ("TOPPADDING", (0, 0), (-1, -1), 12),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 12),
        ("BOX", (0, 0), (-1, -1), 1, multiqc_border),
    ]))
    
    return highlight_table


def _build_status_block(title: str, status: str, confidence: str, interpretation: str, styles: dict) -> Table:
    """
    Build MultiQC-style status block with signal-first layout.
    
    Args:
        title: Block title (e.g., "DATA INTEGRITY")
        status: Status indicator (e.g., "PASS", "MODERATE")
        confidence: Confidence level (e.g., "HIGH", "LOW")
        interpretation: Short interpretation text (max 2 lines)
        styles: Style dictionary.
    
    Returns:
        Table styled as status block.
    """
    # Determine status symbol and color - MultiQC-inspired
    status_upper = str(status).upper()
    if "PASS" in status_upper or "STRONG" in status_upper or "HIGH" in status_upper:
        symbol = "✔"
        bg_color = colors.HexColor("#e8f5e9")  # Light green - success
    elif "MODERATE" in status_upper or "WARNING" in status_upper:
        symbol = "⚠"
        bg_color = colors.HexColor("#fff3e0")  # Light orange - warning (MultiQC brand color tint)
    elif "FAIL" in status_upper or "LIMITED" in status_upper or "WEAK" in status_upper or "LOW" in status_upper:
        symbol = "✗"
        bg_color = colors.HexColor("#ffebee")  # Light red/pink - error
    else:
        symbol = "○"
        bg_color = colors.HexColor("#f5f5f5")  # Light gray - neutral
    
    status_text = f"<b>{title}</b><br/>{symbol} {status}<br/><b>Confidence:</b> {confidence}<br/>{interpretation}"
    
    status_data = [[Paragraph(status_text, styles["body"])]]
    status_table = Table(status_data, colWidths=[6.5 * inch])
    status_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), bg_color),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 12),
        ("RIGHTPADDING", (0, 0), (-1, -1), 12),
        ("TOPPADDING", (0, 0), (-1, -1), 10),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
        ("BOX", (0, 0), (-1, -1), 2, colors.HexColor("#888888")),
    ]))
    
    return status_table


def _build_warning_box(message: str, styles: dict) -> Table:
    """
    Build warning box for power/system issues.
    
    Args:
        message: Warning message text.
        styles: Style dictionary.
    
    Returns:
        Table styled as warning box.
    """
    warning_text = f"<b>⚠ WARNING</b><br/>{message}"
    warning_data = [[Paragraph(warning_text, styles["body"])]]
    warning_table = Table(warning_data, colWidths=[6.5 * inch])
    warning_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#fff3cd")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 12),
        ("RIGHTPADDING", (0, 0), (-1, -1), 12),
        ("TOPPADDING", (0, 0), (-1, -1), 10),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
        ("BOX", (0, 0), (-1, -1), 2, colors.HexColor("#ffc107")),
    ]))
    
    return warning_table


def _get_stability_badge(ci_width: float) -> tuple[str, str]:
    """
    Get stability badge based on CI width.
    
    Args:
        ci_width: Confidence interval width.
    
    Returns:
        Tuple of (badge_symbol, badge_text)
    """
    if ci_width < 0.3:
        return ("🟢", "Stable")
    elif ci_width < 0.6:
        return ("🟡", "Moderate")
    else:
        return ("🔴", "Unstable")


def _style_scientific_table(table: Table, align_numbers_right: bool = True, alternating_rows: bool = True) -> Table:
    """
    Apply MultiQC-style scientific report styling to a table.
    
    Args:
        table: Table object to style.
        align_numbers_right: If True, right-align numeric columns (default: True).
        alternating_rows: If True, apply alternating row shading (default: True).
    
    Returns:
        Styled table (modified in place, also returned).
    """
    # Determine column count and row count
    num_cols = len(table._colWidths) if hasattr(table, "_colWidths") else 2
    num_rows = len(table._data) if hasattr(table, "_data") else 0
    
    # MultiQC-inspired table styling
    multiqc_header_bg = colors.HexColor("#f5f5f5")
    multiqc_dark_text = colors.HexColor("#160F26")
    multiqc_border = colors.HexColor("#e0e0e0")
    multiqc_row_bg = colors.HexColor("#fafafa")
    
    style_list = [
        # Header styling - MultiQC style: light background, dark bold text
        ("BACKGROUND", (0, 0), (-1, 0), multiqc_header_bg),
        ("TEXTCOLOR", (0, 0), (-1, 0), multiqc_dark_text),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 11),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
        ("TOPPADDING", (0, 0), (-1, 0), 8),
        ("LINEBELOW", (0, 0), (-1, 0), 2, multiqc_border),
        
        # Body styling: Times-Roman for body text, clean padding
        ("FONTNAME", (0, 1), (-1, -1), "Times-Roman"),
        ("FONTSIZE", (0, 1), (-1, -1), 10),
        ("TEXTCOLOR", (0, 1), (-1, -1), multiqc_dark_text),
        ("TOPPADDING", (0, 1), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 1), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        
        # Grid lines - MultiQC style: subtle borders
        ("GRID", (0, 0), (-1, -1), 0.5, multiqc_border),
        
        # Alignment
        ("ALIGN", (0, 0), (0, -1), "LEFT"),  # First column left-aligned (metric names)
    ]
    
    # Alternating row shading for body rows - MultiQC style
    if alternating_rows and num_rows > 1:
        for row_idx in range(1, num_rows):
            if row_idx % 2 == 0:  # Even rows (0-indexed, so row 2, 4, 6...)
                style_list.append(("BACKGROUND", (0, row_idx), (-1, row_idx), multiqc_row_bg))
    
    # Right-align numeric columns if requested
    if align_numbers_right and num_cols > 1:
        for col_idx in range(1, num_cols):
            style_list.append(("ALIGN", (col_idx, 0), (col_idx, -1), "RIGHT"))
    
    table.setStyle(TableStyle(style_list))
    return table


def _add_section_divider(story: list) -> None:
    """Add a MultiQC-inspired horizontal divider line after section headers."""
    divider_data = [[""]]
    divider_table = Table(divider_data, colWidths=[6.5 * inch])
    multiqc_border = colors.HexColor("#e0e0e0")
    divider_table.setStyle(TableStyle([
        ("LINEBELOW", (0, 0), (-1, -1), 1.5, multiqc_border),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), SPACE_AFTER_SECTION),
    ]))
    story.append(divider_table)


def _image_flowable(path: Any, max_width_inch: float = 5.5, max_height_inch: float = 4.0):
    """Create an Image flowable that preserves aspect ratio and fits within max dimensions (publication 300 DPI intent)."""
    try:
        if ImageReader is None:
            return Image(str(path), width=max_width_inch * inch, height=max_height_inch * inch)
        reader = ImageReader(str(path))
        img_w, img_h = reader.getSize()
        if not img_w or not img_h:
            return Image(str(path), width=max_width_inch * inch, height=max_height_inch * inch)
        scale_w = (max_width_inch * 72) / img_w
        scale_h = (max_height_inch * 72) / img_h
        scale = min(scale_w, scale_h, 1.0)
        w = img_w * scale
        h = img_h * scale
        return Image(str(path), width=w, height=h)
    except Exception:
        return Image(str(path), width=max_width_inch * inch, height=max_height_inch * inch)


def _build_title_page_flowables(
    dataset_name: str,
    run_date: str,
    version: str,
    overall_approval: str,
    config_summary: List[List[str]],
    styles: dict,
) -> List:
    """Build flowables for the title page: title, dataset, date, version, approval badge, config table."""
    from reportlab.lib.units import inch
    flowables = []
    flowables.append(Spacer(1, 0.5 * inch))
    flowables.append(Paragraph("Data Integrity &amp; Statistical Validation Report", styles["title"]))
    flowables.append(Spacer(1, 0.2 * inch))
    flowables.append(Paragraph(f"<b>Dataset:</b> {dataset_name}", styles["body"]))
    flowables.append(Paragraph(f"<b>Date:</b> {run_date}", styles["body"]))
    flowables.append(Paragraph(f"<b>Version:</b> {version}", styles["body"]))
    flowables.append(Spacer(1, 0.25 * inch))
    approval_upper = (overall_approval or "N/A").upper()
    if approval_upper == "APPROVED":
        badge_color = colors.HexColor("#2e7d32")
        badge_text = "APPROVED"
    elif approval_upper == "CONDITIONAL":
        badge_color = colors.HexColor("#ed6c02")
        badge_text = "CONDITIONAL"
    else:
        badge_color = colors.HexColor("#c62828")
        badge_text = approval_upper if approval_upper in ("REJECTED", "N/A") else overall_approval
    badge_data = [[Paragraph(f'<font size="16" color="white"><b>{badge_text}</b></font>', styles["body"])]]
    badge_table = Table(badge_data, colWidths=[2.5 * inch])
    badge_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), badge_color),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 12),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 12),
        ("BOX", (0, 0), (-1, -1), 2, colors.HexColor("#333333")),
    ]))
    flowables.append(badge_table)
    flowables.append(Spacer(1, 0.3 * inch))
    flowables.append(Paragraph("Pipeline configuration summary", styles["subsection"]))
    flowables.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
    if config_summary:
        config_table = Table(config_summary, colWidths=[2.5 * inch, 4 * inch])
        _style_scientific_table(config_table, align_numbers_right=False, alternating_rows=False)
        flowables.append(config_table)
    flowables.append(PageBreak())
    return flowables


def _build_toc_flowables(section_titles: List[str], styles: dict) -> List:
    """Build Table of Contents flowables (section list + page break)."""
    flowables = []
    flowables.append(Paragraph("Table of Contents", styles["section"]))
    flowables.append(Spacer(1, 0.15 * inch))
    for i, title in enumerate(section_titles, 1):
        flowables.append(Paragraph(f"{i}. {title}", styles["body"]))
        flowables.append(Spacer(1, 4))
    flowables.append(PageBreak())
    return flowables


def _build_certification_flowables(
    dataset_hash: str,
    quality_gates: List[tuple],
    overall_approval: str,
    validation_results_hash: str,
    styles: dict,
) -> List:
    """Build final Certification page: hash, quality gates, statement, signature placeholder."""
    flowables = []
    flowables.append(Spacer(1, 0.2 * inch))
    flowables.append(Paragraph("Certification", styles["section"]))
    _add_section_divider(flowables)
    flowables.append(Paragraph("<b>Dataset hash (SHA-256 prefix):</b>", styles["body"]))
    flowables.append(Paragraph(f'<font face="Courier" size="10">{dataset_hash}</font>', styles["body"]))
    flowables.append(Spacer(1, 0.15 * inch))
    flowables.append(Paragraph("<b>Quality gate results</b>", styles["subsection"]))
    if quality_gates:
        qg_data = [["Module", "Status", "Details"]]
        for mod, status, detail in quality_gates:
            qg_data.append([mod, status, (detail or "")[:60]])
        qg_table = Table(qg_data, colWidths=[1.5 * inch, 1.2 * inch, 3.8 * inch])
        _style_scientific_table(qg_table, align_numbers_right=False, alternating_rows=True)
        flowables.append(qg_table)
    flowables.append(Spacer(1, 0.2 * inch))
    approval = (overall_approval or "N/A").upper()
    if approval == "APPROVED":
        statement = "This dataset has been validated and is <b>approved</b> for downstream power-law modeling (Week 3)."
    elif approval == "CONDITIONAL":
        statement = "This dataset has been validated and is <b>conditional</b> for downstream power-law modeling (Week 3). Review quality gate details."
    else:
        statement = "This dataset has been validated and is <b>rejected</b> for downstream power-law modeling (Week 3) based on quality gate results."
    flowables.append(Paragraph(statement, styles["body"]))
    flowables.append(Spacer(1, 0.2 * inch))
    flowables.append(Paragraph("<b>Digital signature placeholder (hash of validation results):</b>", styles["body"]))
    flowables.append(Paragraph(f'<font face="Courier" size="9">{validation_results_hash}</font>', styles["body"]))
    flowables.append(Spacer(1, 0.3 * inch))
    flowables.append(Paragraph("<i>This report was generated by the Data Integrity &amp; Statistical Validation layer. The hash above binds the certification to the validation outputs.</i>", styles["body"]))
    return flowables


def _get_status_color_tint(status: str) -> Optional[colors.HexColor]:
    """
    Get light color tint for status values.
    
    Args:
        status: Status string (PASS, MODERATE, LIMITED, etc.)
    
    Returns:
        HexColor for background tint, or None if no tint should be applied.
    """
    status_upper = str(status).upper()
    if "PASS" in status_upper or "STRONG" in status_upper or "HIGH" in status_upper:
        return colors.HexColor("#e8f5e9")  # Light green - MultiQC success color
    elif "MODERATE" in status_upper or "WARNING" in status_upper:
        return colors.HexColor("#fff3e0")  # Light orange tint - MultiQC warning (matches brand orange)
    elif "LIMITED" in status_upper or "WEAK" in status_upper or "FAIL" in status_upper or "LOW" in status_upper:
        return colors.HexColor("#ffebee")  # Light red/pink - MultiQC error color
    return None


def add_histogram_section(
    story: list,
    histogram_paths: Dict[str, Optional[Path]],
    styles: dict,
    figure_number: int = 1,
) -> None:
    """
    Add histogram visualization section to PDF report with standardized figure block layout.
    
    Args:
        story: List of PDF elements (modified in place).
        histogram_paths: Dictionary with keys "linear" and/or "log" mapping to Path or None.
        styles: Style dictionary from _build_scientific_styles().
        figure_number: Starting figure number for this section (not used - PHASE 6 removes figure numbers).
    """
    # PHASE 6: Grouped Header Instead of Figure Numbers
    story.append(Paragraph("Distribution Diagnostics", styles["subsection"]))
    story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
    
    if not histogram_paths:
        story.append(Paragraph("Histogram visualizations unavailable.", styles["body"]))
        return
    # Add linear histogram if available
    linear_path = histogram_paths.get("linear")
    if linear_path and linear_path.exists():
        try:
            # PHASE 5: Standardized Figure Block (Title, Image, 3-line interpretation)
            story.append(Paragraph("<b>Protein Length Distribution Histogram</b>", styles["body"]))
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
            img = _image_flowable(linear_path, max_width_inch=5.5, max_height_inch=4)
            story.append(img)
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
            # PHASE 3: Max 3 lines interpretation
            interpretation = (
                "Histogram confirms distribution shape. Allows inspection of tail behavior, outliers, and symmetry."
            )
            story.append(Paragraph(interpretation, styles["body"]))
            story.append(Spacer(1, FIGURE_SPACE / 72.0 * inch))
        except Exception:
            story.append(Paragraph("Linear histogram unavailable.", styles["body"]))
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
    
    # Add log histogram if available
    log_path = histogram_paths.get("log")
    if log_path and log_path.exists():
        try:
            # PHASE 5: Standardized Figure Block
            story.append(Paragraph("<b>Protein Length Distribution Histogram (Log Scale)</b>", styles["body"]))
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
            img = _image_flowable(log_path, max_width_inch=5.5, max_height_inch=4)
            story.append(img)
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
            # PHASE 3: Max 3 lines interpretation
            interpretation = (
                "Log-scale enhances tail visualization across multiple orders of magnitude."
            )
            story.append(Paragraph(interpretation, styles["body"]))
            story.append(Spacer(1, FIGURE_SPACE / 72.0 * inch))
        except Exception:
            story.append(Paragraph("Log histogram unavailable.", styles["body"]))
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
    
    # If both unavailable, show message
    if not (linear_path and linear_path.exists()) and not (log_path and log_path.exists()):
        story.append(Paragraph("Histogram visualizations unavailable.", styles["body"]))


def generate_week2_pdf_report(
    results_json_path: Path,
    output_pdf_path: Path,
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
    Generate a comprehensive PDF report from Data Integrity & Statistical Validation results.

    Args:
        results_json_path: Path to validation_status_{stem}.json file (or legacy week2_status_{stem}.json).
        output_pdf_path: Path where PDF should be written.
        run_metadata: Dictionary with runtime metadata (optional).

    Returns:
        Path to generated PDF file, or None if generation failed.
    """
    if not REPORTLAB_AVAILABLE:
        print("Warning: ReportLab not available. PDF generation skipped.", file=sys.stderr)
        return None

    # Load status JSON
    status_data = _load_json_file(results_json_path)
    if not status_data:
        print(f"Warning: Could not load status JSON from {results_json_path}", file=sys.stderr)
        return None

    # Try to load diagnostic results JSON (optional)
    output_dir = results_json_path.parent
    # Extract dataset_stem from filename (validation_status_{stem}.json or legacy week2_status_{stem}.json)
    filename = results_json_path.name
    if filename.startswith("validation_status_") and filename.endswith(".json"):
        dataset_stem = filename[18:-5]
    elif filename.startswith("week2_status_") and filename.endswith(".json"):
        dataset_stem = filename[13:-5]
    else:
        dataset_stem = status_data.get("dataset_hash", "unknown")
    # Diagnostic results use dataset_stem (same as status filename)
    from week2_validation.output_names import DIAGNOSTIC_RESULTS_FILENAME_PATTERN
    diagnostic_path = output_dir / DIAGNOSTIC_RESULTS_FILENAME_PATTERN.format(stem=dataset_stem)
    diagnostic_data = _load_json_file(diagnostic_path)
    if diagnostic_data is None:
        legacy_diag = output_dir / f"week2_diagnostic_results_{dataset_stem}.json"
        diagnostic_data = _load_json_file(legacy_diag)
    
    # Load advanced inference data from diagnostic_data if not provided as parameters
    # Note: diagnostic_data structure is top-level keys (distribution, benford, cosmic, advanced_inference)
    # OR nested under "diagnostic_results" (for backward compatibility)
    advanced_inference = None
    if diagnostic_data:
        # Try top-level first (current structure)
        if "advanced_inference" in diagnostic_data:
            advanced_inference = diagnostic_data["advanced_inference"]
        # Fallback to nested structure (backward compatibility)
        elif "diagnostic_results" in diagnostic_data and "advanced_inference" in diagnostic_data["diagnostic_results"]:
            advanced_inference = diagnostic_data["diagnostic_results"]["advanced_inference"]
    
    # Override with provided parameters if available, otherwise use loaded data
    if null_model_data is None:
        if advanced_inference and advanced_inference.get("null_model"):
            null_model_data = advanced_inference["null_model"]
    
    if external_validity_data is None:
        if advanced_inference and advanced_inference.get("stability"):
            stability_data = advanced_inference["stability"]
            # Convert stability CI data to external validity format
            if stability_data and isinstance(stability_data, dict):
                external_validity_data = {
                    "stability_ci_lower": stability_data.get("stability_ci_lower"),
                    "stability_ci_upper": stability_data.get("stability_ci_upper"),
                    "bootstrap_iterations": stability_data.get("bootstrap_iterations"),
                }
    
    if bias_analysis_data is None:
        if advanced_inference and advanced_inference.get("bias"):
            bias_analysis_data = advanced_inference["bias"]
    
    # Load claim strength data from advanced_inference if not provided as parameter
    if claim_strength_data is None:
        if advanced_inference and advanced_inference.get("claim_strength"):
            claim_strength_data = advanced_inference["claim_strength"]

    # Try to load certification JSON (optional)
    from week2_validation.output_names import CERTIFICATION_FILENAME_PATTERN
    cert_path = output_dir / CERTIFICATION_FILENAME_PATTERN.format(stem=dataset_stem)
    cert_data = _load_json_file(cert_path)
    if cert_data is None:
        legacy_cert = output_dir / f"week2_dataset_certification_{dataset_stem}.json"
        cert_data = _load_json_file(legacy_cert)

    try:
        # Build scientific styles
        custom_styles = _build_scientific_styles()
        
        try:
            from week2_validation import __version__ as pipeline_version
        except ImportError:
            pipeline_version = "unknown"

        run_timestamp = status_data.get("run_timestamp_utc", "N/A")
        dataset_hash = status_data.get("dataset_hash", "N/A")
        cosmic_version = "N/A"
        if diagnostic_data:
            cosmic = diagnostic_data.get("cosmic", {})
            cosmic_version = cosmic.get("cosmic_reference_version", "N/A")

        quality_gates = status_data.get("quality_gates") or {}
        overall_approval = quality_gates.get("overall_approval", "N/A")
        gates_list = quality_gates.get("gates", [])
        config_summary = [
            ["Dataset hash", dataset_hash],
            ["Run timestamp", run_timestamp],
            ["Validation version", status_data.get("validation_version", "N/A")],
            ["Diagnostics run", ", ".join(status_data.get("diagnostics_run", [])) or "—"],
            ["Frozen required", str(status_data.get("dataset_frozen_required", "—"))],
        ]
        qg_for_cert = [(g.get("module_name", ""), g.get("status", "N/A"), (g.get("reason") or "")[:80]) for g in gates_list]
        try:
            payload = json.dumps({"status": status_data, "diagnostic": diagnostic_data or {}}, sort_keys=True, default=str)
            validation_results_hash = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:32]
        except Exception:
            validation_results_hash = hashlib.sha256(str(status_data).encode("utf-8")).hexdigest()[:32]

        title_flowables = _build_title_page_flowables(
            dataset_name=dataset_stem,
            run_date=run_timestamp,
            version=status_data.get("validation_version", "unknown"),
            overall_approval=overall_approval,
            config_summary=config_summary,
            styles=custom_styles,
        )
        toc_section_titles = [
            "Executive Summary & Quality Gates",
            "Data Integrity Validation",
            "Distribution & Statistical Testing",
            "Benford Analysis",
            "COSMIC Cross-Validation",
            "Advanced Inference",
            "Reproducibility",
            "Certification",
        ]
        toc_flowables = _build_toc_flowables(toc_section_titles, custom_styles)

        page_count_ref = [0]

        def counting_footer(canvas, doc):
            page_count_ref[0] = canvas.getPageNumber()

        def full_header_footer(canvas, doc):
            total_pages = page_count_ref[0]
            page_num = canvas.getPageNumber()
            canvas.saveState()
            multiqc_gray = colors.HexColor("#666666")
            canvas.setFont("Helvetica", 9)
            canvas.setFillColor(multiqc_gray)
            if page_num > 1:
                header_text = f"Week 2: Data Integrity & Statistical Validation | {dataset_stem}"
                canvas.drawString(0.75 * inch, letter[1] - 0.5 * inch, header_text)
            left_text = "Data Integrity & Statistical Validation"
            canvas.drawString(0.75 * inch, 0.5 * inch, left_text)
            center_text = f"COSMIC: {cosmic_version}"
            if pipeline_version != "unknown":
                center_text += f" | v{pipeline_version}"
            cw = canvas.stringWidth(center_text, "Helvetica", 9)
            canvas.drawString((letter[0] - cw) / 2, 0.5 * inch, center_text)
            page_text = f"Page {page_num} of {total_pages}"
            pw = canvas.stringWidth(page_text, "Helvetica", 9)
            canvas.drawString(letter[0] - 0.75 * inch - pw, 0.5 * inch, page_text)
            if (overall_approval or "").upper() != "APPROVED":
                canvas.setFont("Helvetica-Bold", 52)
                canvas.setFillColor(colors.HexColor("#cccccc"))
                canvas.saveState()
                canvas.translate(letter[0] / 2, letter[1] / 2)
                canvas.rotate(45)
                canvas.drawString(-2.5 * inch, -0.4 * inch, "DRAFT")
                canvas.restoreState()
            canvas.restoreState()
            if page_num == 2:
                try:
                    canvas.bookmarkPage("toc")
                    canvas.addOutlineEntry("Table of Contents", "toc", 0)
                except Exception:
                    pass
            if page_num == total_pages and total_pages > 0:
                try:
                    canvas.bookmarkPage("cert")
                    canvas.addOutlineEntry("Certification", "cert", 0)
                except Exception:
                    pass

        story = []
        story.extend(title_flowables)
        story.extend(toc_flowables)

        # =====================================================================
        # Executive Summary (first content page after TOC)
        # =====================================================================
        title_info = [
            f"<b>Run:</b> {run_timestamp} | <b>Dataset:</b> {dataset_hash} | <b>COSMIC:</b> {cosmic_version}",
        ]
        story.append(_build_highlight_box(title_info, custom_styles))
        story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))

        # =====================================================================
        # QC Dashboard Table - PRIMARY VISUAL OBJECT ON PAGE 1
        # =====================================================================
        # Extract classifications from existing computed data
        # Get exit_code first
        exit_code = status_data.get("exit_code", -1)
        data_integrity_status = "PASS" if exit_code == 0 else "FAIL"
        
        # Distribution validity - check if data exists and passes basic checks
        distribution_status = "PASS"
        if diagnostic_data:
            dist = diagnostic_data.get("distribution", {})
            if not dist or len(dist) == 0:
                distribution_status = "LIMITED"
        
        # Benford applicability
        benford_status = "N/A"
        if diagnostic_data:
            benford = diagnostic_data.get("benford", {})
            applicability = benford.get("applicability")
            if applicability is True:
                benford_status = "APPLICABLE"
            elif applicability is False:
                benford_status = "NOT APPLICABLE"
        
        # COSMIC biological agreement
        cosmic_status = "N/A"
        cosmic_confidence = "N/A"
        if diagnostic_data:
            cosmic = diagnostic_data.get("cosmic", {})
            validation_class = cosmic.get("cosmic_validation_classification", "N/A")
            if validation_class != "N/A":
                cosmic_status = validation_class
                # Extract confidence from classification
                if "STRONG" in str(validation_class).upper():
                    cosmic_confidence = "HIGH"
                elif "MODERATE" in str(validation_class).upper():
                    cosmic_confidence = "MODERATE"
                elif "LIMITED" in str(validation_class).upper() or "WEAK" in str(validation_class).upper():
                    cosmic_confidence = "LOW"
        
        # External validity
        external_validity_status = "N/A"
        if external_validity_data:
            ci_lower = external_validity_data.get("stability_ci_lower")
            ci_upper = external_validity_data.get("stability_ci_upper")
            if ci_lower is not None and ci_upper is not None:
                ci_width = ci_upper - ci_lower
                if ci_width < 0.3:
                    external_validity_status = "STABLE"
                elif ci_width < 0.6:
                    external_validity_status = "MODERATE"
                else:
                    external_validity_status = "LIMITED"
        
        # Overall scientific claim strength
        overall_status = "N/A"
        overall_confidence = "N/A"
        if claim_strength_data:
            classification = claim_strength_data.get("classification", "N/A")
            confidence_level = claim_strength_data.get("confidence_level", "N/A")
            if classification != "N/A":
                overall_status = classification
                overall_confidence = confidence_level if confidence_level != "N/A" else "N/A"
        
        # Build executive dashboard table
        dashboard_data = [
            ["QC Layer", "Status", "Confidence"],
            ["Data Integrity", data_integrity_status, "HIGH" if data_integrity_status == "PASS" else "LOW"],
            ["Distribution Validity", distribution_status, "HIGH" if distribution_status == "PASS" else "MODERATE"],
            ["Benford Applicability", benford_status, "N/A"],
            ["COSMIC Biological Agreement", cosmic_status, cosmic_confidence],
            ["External Validity", external_validity_status, "HIGH" if external_validity_status == "STABLE" else "MODERATE" if external_validity_status == "MODERATE" else "LOW"],
            ["Overall Scientific Claim Strength", overall_status, overall_confidence],
        ]
        
        # Convert table data to Paragraphs for proper text wrapping
        dashboard_data_para = []
        for row_idx, row in enumerate(dashboard_data):
            para_row = []
            for cell in row:
                # Use bold style for header row
                if row_idx == 0:
                    para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                else:
                    para_row.append(Paragraph(str(cell), custom_styles["body"]))
            dashboard_data_para.append(para_row)
        
        # Adjust column widths: wider Status column to prevent truncation
        dashboard_table = Table(dashboard_data_para, colWidths=[2.5 * inch, 2.5 * inch, 1.5 * inch])
        
        # Build style list with color tints
        dashboard_style_list = [
            # Header styling
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8e8e8")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#1a1a1a")),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 11),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
            ("TOPPADDING", (0, 0), (-1, 0), 8),
            ("LINEBELOW", (0, 0), (-1, 0), 2, colors.HexColor("#888888")),
            
            # Body styling
            ("FONTSIZE", (0, 1), (-1, -1), 10),
            ("TOPPADDING", (0, 1), (-1, -1), 6),
            ("BOTTOMPADDING", (0, 1), (-1, -1), 6),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cccccc")),
            ("ALIGN", (0, 0), (0, -1), "LEFT"),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),  # Top align for wrapped text
        ]
        
        # Add alternating row shading and color tints
        for row_idx in range(1, len(dashboard_data_para)):
            if row_idx % 2 == 0:
                dashboard_style_list.append(("BACKGROUND", (0, row_idx), (-1, row_idx), colors.HexColor("#f8f8f8")))
            
            # Add color tint to status column (column 1)
            status_cell_value = dashboard_data[row_idx][1]
            tint_color = _get_status_color_tint(status_cell_value)
            if tint_color:
                # Override background for status cell
                dashboard_style_list.append(("BACKGROUND", (1, row_idx), (1, row_idx), tint_color))
        
        dashboard_table.setStyle(TableStyle(dashboard_style_list))
        
        # Dashboard is PRIMARY VISUAL OBJECT - no section header needed
        story.append(dashboard_table)
        story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))

        # =====================================================================
        # Data Provenance: step → status → details table
        # =====================================================================
        provenance_chain = (diagnostic_data or {}).get("provenance") or []
        if provenance_chain:
            _add_cond_break(story, MIN_METRIC_TABLE)
            story.append(Paragraph("Data Provenance (Chain of Custody)", custom_styles["section"]))
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
            prov_table_data = [["Step", "Status", "Details"]]
            for s in provenance_chain:
                raw_details = str(s.get("details", "") or "")
                details = raw_details[:80] + ("..." if len(raw_details) > 80 else "")
                prov_table_data.append([
                    str(s.get("step", "")),
                    str(s.get("status", "N/A")),
                    details,
                ])
            prov_table = Table(prov_table_data, colWidths=[1.8 * inch, 1.2 * inch, 3.0 * inch])
            _style_scientific_table(prov_table, align_numbers_right=False, alternating_rows=True)
            story.append(prov_table)
            story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))
        
        # =====================================================================
        # PHASE 2: Status Blocks Replace Executive Paragraphs
        # =====================================================================
        # Extract data for status blocks
        status = status_data.get("status", "UNKNOWN")
        dq = status_data.get("data_quality", {})
        total_rows = dq.get("total_rows_original", "N/A")
        rows_used = dq.get("rows_used_for_analysis", "N/A")
        
        approved = False
        if cert_data:
            approved = cert_data.get("approved_for_powerlaw_modeling", False)
        pass_fail = "PASS" if approved else "Not Approved"
        
        # Build status blocks (signal-first layout)
        story.append(_build_status_block(
            "DATA INTEGRITY",
            data_integrity_status,
            "HIGH" if data_integrity_status == "PASS" else "LOW",
            f"Schema validation: {status}. Rows: {total_rows} total, {rows_used} used.",
            custom_styles
        ))
        story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        
        # Distribution status block
        dist_interpretation = "Heavy-tailed biologically plausible distribution."
        if diagnostic_data:
            dist = diagnostic_data.get("distribution", {})
            skewness = dist.get("skewness")
            if skewness is not None:
                dist_interpretation = f"Skewness {_safe_format(skewness, 3)} indicates biologically plausible distribution."
        
        story.append(_build_status_block(
            "DISTRIBUTION VALIDITY",
            distribution_status,
            "HIGH" if distribution_status == "PASS" else "MODERATE",
            dist_interpretation,
            custom_styles
        ))
        story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        
        # COSMIC status block
        cosmic_interpretation = "Biological agreement assessment."
        if diagnostic_data:
            cosmic = diagnostic_data.get("cosmic", {})
            overlap_count = cosmic.get("overlap_count", 0)
            if isinstance(overlap_count, (int, float)) and overlap_count < 10:
                cosmic_interpretation = f"Low overlap (n={int(overlap_count)}) limits statistical power."
        
        story.append(_build_status_block(
            "COSMIC BIOLOGICAL AGREEMENT",
            cosmic_status if cosmic_status != "N/A" else "N/A",
            cosmic_confidence if cosmic_confidence != "N/A" else "N/A",
            cosmic_interpretation,
            custom_styles
        ))
        
        # =====================================================================
        # PHASE 10: Interpretation Rule Micro-Box
        # =====================================================================
        interpretation_rules = [
            "<b>Interpretation Rules</b><br/>"
            "• Effect size evaluated before p-value<br/>"
            "• Low overlap reduces statistical power<br/>"
            "• Wide CI = unstable biological inference"
        ]
        story.append(_build_highlight_box(interpretation_rules, custom_styles))
        story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))

        # =====================================================================
        # Quality Gates (per-module PASS/WARN/FAIL and overall approval)
        # =====================================================================
        quality_gates = status_data.get("quality_gates") or {}
        gates_list = quality_gates.get("gates", [])
        overall_approval = quality_gates.get("overall_approval", "N/A")
        if gates_list:
            story.append(Paragraph("Quality Gates", custom_styles["section"]))
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
            qg_table_data = [["Module", "Status", "Metric", "Value", "Threshold", "Reason"]]
            for g in gates_list:
                status = str(g.get("status", ""))
                val = g.get("metric_value")
                thresh = g.get("threshold")
                val_str = _safe_format(val, 4) if val is not None else "—"
                thresh_str = _safe_format(thresh, 4) if thresh is not None else "—"
                reason = (g.get("reason") or "")[:60] + ("..." if len((g.get("reason") or "")) > 60 else "")
                qg_table_data.append([
                    str(g.get("module_name", "")),
                    status,
                    str(g.get("metric_name", "")),
                    val_str,
                    thresh_str,
                    reason,
                ])
            qg_para = [[Paragraph(str(cell), custom_styles["body"]) for cell in row] for row in qg_table_data]
            qg_table = Table(qg_para, colWidths=[1.0 * inch, 0.9 * inch, 1.0 * inch, 0.6 * inch, 0.6 * inch, 1.8 * inch])
            _style_scientific_table(qg_table, align_numbers_right=True, alternating_rows=True)
            story.append(qg_table)
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
            story.append(Paragraph(
                f"<b>Overall approval:</b> {overall_approval}. "
                "APPROVED = no FAIL, ≥2 PASS; CONDITIONAL = 1+ WARN, no FAIL; REJECTED = any FAIL.",
                custom_styles["body"]
            ))
            story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))

        # =====================================================================
        # Section 2: Data Integrity Validation
        # =====================================================================
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))
        story.append(Paragraph("2. Data Integrity Validation", custom_styles["section"]))
        _add_section_divider(story)
        
        # PHASE 2: Signal Block FIRST
        schema_status = "PASSED" if exit_code == 0 else "FAILED"
        excluded_fraction = dq.get("excluded_fraction", 0.0)
        warning_flag = dq.get("warning_flag", False)
        high_risk_flag = dq.get("high_risk_flag", False)
        
        integrity_status = "PASS" if exit_code == 0 else "FAIL"
        integrity_confidence = "HIGH" if exit_code == 0 and not warning_flag and not high_risk_flag else "LOW"
        integrity_interp = f"Schema: {schema_status}. Excluded: {_safe_format(excluded_fraction * 100, 2)}%."
        if warning_flag:
            integrity_interp += " Warning flags present."
        if high_risk_flag:
            integrity_interp += " High risk detected."
        
        story.append(_build_status_block(
            "DATA INTEGRITY STATUS",
            integrity_status,
            integrity_confidence,
            integrity_interp,
            custom_styles
        ))
        story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
        
        # PHASE 4: Table-First Evidence Display
        _add_cond_break(story, MIN_METRIC_TABLE)
        integrity_table_data = [
            ["Metric", "Value", "Interpretation"],
            ["Schema Validation", schema_status, "PASS" if exit_code == 0 else "FAIL"],
            ["Excluded Fraction", f"{_safe_format(excluded_fraction * 100, 2)}%", "Low exclusion" if excluded_fraction < 0.1 else "Moderate exclusion" if excluded_fraction < 0.3 else "High exclusion"],
            ["Warning Flag", "True" if warning_flag else "False", "Data quality concerns" if warning_flag else "No warnings"],
            ["High Risk Flag", "True" if high_risk_flag else "False", "Critical issues detected" if high_risk_flag else "No critical issues"],
        ]
        
        integrity_table_para = []
        for row_idx, row in enumerate(integrity_table_data):
            para_row = []
            for cell in row:
                if row_idx == 0:
                    para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                else:
                    para_row.append(Paragraph(str(cell), custom_styles["body"]))
            integrity_table_para.append(para_row)
        
        integrity_table = Table(integrity_table_para, colWidths=[2 * inch, 1.5 * inch, 3 * inch])
        _style_scientific_table(integrity_table, align_numbers_right=False, alternating_rows=True)
        story.append(integrity_table)
        story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))

        # =====================================================================
        # Section 3: Distribution & Statistical Testing
        # =====================================================================
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))
        story.append(Paragraph("3. Distribution & Statistical Testing", custom_styles["section"]))
        _add_section_divider(story)
        
        # PHASE 2: Signal Block FIRST
        dist_status_signal = "PASS"
        dist_confidence_signal = "HIGH"
        if diagnostic_data:
            dist = diagnostic_data.get("distribution", {})
            if not dist or len(dist) == 0:
                dist_status_signal = "LIMITED"
                dist_confidence_signal = "LOW"
        
        story.append(_build_status_block(
            "DISTRIBUTION VALIDITY STATUS",
            dist_status_signal,
            dist_confidence_signal,
            "Tests whether data follows log-normal distribution (common in biology). Smaller test statistics = better fit to expected pattern.",
            custom_styles
        ))
        story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))

        if diagnostic_data:
            dist = diagnostic_data.get("distribution", {})
            log_norm = diagnostic_data.get("log_normality", {})

            # PHASE 4: Table-First Evidence Display
            _add_cond_break(story, MIN_METRIC_TABLE)
            stats_data = [
                ["Metric", "Value", "Interpretation"],
                ["Mean", _safe_format(dist.get("mean")), "Central tendency"],
                ["Median", _safe_format(dist.get("median")), "Robust central value"],
                ["Std Dev", _safe_format(dist.get("std")), "Dispersion measure"],
                ["Skewness", _safe_format(dist.get("skewness")), "Distribution asymmetry"],
                ["Minimum", _safe_format(dist.get("minimum")), "Lower bound"],
                ["Maximum", _safe_format(dist.get("maximum")), "Upper bound"],
            ]

            stats_table_para = []
            for row_idx, row in enumerate(stats_data):
                para_row = []
                for cell in row:
                    if row_idx == 0:
                        para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                    else:
                        para_row.append(Paragraph(str(cell), custom_styles["body"]))
                stats_table_para.append(para_row)

            stats_table = Table(stats_table_para, colWidths=[2 * inch, 1.5 * inch, 3 * inch])
            _style_scientific_table(stats_table, align_numbers_right=True, alternating_rows=True)
            story.append(stats_table)
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))

            # PHASE 4: Table-First for Test Results
            _add_cond_break(story, MIN_METRIC_TABLE)
            ks_stat = log_norm.get("ks_statistic")
            ad_stat = log_norm.get("ad_statistic")
            scipy_available = log_norm.get("scipy_available", False)
            
            test_results_data = [
                ["Test", "Statistic", "Interpretation"],
                ["KS Test", _safe_format(ks_stat, 6) if ks_stat is not None else "N/A", "Smaller = better log-normal fit"],
                ["AD Test", _safe_format(ad_stat, 6) if ad_stat is not None else "N/A", "Weighted tail assessment"],
                ["SciPy Available", "Yes" if scipy_available else "No", "Advanced statistics enabled" if scipy_available else "Basic stats only"],
            ]
            
            test_table_para = []
            for row_idx, row in enumerate(test_results_data):
                para_row = []
                for cell in row:
                    if row_idx == 0:
                        para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                    else:
                        para_row.append(Paragraph(str(cell), custom_styles["body"]))
                test_table_para.append(para_row)
            
            test_table = Table(test_table_para, colWidths=[2 * inch, 1.5 * inch, 3 * inch])
            _style_scientific_table(test_table, align_numbers_right=True, alternating_rows=True)
            story.append(test_table)
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
        else:
            story.append(Paragraph("No statistical validation data available.", custom_styles["body"]))

        # PHASE 3: Compress Interpretation (max 4 lines)
        if interpretation_text:
            # Split long text into bullets if needed
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            story.append(Paragraph("<b>Interpretation:</b>", custom_styles["subsection"]))
            # Limit to 4 lines - truncate if longer
            interp_lines = interpretation_text.split('\n')[:4]
            story.append(Paragraph("<br/>".join(interp_lines), custom_styles["body"]))
            story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        
        # PHASE 6: Grouped Header Instead of Figure Numbers
        story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
        story.append(Paragraph("Distribution Diagnostics", custom_styles["subsection"]))
        story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        
        # Find skewness_diagnostic.png - always check output directory first (most reliable)
        final_skewness_path = None
        output_dir = results_json_path.parent if results_json_path else None
        
        # Check output directory first (where file is generated)
        if output_dir:
            skewness_diagnostic_path = output_dir / "skewness_diagnostic.png"
            if skewness_diagnostic_path.exists():
                final_skewness_path = skewness_diagnostic_path
        
        # If not found in output directory, try the provided path
        if not final_skewness_path and skewness_image_path:
            try:
                if not isinstance(skewness_image_path, Path):
                    skewness_image_path = Path(skewness_image_path)
                if skewness_image_path.exists():
                    final_skewness_path = skewness_image_path
            except Exception:
                pass
        
        # Add image if we found a valid path (image page: allow isolation, do not compress)
        if final_skewness_path and final_skewness_path.exists():
            try:
                _add_cond_break(story, MIN_FIGURE_BLOCK)  # Ensure room for figure block
                # PHASE 5: Standardized Figure Block (Title, Image, 3-line interpretation)
                skewness_val = None
                if diagnostic_data:
                    dist = diagnostic_data.get("distribution", {})
                    skewness_val = dist.get("skewness")
                
                if skewness_val is not None and isinstance(skewness_val, (int, float)):
                    skew_val = float(skewness_val)
                    story.append(Paragraph("<b>Skewness Diagnostic</b>", custom_styles["body"]))
                    story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
                    # Add the image (compact size to fit page)
                    img_path_str = str(final_skewness_path.resolve())
                    img = _image_flowable(img_path_str, max_width_inch=6.5, max_height_inch=4.25)
                    story.append(img)
                    story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
                    # PHASE 3: Max 3 lines interpretation
                    interpretation_text_skew = (
                        f"Skewness {skew_val:.4f} indicates {'symmetrical' if abs(skew_val) < 0.5 else 'positive skew' if skew_val > 0 else 'negative skew'} distribution. "
                        f"Gauge confirms structural characteristics."
                    )
                    story.append(Paragraph(interpretation_text_skew, custom_styles["body"]))
                    story.append(Spacer(1, FIGURE_SPACE / 72.0 * inch))
                else:
                    story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
            except Exception as e:
                # Log error but don't fail PDF generation
                print(f"Warning: Could not add skewness image to PDF: {e}", file=sys.stderr)
                story.append(Paragraph("Skewness visualization unavailable.", custom_styles["body"]))
                story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
        else:
            # Image not found - show message
            story.append(Paragraph("Skewness visualization unavailable.", custom_styles["body"]))
            story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))

        # Add histogram visualization section if paths provided (image page: allow isolation)
        if histogram_paths:
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            _add_cond_break(story, MIN_FIGURE_BLOCK)
            add_histogram_section(story, histogram_paths, custom_styles, figure_number=1)

        # =====================================================================
        # Section 4: Benford Analysis
        # =====================================================================
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))
        story.append(Paragraph("4. Benford Analysis", custom_styles["section"]))
        _add_section_divider(story)
        
        # Clear explanation for both statistical and non-statistical audiences
        story.append(Paragraph(
            "<b>What is Benford's Law?</b> In many real-world datasets, numbers starting with 1 appear more frequently (about 30%) than numbers starting with 9 (about 5%). This pattern helps detect data anomalies. <b>Note:</b> This test only applies to datasets spanning multiple orders of magnitude (e.g., values from 1 to 1000+).",
            custom_styles["body"]
        ))
        story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))

        if diagnostic_data:
            benford = diagnostic_data.get("benford", {})
            if benford:
                # Applicability summary
                applicability = benford.get("applicability")
                reason = benford.get("reason_if_not_applicable", "")
                scale_span = benford.get("scale_span_orders_of_magnitude")

                benford_summary_lines = [f"<b>Benford Applicable:</b> {applicability}"]
                if reason:
                    benford_summary_lines.append(f"<b>Reason:</b> {reason}")
                if scale_span is not None:
                    benford_summary_lines.append(f"<b>Scale Span:</b> {_safe_format(scale_span)} orders of magnitude")
                chi_sq = benford.get("chi_squared_statistic")
                df_val = benford.get("degrees_of_freedom")
                p_val = benford.get("p_value")
                if chi_sq is not None:
                    benford_summary_lines.append(f"<b>Chi-squared statistic:</b> {_safe_format(chi_sq, 4)}")
                if df_val is not None:
                    benford_summary_lines.append(f"<b>Degrees of freedom:</b> {df_val}")
                if p_val is not None:
                    benford_summary_lines.append(f"<b>p-value:</b> {_safe_format_p_value(p_val)}")
                if p_val is not None and isinstance(p_val, (int, float)):
                    try:
                        if float(p_val) > 0.05:
                            benford_summary_lines.append("<i>Diagnostic note (descriptive only):</i> p &gt; 0.05 indicates consistency with Benford distribution.")
                        else:
                            benford_summary_lines.append("<i>Diagnostic note (descriptive only):</i> p &lt; 0.05 indicates deviation from Benford distribution.")
                    except (TypeError, ValueError):
                        pass
                
                story.append(_build_highlight_box(benford_summary_lines, custom_styles))
                story.append(Spacer(1, 0.08 * inch))

                # Observed vs Expected frequencies table
                obs_freq = benford.get("observed_frequencies", {})
                exp_freq = benford.get("expected_frequencies", {})

                if obs_freq or exp_freq:
                    _add_cond_break(story, MIN_METRIC_TABLE)
                    freq_data = [["Digit", "Observed", "Expected"]]
                    for digit in range(1, 10):
                        obs = obs_freq.get(digit, obs_freq.get(str(digit), 0))
                        exp = exp_freq.get(digit, exp_freq.get(str(digit), 0))
                        freq_data.append([
                            str(digit),
                            _safe_format(obs, 4),
                            _safe_format(exp, 4),
                        ])

                    freq_table = Table(freq_data, colWidths=[1.5 * inch, 2.5 * inch, 2.5 * inch])
                    _style_scientific_table(freq_table, align_numbers_right=True, alternating_rows=True)
                    story.append(Spacer(1, SPACE_BEFORE_TABLE / 72.0 * inch))
                    story.append(freq_table)
                    story.append(Spacer(1, SPACE_AFTER_TABLE / 72.0 * inch))
                
                # Add Benford visualization if available
                # Try to find benford_analysis.png in the same directory as results_json_path
                results_dir = results_json_path.parent if results_json_path else None
                benford_img_path = None
                if results_dir:
                    benford_img_path = results_dir / "benford_analysis.png"
                
                if benford_img_path and benford_img_path.exists():
                    try:
                        _add_cond_break(story, MIN_FIGURE_BLOCK)  # Image page: allow isolation
                        # PHASE 5: Standardized Figure Block (Title, Image, 3-line interpretation)
                        story.append(Paragraph("<b>Benford Analysis Visualization</b>", custom_styles["body"]))
                        story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
                        # Use absolute path to ensure ReportLab can find it
                        benford_img_path_str = str(benford_img_path.resolve())
                        img = _image_flowable(benford_img_path_str, max_width_inch=5.5, max_height_inch=3.5)
                        story.append(img)
                        story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
                        # PHASE 3: Max 3 lines interpretation
                        interpretation_text_benford = (
                            "Chart shows deviation between observed (blue) and expected (red) first digits. "
                            "Visual mismatch confirms statistical results. Dataset lacks scale span for reliable Benford analysis."
                        )
                        story.append(Paragraph(interpretation_text_benford, custom_styles["body"]))
                        story.append(Spacer(1, FIGURE_SPACE / 72.0 * inch))
                    except Exception as e:
                        # Silently skip if image can't be loaded
                        pass
            else:
                story.append(Paragraph("No Benford analysis data available.", custom_styles["body"]))
        else:
            story.append(Paragraph("No Benford analysis data available.", custom_styles["body"]))

        # =====================================================================
        # Section 5: COSMIC Cross-Validation
        # =====================================================================
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))
        story.append(Paragraph("5. COSMIC Cross-Validation", custom_styles["section"]))
        _add_section_divider(story)
        
        if diagnostic_data:
            cosmic = diagnostic_data.get("cosmic", {})
            if cosmic:
                # PHASE 2: Signal Block FIRST
                cosmic_signal_status = "N/A"
                cosmic_signal_power = "N/A"
                cosmic_signal_uncertainty = "N/A"
                
                validation_class = cosmic.get("cosmic_validation_classification", "N/A")
                overlap_count = cosmic.get("overlap_count", 0)
                
                if validation_class != "N/A":
                    cosmic_signal_status = validation_class
                    if isinstance(overlap_count, (int, float)):
                        if overlap_count < 10:
                            cosmic_signal_power = "LOW"
                            cosmic_signal_uncertainty = "HIGH"
                        elif overlap_count < 30:
                            cosmic_signal_power = "MODERATE"
                            cosmic_signal_uncertainty = "MODERATE"
                        else:
                            cosmic_signal_power = "HIGH"
                            cosmic_signal_uncertainty = "LOW"
                
                signal_text = f"<b>COSMIC CROSS VALIDATION</b><br/>STATUS: {cosmic_signal_status}<br/>POWER: {cosmic_signal_power}<br/>UNCERTAINTY: {cosmic_signal_uncertainty}"
                signal_data = [[Paragraph(signal_text, custom_styles["body"])]]
                signal_table = Table(signal_data, colWidths=[6.5 * inch])
                signal_table.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#fff9c4") if cosmic_signal_power == "LOW" else colors.HexColor("#e8f5e9")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 12),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 12),
                    ("TOPPADDING", (0, 0), (-1, -1), 10),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
                    ("BOX", (0, 0), (-1, -1), 2, colors.HexColor("#888888")),
                ]))
                story.append(signal_table)
                story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
                
                # PHASE 8: Power Warning System
                if isinstance(overlap_count, (int, float)) and overlap_count < 10:
                    story.append(_build_warning_box(
                        f"COSMIC overlap < 10 fusion pairs (n={int(overlap_count)}). Statistical uncertainty expected.",
                        custom_styles
                    ))
                    story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
                
                # PHASE 4: Table-First - Overlap Statistics
                _add_cond_break(story, MIN_METRIC_TABLE)
                overlap_table_data = [
                    ["Statistic", "Value", "Interpretation"],
                    ["Total Fusions (Ours)", str(cosmic.get("total_fusions_ours", "N/A")), "Dataset size"],
                    ["Total Fusions (COSMIC)", str(cosmic.get("total_fusions_cosmic", "N/A")), "Reference database size"],
                    ["Overlap Count", str(cosmic.get("overlap_count", "N/A")), "Common fusion pairs"],
                    ["Only in Ours", str(cosmic.get("only_in_ours_count", "N/A")), "Dataset-specific"],
                    ["Only in COSMIC", str(cosmic.get("only_in_cosmic_count", "N/A")), "Reference-specific"],
                ]
                
                overlap_table_para = []
                for row_idx, row in enumerate(overlap_table_data):
                    para_row = []
                    for cell in row:
                        if row_idx == 0:
                            para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                        else:
                            para_row.append(Paragraph(str(cell), custom_styles["body"]))
                    overlap_table_para.append(para_row)
                
                overlap_table = Table(overlap_table_para, colWidths=[2 * inch, 1.5 * inch, 3 * inch])
                _style_scientific_table(overlap_table, align_numbers_right=True, alternating_rows=True)
                story.append(overlap_table)
                story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
                
                # Add brief explanation before table for clarity
                story.append(Paragraph(
                    "<b>Spearman Correlation:</b> Measures how well the rank order of fusion recurrence matches between our data and COSMIC. Values range from -1 (opposite ranking) to +1 (perfect matching). A value near 0 indicates no relationship.",
                    custom_styles["body"]
                ))
                story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
                
                # PHASE 4: Table-First - Spearman correlation
                _add_cond_break(story, MIN_METRIC_TABLE)
                spearman_rho = cosmic.get("spearman_rho")
                spearman_p = cosmic.get("spearman_p_value")
                enrichment_p = cosmic.get("enrichment_p_value")
                overlap_count = cosmic.get("overlap_count", 0)
                
                # Determine interpretation
                rho_interp = "N/A"
                if spearman_rho is not None and isinstance(spearman_rho, (int, float)):
                    abs_rho = abs(float(spearman_rho))
                    if abs_rho >= 0.7:
                        rho_interp = "Strong rank agreement"
                    elif abs_rho >= 0.4:
                        rho_interp = "Moderate rank agreement"
                    else:
                        rho_interp = "Weak rank agreement"
                    if isinstance(overlap_count, (int, float)) and overlap_count < 10:
                        rho_interp += " (low power)"
                
                spearman_data = [
                    ["Metric", "Value", "Interpretation"],
                    ["Spearman Rho", _safe_format_correlation(spearman_rho), rho_interp],
                    ["Spearman P Value", _safe_format_p_value(spearman_p), "Significance test"],
                    ["Enrichment P Value", _safe_format_p_value(enrichment_p), "Top fusion overlap test"],
                ]

                spearman_table_para = []
                for row_idx, row in enumerate(spearman_data):
                    para_row = []
                    for cell in row:
                        if row_idx == 0:
                            para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                        else:
                            para_row.append(Paragraph(str(cell), custom_styles["body"]))
                    spearman_table_para.append(para_row)

                spearman_table = Table(spearman_table_para, colWidths=[2 * inch, 1.5 * inch, 3 * inch])
                _style_scientific_table(spearman_table, align_numbers_right=True, alternating_rows=True)
                story.append(spearman_table)
                story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))

                # Quality gate score - status block format
                validation_score = cosmic.get("cosmic_validation_score")
                validation_class = cosmic.get("cosmic_validation_classification", "N/A")
                
                if validation_class != "N/A":
                    story.append(_build_status_block(
                        "COSMIC VALIDATION CLASSIFICATION",
                        validation_class,
                        cosmic_confidence if cosmic_confidence != "N/A" else "N/A",
                        f"Quality gate score: {_safe_format(validation_score)}",
                        custom_styles
                    ))
                    story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            else:
                story.append(Paragraph("No COSMIC validation data available.", custom_styles["body"]))
        else:
            story.append(Paragraph("No COSMIC validation data available.", custom_styles["body"]))

        # =====================================================================
        # Section 6: Robustness & Sensitivity Analysis
        # =====================================================================
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SECTION_SPACE / 72.0 * inch))
        story.append(Paragraph("6. Robustness & Sensitivity Analysis", custom_styles["section"]))
        _add_section_divider(story)

        if robustness_data:
            # PHASE 2: Signal Block FIRST
            robustness_status = "AVAILABLE"
            robustness_confidence = "MODERATE"
            robustness_interp = "Bootstrap and jackknife analyses assess correlation stability."
            story.append(_build_status_block(
                "ROBUSTNESS ANALYSIS STATUS",
                robustness_status,
                robustness_confidence,
                robustness_interp,
                custom_styles
            ))
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            # PHASE 4: Table-First - Stability vs Overlap Table
            stability_df = robustness_data.get("stability_table")
            if stability_df is not None and len(stability_df) > 0:
                _add_cond_break(story, MIN_METRIC_TABLE)
                stability_data_table = [["Overlap Size", "Mean ρ", "Std Dev", "CI Width"]]
                # Support both DataFrame (in-memory) and list-of-dicts (e.g. from JSON)
                if hasattr(stability_df, "iterrows"):
                    rows_iter = (row for _, row in stability_df.iterrows())
                else:
                    rows_iter = iter(stability_df)
                for row in rows_iter:
                    overlap_val = row.get("overlap_size") if hasattr(row, "get") else (row["overlap_size"] if isinstance(row, dict) else None)
                    if overlap_val is None:
                        continue
                    stability_data_table.append([
                        str(int(overlap_val)),
                        _safe_format(row.get("mean_rho") if hasattr(row, "get") else None),
                        _safe_format(row.get("rho_std") if hasattr(row, "get") else None),
                        _safe_format(row.get("ci_width") if hasattr(row, "get") else None),
                    ])
                
                stability_table_para = []
                for row_idx, row in enumerate(stability_data_table):
                    para_row = []
                    for cell in row:
                        if row_idx == 0:
                            para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                        else:
                            para_row.append(Paragraph(str(cell), custom_styles["body"]))
                    stability_table_para.append(para_row)
                
                stability_table = Table(stability_table_para, colWidths=[1.5 * inch, 1.5 * inch, 1.5 * inch, 2 * inch])
                _style_scientific_table(stability_table, align_numbers_right=True, alternating_rows=True)
                story.append(stability_table)
                story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            
            # PHASE 5: Standardized Figure Block - Bootstrap Distribution (image page: allow isolation)
            bootstrap_plot_path = robustness_data.get("bootstrap_plot_path")
            if bootstrap_plot_path and Path(bootstrap_plot_path).exists():
                try:
                    _add_cond_break(story, MIN_FIGURE_BLOCK)
                    story.append(Paragraph("<b>Bootstrap Correlation Distribution</b>", custom_styles["body"]))
                    story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
                    img = _image_flowable(bootstrap_plot_path, max_width_inch=5.5, max_height_inch=3.5)
                    story.append(img)
                    story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
                    # PHASE 3: Max 3 lines interpretation
                    story.append(Paragraph("Bootstrap resampling shows correlation distribution. Wider distribution indicates higher uncertainty.", custom_styles["body"]))
                    story.append(Spacer(1, FIGURE_SPACE / 72.0 * inch))
                except Exception:
                    story.append(Paragraph("Bootstrap distribution plot unavailable.", custom_styles["body"]))
                    story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
            
            # PHASE 4: Table-First - Jackknife Sensitivity Summary
            jackknife_df = robustness_data.get("jackknife_table")
            if jackknife_df is not None and len(jackknife_df) > 0:
                _add_cond_break(story, MIN_METRIC_TABLE)
                # Show top 10 most sensitive pairs (by absolute delta). Support DataFrame or list-of-dicts.
                if hasattr(jackknife_df, "copy"):
                    jackknife_sorted = jackknife_df.copy()
                    jackknife_sorted["abs_delta"] = jackknife_sorted["delta_from_full_rho"].abs()
                    jackknife_rows = list(jackknife_sorted.sort_values("abs_delta", ascending=False).head(10).iterrows())
                    jackknife_rows = [r for _, r in jackknife_rows]
                else:
                    def _abs_delta(r):
                        v = r.get("delta_from_full_rho") if hasattr(r, "get") else None
                        try:
                            return abs(float(v)) if v is not None else 0.0
                        except (TypeError, ValueError):
                            return 0.0
                    jackknife_rows = sorted(jackknife_df, key=_abs_delta, reverse=True)[:10]
                
                jackknife_data_table = [["Fusion Pair", "ρ After Removal", "Δ from Full"]]
                for row in jackknife_rows:
                    pair = row.get("fusion_pair_removed", "—") if hasattr(row, "get") else "—"
                    jackknife_data_table.append([
                        str(pair),
                        _safe_format(row.get("rho_after_removal") if hasattr(row, "get") else None),
                        _safe_format(row.get("delta_from_full_rho") if hasattr(row, "get") else None),
                    ])
                
                jackknife_table_para = []
                for row_idx, row in enumerate(jackknife_data_table):
                    para_row = []
                    for cell in row:
                        if row_idx == 0:
                            para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                        else:
                            para_row.append(Paragraph(str(cell), custom_styles["body"]))
                    jackknife_table_para.append(para_row)
                
                jackknife_table = Table(jackknife_table_para, colWidths=[2.5 * inch, 2 * inch, 2 * inch])
                _style_scientific_table(jackknife_table, align_numbers_right=True, alternating_rows=True)
                story.append(jackknife_table)
                story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            
            # PHASE 3: Compress Effect Size Interpretation (max 4 lines)
            effect_interpretation = robustness_data.get("effect_interpretation")
            if effect_interpretation:
                # Split long text into bullets if needed
                interp_lines = effect_interpretation.split('\n')[:4]
                story.append(Paragraph("<b>Biological Interpretation:</b> " + "<br/>".join(interp_lines), custom_styles["body"]))
                story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        else:
            story.append(Paragraph(
                "Robustness analysis unavailable due to insufficient COSMIC overlap.",
                custom_styles["body"]
            ))

        # =====================================================================
        # Section 7: External Validity Assessment
        # =====================================================================
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))
        story.append(Paragraph("7. External Validity Assessment", custom_styles["section"]))
        _add_section_divider(story)
        story.append(Spacer(1, 0.08 * inch))
        
        # 7.1 Null Model Falsification Test
        _add_cond_break(story, MIN_SECTION_HEADER)  # Subsection orphan protection
        story.append(Paragraph("7.1 Null Model Falsification Test", custom_styles["subsection"]))
        story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        
        if null_model_data and isinstance(null_model_data, dict) and len(null_model_data) > 0:
            # Clear explanation for both audiences
            story.append(Paragraph(
                "<b>What this test does:</b> We randomly shuffle the gene pairs 1,000 times and count how often we see similar overlaps. If our observed overlap is much higher than random chance, it suggests a real biological relationship rather than coincidence. <b>Interpretation:</b> Lower p-value (< 0.05) = overlap unlikely due to chance = biological relationship likely.",
                custom_styles["body"]
            ))
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            
            # PHASE 4: Table-First with Interpretation
            _add_cond_break(story, MIN_METRIC_TABLE)
            observed_overlap = null_model_data.get("observed_overlap")
            null_p = null_model_data.get("empirical_p_value")
            n_sim = null_model_data.get("n_simulations", 1000)
            null_mean = null_model_data.get("null_mean_overlap")
            null_std = null_model_data.get("null_std_overlap")
            
            null_interp = "N/A"
            if null_p is not None:
                null_interp = "Significant" if null_p < 0.05 else "Not significant"
            
            null_data_table = [
                ["Metric", "Value", "Interpretation"],
                ["Observed Overlap", _safe_format(observed_overlap), "Actual overlap count"],
                ["Null Mean Overlap", _safe_format(null_mean), "Random expectation"],
                ["Null Std Overlap", _safe_format(null_std), "Random variability"],
                ["Empirical P-Value", _safe_format_p_value(null_p), null_interp],
                ["N Simulations", str(n_sim), "Resampling count"],
            ]
            
            null_table_para = []
            for row_idx, row in enumerate(null_data_table):
                para_row = []
                for cell in row:
                    if row_idx == 0:
                        para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                    else:
                        para_row.append(Paragraph(str(cell), custom_styles["body"]))
                null_table_para.append(para_row)
            
            null_table = Table(null_table_para, colWidths=[2 * inch, 1.5 * inch, 3 * inch])
            _style_scientific_table(null_table, align_numbers_right=True, alternating_rows=True)
            story.append(null_table)
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            
            # PHASE 3: Compress interpretation (max 4 lines)
            if null_p is not None:
                interpretation = (
                    f"After {n_sim} simulations, p={null_p:.4f}. Null hypothesis {'rejected' if null_p < 0.05 else 'accepted'}. "
                    f"{'Overlap unlikely due to chance.' if null_p < 0.05 else 'Overlap could occur by chance.'}"
                )
                story.append(Paragraph(interpretation, custom_styles["body"]))
                story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        else:
            story.append(Paragraph(
                "Null model falsification test unavailable due to insufficient data or missing COSMIC overlap.",
                custom_styles["body"]
            ))
        
        story.append(Spacer(1, 6 / 72.0 * inch))
        
        # 7.2 External Validity Stability
        _add_cond_break(story, MIN_SECTION_HEADER)  # Subsection orphan protection
        story.append(Paragraph("7.2 External Validity Stability", custom_styles["subsection"]))
        story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
        
        if external_validity_data and isinstance(external_validity_data, dict) and len(external_validity_data) > 0:
            # Clear explanation for both audiences
            story.append(Paragraph(
                "<b>What this test does:</b> We resample the data 100 times (bootstrap analysis) to estimate a confidence interval (CI) for the correlation. <b>Interpretation:</b> Narrow CI (< 0.3 width) = stable and reliable results. Wide CI (> 0.6 width) = results may vary with different samples = less confidence in generalizing findings.",
                custom_styles["body"]
            ))
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            
            ci_lower = external_validity_data.get("stability_ci_lower")
            ci_upper = external_validity_data.get("stability_ci_upper")
            n_bootstrap = external_validity_data.get("bootstrap_iterations", 100)
            
            if ci_lower is not None and ci_upper is not None:
                ci_width = ci_upper - ci_lower
                badge_symbol, badge_text = _get_stability_badge(ci_width)
                
                # PHASE 9: Stability Signal Badge
                stability_badge_text = f"<b>STABILITY: {badge_symbol} {badge_text}</b>"
                badge_data = [[Paragraph(stability_badge_text, custom_styles["body"])]]
                badge_table = Table(badge_data, colWidths=[6.5 * inch])
                badge_color = colors.HexColor("#e8f5e9") if ci_width < 0.3 else colors.HexColor("#fff9c4") if ci_width < 0.6 else colors.HexColor("#ffebee")
                badge_table.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, -1), badge_color),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 12),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 12),
                    ("TOPPADDING", (0, 0), (-1, -1), 10),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
                    ("BOX", (0, 0), (-1, -1), 2, colors.HexColor("#888888")),
                ]))
                story.append(badge_table)
                story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
                
                # PHASE 4: Table-First with Interpretation
                _add_cond_break(story, MIN_METRIC_TABLE)
                validity_data_table = [
                    ["Metric", "Value", "Interpretation"],
                    ["95% CI Lower", _safe_format_ci(ci_lower), "Lower bound"],
                    ["95% CI Upper", _safe_format_ci(ci_upper), "Upper bound"],
                    ["CI Width", _safe_format_ci(ci_width), badge_text],
                    ["Bootstrap Iterations", str(n_bootstrap), "Resampling count"],
                ]
                
                validity_table_para = []
                for row_idx, row in enumerate(validity_data_table):
                    para_row = []
                    for cell in row:
                        if row_idx == 0:
                            para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                        else:
                            para_row.append(Paragraph(str(cell), custom_styles["body"]))
                    validity_table_para.append(para_row)
                
                validity_table = Table(validity_table_para, colWidths=[2 * inch, 1.5 * inch, 3 * inch])
                _style_scientific_table(validity_table, align_numbers_right=True, alternating_rows=True)
                story.append(validity_table)
                story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            else:
                story.append(Paragraph(
                    "Stability confidence interval could not be computed due to insufficient data.",
                    custom_styles["body"]
                ))
        else:
            story.append(Paragraph(
                "External validity stability test unavailable due to insufficient data or missing COSMIC overlap.",
                custom_styles["body"]
            ))
        
        story.append(Spacer(1, 6 / 72.0 * inch))
        
        # 7.3 COSMIC Sampling Bias Characterization
        _add_cond_break(story, MIN_SECTION_HEADER)  # Subsection orphan protection
        story.append(Paragraph("7.3 COSMIC Sampling Bias Characterization", custom_styles["subsection"]))
        story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
        
        if bias_analysis_data and isinstance(bias_analysis_data, dict) and len(bias_analysis_data) > 0:
            # Clear explanation for both audiences
            story.append(Paragraph(
                "<b>What this test does:</b> We compare the distribution of recurrence counts between our full dataset and the subset that overlaps with COSMIC. <b>Interpretation:</b> If these distributions are significantly different, COSMIC may be biased toward certain fusion types (e.g., more common or well-studied ones), which should be considered when interpreting results.",
                custom_styles["body"]
            ))
            story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
            
            mw_p = bias_analysis_data.get("mannwhitney_p_value")
            our_median = bias_analysis_data.get("our_recurrence_median")
            cosmic_median = bias_analysis_data.get("cosmic_overlap_median")
            bias_detected = bias_analysis_data.get("bias_detected")
            
            bias_data_table = [
                ["Metric", "Value"],
                ["Mann-Whitney P-Value", _safe_format_p_value(mw_p)],
                ["Our Recurrence Median", _safe_format(our_median)],
                ["COSMIC Overlap Median", _safe_format(cosmic_median)],
                ["Bias Detected", "Yes" if bias_detected else "No" if bias_detected is not None else "N/A"],
            ]
            
            _add_cond_break(story, MIN_METRIC_TABLE)
            bias_table = Table(bias_data_table, colWidths=[3.5 * inch, 3 * inch])
            _style_scientific_table(bias_table, align_numbers_right=True, alternating_rows=True)
            story.append(Spacer(1, SPACE_BEFORE_TABLE / 72.0 * inch))
            story.append(bias_table)
            story.append(Spacer(1, SPACE_AFTER_TABLE / 72.0 * inch))
            
            # PHASE 3: Compress interpretation (max 4 lines)
            if mw_p is not None:
                interpretation = (
                    f"Sampling bias {'detected' if bias_detected else 'not detected'} (Mann-Whitney p={mw_p:.4f}). "
                    f"{'COSMIC may be biased toward certain fusion types.' if bias_detected else 'No significant sampling bias detected.'}"
                )
                story.append(Paragraph(interpretation, custom_styles["body"]))
                story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        else:
            story.append(Paragraph(
                "COSMIC bias quantification unavailable due to missing COSMIC data.",
                custom_styles["body"]
            ))
        
        # =====================================================================
        # Section 8: Final Evidence Synthesis & Claim Strength
        # =====================================================================
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))
        story.append(Paragraph("8. Final Evidence Synthesis & Claim Strength", custom_styles["section"]))
        _add_section_divider(story)
        story.append(Spacer(1, 0.08 * inch))
        
        # PHASE 3: Compress paragraph (max 4 lines)
        story.append(Paragraph(
            "<b>Evidence Synthesis:</b> Combines correlation strength, significance testing, and stability into single classification.",
            custom_styles["body"]
        ))
        story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        
        # PHASE 7: Evidence Synthesis Matrix (before claim strength)
        _add_cond_break(story, MIN_METRIC_TABLE)
        evidence_synthesis_data = [["Evidence Layer", "Signal", "Confidence", "Limitation"]]
        
        # Data Integrity
        data_integrity_result = "PASS" if exit_code == 0 else "FAIL"
        data_integrity_confidence = "HIGH" if exit_code == 0 else "LOW"
        data_integrity_limitation = "None" if exit_code == 0 else "Schema validation failed"
        evidence_synthesis_data.append(["Data Integrity", data_integrity_result, data_integrity_confidence, data_integrity_limitation])
        
        # Distribution Validity
        dist_result = "STRONG"
        dist_confidence = "HIGH"
        dist_limitation = "None"
        if diagnostic_data:
            dist = diagnostic_data.get("distribution", {})
            if not dist or len(dist) == 0:
                dist_result = "LIMITED"
                dist_confidence = "LOW"
                dist_limitation = "Insufficient distribution data"
        evidence_synthesis_data.append(["Distribution", dist_result, dist_confidence, dist_limitation])
        
        # COSMIC Correlation
        cosmic_result = "N/A"
        cosmic_confidence = "N/A"
        cosmic_limitation = "N/A"
        if diagnostic_data:
            cosmic = diagnostic_data.get("cosmic", {})
            spearman_rho = cosmic.get("spearman_rho")
            overlap_count = cosmic.get("overlap_count", 0)
            if spearman_rho is not None:
                abs_rho = abs(float(spearman_rho))
                if abs_rho >= 0.7:
                    cosmic_result = "STRONG"
                    cosmic_confidence = "HIGH" if overlap_count >= 30 else "MEDIUM"
                    cosmic_limitation = "Low N" if overlap_count < 30 else "None"
                elif abs_rho >= 0.4:
                    cosmic_result = "MODERATE"
                    cosmic_confidence = "MEDIUM"
                    cosmic_limitation = "Low N" if overlap_count < 30 else "Moderate correlation"
                else:
                    cosmic_result = "WEAK"
                    cosmic_confidence = "LOW"
                    cosmic_limitation = "Weak correlation or small sample size"
            else:
                cosmic_result = "N/A"
                cosmic_confidence = "N/A"
                cosmic_limitation = "No COSMIC overlap data"
        evidence_synthesis_data.append(["COSMIC Overlap", cosmic_result, cosmic_confidence, cosmic_limitation])
        
        # Null Model
        null_result = "N/A"
        null_confidence = "N/A"
        null_limitation = "N/A"
        if null_model_data:
            null_p = null_model_data.get("empirical_p_value")
            if null_p is not None:
                null_result = "STRONG" if null_p < 0.05 else "WEAK"
                null_confidence = "HIGH" if null_p < 0.05 else "LOW"
                null_limitation = "Simulation assumptions" if null_p < 0.05 else "Not significant"
            else:
                null_result = "N/A"
                null_confidence = "N/A"
                null_limitation = "Insufficient data"
        evidence_synthesis_data.append(["Null Model", null_result, null_confidence, null_limitation])
        
        # External Validity
        ext_validity_result = "N/A"
        ext_validity_confidence = "N/A"
        ext_validity_limitation = "N/A"
        if external_validity_data:
            ci_lower = external_validity_data.get("stability_ci_lower")
            ci_upper = external_validity_data.get("stability_ci_upper")
            if ci_lower is not None and ci_upper is not None:
                ci_width = ci_upper - ci_lower
                if ci_width < 0.3:
                    ext_validity_result = "STRONG"
                    ext_validity_confidence = "HIGH"
                    ext_validity_limitation = "None"
                elif ci_width < 0.6:
                    ext_validity_result = "MODERATE"
                    ext_validity_confidence = "MEDIUM"
                    ext_validity_limitation = "Wide CI"
                else:
                    ext_validity_result = "WEAK"
                    ext_validity_confidence = "LOW"
                    ext_validity_limitation = "Very wide CI"
            else:
                ext_validity_result = "N/A"
                ext_validity_confidence = "N/A"
                ext_validity_limitation = "Insufficient data"
        evidence_synthesis_data.append(["External Stability", ext_validity_result, ext_validity_confidence, ext_validity_limitation])
        
        # Convert to Paragraphs
        evidence_table_para = []
        for row_idx, row in enumerate(evidence_synthesis_data):
            para_row = []
            for cell in row:
                if row_idx == 0:
                    para_row.append(Paragraph(f"<b>{str(cell)}</b>", custom_styles["body"]))
                else:
                    para_row.append(Paragraph(str(cell), custom_styles["body"]))
            evidence_table_para.append(para_row)
        
        evidence_table = Table(evidence_table_para, colWidths=[2 * inch, 1.5 * inch, 1.5 * inch, 1.5 * inch])
        _style_scientific_table(evidence_table, align_numbers_right=False, alternating_rows=True)
        story.append(evidence_table)
        story.append(Spacer(1, TABLE_SPACE / 72.0 * inch))
        
        # Power Analysis
        if diagnostic_data:
            cosmic = diagnostic_data.get("cosmic", {})
            overlap_count = cosmic.get("overlap_count", 0)
            if isinstance(overlap_count, (int, float)) and overlap_count > 0:
                story.append(Paragraph(
                    "<b>8.1 Statistical Power Analysis:</b>",
                    custom_styles["subsection"]
                ))
                story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
                # PHASE 3: Compress to max 4 lines
                power_text = (
                    f"Overlap n={int(overlap_count)}: <b>extremely low</b> power (<20% for ρ>0.7). "
                    f"High Type II error risk. Correlations highly uncertain. "
                    f"Minimum 30-50 pairs recommended for reliable analysis."
                )
                story.append(Paragraph(power_text, custom_styles["body"]))
                story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        
        if claim_strength_data and isinstance(claim_strength_data, dict) and len(claim_strength_data) > 0:
            classification = claim_strength_data.get("classification", "N/A")
            confidence_level = claim_strength_data.get("confidence_level", "N/A")
            interpretation = claim_strength_data.get("interpretation", "")
            
            # Highlight box for classification
            claim_lines = [
                f"<b>Classification:</b> {classification}",
                f"<b>Confidence Level:</b> {confidence_level}",
            ]
            story.append(Paragraph("8.2 Claim Strength Classification", custom_styles["subsection"]))
            story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
            story.append(_build_highlight_box(claim_lines, custom_styles))
            story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
            
            # Interpretation paragraph
            if interpretation:
                story.append(Paragraph(interpretation, custom_styles["body"]))
            else:
                story.append(Paragraph(
                    "Scientific claim strength assessment completed. See classification above.",
                    custom_styles["body"]
                ))
            
            # What Would Change Conclusion Section
            story.append(Spacer(1, 6 / 72.0 * inch))
            story.append(Paragraph(
                "<b>8.3 What Evidence Would Overturn This Conclusion?</b>",
                custom_styles["subsection"]
            ))
            story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
            
            if diagnostic_data:
                cosmic = diagnostic_data.get("cosmic", {})
                overlap_count = cosmic.get("overlap_count", 0)
                spearman_rho = cosmic.get("spearman_rho")
                null_p = None
                if advanced_inference and advanced_inference.get("null_model"):
                    null_p = advanced_inference["null_model"].get("empirical_p_value")
                
                # PHASE 3: Compress to bullet summary (max 4 lines per item)
                reversal_text = (
                    "<b>Evidence that would overturn conclusion:</b><br/>"
                    "• Increased sample size (n≥30-50) with decreased correlation<br/>"
                    "• Replication failure in independent datasets<br/>"
                    "• Systematic COSMIC bias discovery<br/>"
                    "• Technical artifacts driving correlation<br/>"
                    "• Null model failure (p>0.05) with larger samples"
                )
                story.append(Paragraph(reversal_text, custom_styles["body"]))
                story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
                
                # PHASE 3: Compress replication text
                replication_text = (
                    "<b>Ideal replication:</b> Independent dataset, larger COSMIC subset (50+ pairs), "
                    "alternative databases (OncoKB, cBioPortal), stratified analysis, prospective validation."
                )
                story.append(Paragraph(replication_text, custom_styles["body"]))
                story.append(Spacer(1, PARAGRAPH_SPACE / 72.0 * inch))
        else:
            story.append(Paragraph(
                "Scientific claim strength classification unavailable due to missing statistical metrics.",
                custom_styles["body"]
            ))
        
        # =====================================================================
        # Section 9: Limitations
        # =====================================================================
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SUBSECTION_SPACE / 72.0 * inch))
        story.append(Paragraph("9. Limitations", custom_styles["section"]))
        _add_section_divider(story)
        story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))

        if limitations_text:
            story.append(Paragraph(limitations_text, custom_styles["body"]))
        else:
            story.append(Paragraph(
                "Limitations section unavailable.",
                custom_styles["body"]
            ))

        # Flow to Section 10 on same page if space allows
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SECTION_SPACE / 72.0 * inch))

        # =====================================================================
        # Section 10: Methods Transparency
        # =====================================================================
        story.append(Paragraph("10. Methods Transparency", custom_styles["section"]))
        _add_section_divider(story)
        story.append(Spacer(1, 4 / 72.0 * inch))

        if scientific_narrative:
            story.append(Paragraph(scientific_narrative, custom_styles["body"]))
            story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))
        
        if confidence_statement:
            story.append(Paragraph(confidence_statement, custom_styles["body"]))
            story.append(Spacer(1, SPACE_BETWEEN_PARAGRAPHS / 72.0 * inch))

        # Flow to Section 11 on same page if space allows
        _add_cond_break(story, MIN_SECTION_HEADER)  # Orphan header protection
        story.append(Spacer(1, SECTION_SPACE / 72.0 * inch))

        # =====================================================================
        # Section 11: Reproducibility Metadata
        # =====================================================================
        story.append(Paragraph("11. Reproducibility Metadata", custom_styles["section"]))
        _add_section_divider(story)
        story.append(Spacer(1, 4 / 72.0 * inch))

        # COSMIC source
        cosmic_source = "N/A"
        cosmic_version_prov = "N/A"
        if diagnostic_data:
            cosmic = diagnostic_data.get("cosmic", {})
            cosmic_source = cosmic.get("cosmic_reference_source", "N/A")
            cosmic_version_prov = cosmic.get("cosmic_reference_version", "N/A")

        # Python version
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

        # Library versions
        numpy_version = "N/A"
        pandas_version = "N/A"
        scipy_version = "N/A"
        try:
            import numpy
            numpy_version = numpy.__version__
        except ImportError:
            pass
        try:
            import pandas
            pandas_version = pandas.__version__
        except ImportError:
            pass
        try:
            import scipy
            scipy_version = scipy.__version__
        except ImportError:
            pass

        # Pipeline version
        try:
            from week2_validation import __version__ as pipeline_version
        except ImportError:
            pipeline_version = "unknown"

        # Build provenance table
        prov_data = [
            ["Field", "Value"],
            ["COSMIC Source", str(cosmic_source)],
            ["COSMIC Version", str(cosmic_version_prov)],
            ["Python Version", python_version],
            ["NumPy Version", numpy_version],
            ["Pandas Version", pandas_version],
            ["SciPy Version", scipy_version],
            ["Pipeline Version", pipeline_version],
        ]

        prov_table = Table(prov_data, colWidths=[3 * inch, 3.5 * inch])
        _style_scientific_table(prov_table, align_numbers_right=False, alternating_rows=True)
        story.append(Spacer(1, SPACE_BEFORE_TABLE / 72.0 * inch))
        story.append(prov_table)
        story.append(Spacer(1, SPACE_AFTER_TABLE / 72.0 * inch))

        story.append(Spacer(1, 0.15 * inch))
        story.append(Paragraph(
            "<i>Generated by Data Integrity &amp; Statistical Validation layer</i>",
            custom_styles["body"]
        ))

        # Certification page
        story.append(PageBreak())
        story.extend(_build_certification_flowables(
            dataset_hash=dataset_hash,
            quality_gates=qg_for_cert,
            overall_approval=overall_approval,
            validation_results_hash=validation_results_hash,
            styles=custom_styles,
        ))

        
        try:
            # Two-pass build: first pass to get total page count
            with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
                temp_path = tmp.name

            doc_count = SimpleDocTemplate(
                temp_path,
                pagesize=letter,
                rightMargin=0.75 * inch,
                leftMargin=0.75 * inch,
                topMargin=0.75 * inch,
                bottomMargin=0.75 * inch,
            )
            doc_count.onFirstPage = counting_footer
            doc_count.onLaterPages = counting_footer
            doc_count.build(copy.deepcopy(story))
        finally:
            try:
                Path(temp_path).unlink(missing_ok=True)
            except Exception:
                pass

        # Second pass: final PDF
        doc_final = SimpleDocTemplate(
            str(output_pdf_path),
            pagesize=letter,
            rightMargin=0.75 * inch,
            leftMargin=0.75 * inch,
            topMargin=0.75 * inch,
            bottomMargin=0.75 * inch,
        )
        doc_final.onFirstPage = full_header_footer
        doc_final.onLaterPages = full_header_footer
        doc_final.build(copy.deepcopy(story))

        return output_pdf_path

    except Exception as e:
        print(f"Warning: PDF generation failed: {e}", file=sys.stderr)
        return None
