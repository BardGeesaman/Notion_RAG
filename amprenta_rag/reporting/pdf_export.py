"""
PDF export utilities for evidence reports and markdown content.

This module converts markdown evidence reports into styled PDF documents using
WeasyPrint, with a professional layout suitable for sharing with stakeholders.
"""

from __future__ import annotations

from io import BytesIO
from pathlib import Path
from typing import Optional
from datetime import datetime

import markdown
from weasyprint import HTML

from amprenta_rag.reporting.evidence_report import format_evidence_report

CSS_STYLE = """
@page {
    size: A4;
    margin: 1in;
    @bottom-center {
        content: "Page " counter(page) " of " counter(pages);
        font-family: Arial, sans-serif;
        font-size: 10px;
        color: #666;
    }
}

body {
    font-family: Georgia, serif;
    font-size: 12px;
    line-height: 1.6;
    color: #222;
}

h1, h2, h3, h4 {
    font-family: Arial, sans-serif;
    color: #1f3b57;
    margin-top: 0.4in;
    margin-bottom: 0.15in;
}

h1 { font-size: 24px; border-bottom: 2px solid #1f3b57; padding-bottom: 6px; }
h2 { font-size: 20px; }
h3 { font-size: 16px; }
h4 { font-size: 14px; }

p {
    margin: 0 0 0.15in 0;
}

table {
    width: 100%;
    border-collapse: collapse;
    margin: 0.2in 0;
    font-size: 11px;
}

th, td {
    border: 1px solid #ccc;
    padding: 6px 8px;
    text-align: left;
}

th {
    background-color: #f0f4f8;
    font-weight: bold;
}

tr:nth-child(even) {
    background: #fafafa;
}

.section {
    border: 1px solid #e0e0e0;
    padding: 12px;
    margin: 0.2in 0;
    border-radius: 4px;
    background: #fdfdfd;
}

.meta {
    font-size: 11px;
    color: #555;
    margin-bottom: 0.1in;
}

.timestamp {
    font-size: 10px;
    color: #888;
}
"""


def markdown_to_html(md_content: str) -> str:
    """Convert markdown text to a styled HTML document."""
    html_body = markdown.markdown(md_content, extensions=["tables", "fenced_code"])
    timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    html = f"""
    <html>
    <head>
        <meta charset="utf-8">
        <style>{CSS_STYLE}</style>
    </head>
    <body>
        <div class="timestamp">Generated: {timestamp}</div>
        {html_body}
    </body>
    </html>
    """
    return html


def html_to_pdf(html_content: str) -> bytes:
    """Render HTML to PDF bytes using WeasyPrint."""
    pdf_bytes = HTML(string=html_content).write_pdf()
    return pdf_bytes


def markdown_to_pdf(md_content: str, output_path: Optional[Path] = None) -> bytes:
    """Convert markdown content to PDF bytes; optionally write to a file."""
    html = markdown_to_html(md_content)
    pdf_bytes = html_to_pdf(html)
    if output_path:
        output_path.write_bytes(pdf_bytes)
    return pdf_bytes


def export_evidence_report_to_pdf(report, output_path: Optional[Path] = None) -> bytes:
    """
    Format an EvidenceReport to markdown and export as PDF.
    """
    md_content = format_evidence_report(report, include_metadata=True)
    return markdown_to_pdf(md_content, output_path=output_path)

