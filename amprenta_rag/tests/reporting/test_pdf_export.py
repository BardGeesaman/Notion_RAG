from __future__ import annotations

import sys
from types import ModuleType
import pytest
from pathlib import Path

# Mock weasyprint BEFORE importing the module under test
fake_weasyprint = ModuleType("weasyprint")
fake_weasyprint.HTML = None  # Will be mocked in tests
sys.modules["weasyprint"] = fake_weasyprint

from amprenta_rag.reporting import pdf_export as pe


def test_markdown_to_html():
    md = "# Hello\n\n- item 1"
    html = pe.markdown_to_html(md)
    assert "<html>" in html
    assert "<h1>Hello</h1>" in html
    assert "<li>item 1</li>" in html
    assert "Generated:" in html


def test_html_to_pdf(monkeypatch):
    class FakeHTML:
        def __init__(self, string=None):
            self.string = string

        def write_pdf(self):
            return b"PDF_BYTES"

    # Patch the imported HTML class in pdf_export module
    monkeypatch.setattr(pe, "HTML", FakeHTML)
    
    pdf = pe.html_to_pdf("<html></html>")
    assert pdf == b"PDF_BYTES"


def test_markdown_to_pdf(monkeypatch, tmp_path):
    # Mock lower-level functions
    monkeypatch.setattr(pe, "markdown_to_html", lambda md: "<html>" + md + "</html>")
    monkeypatch.setattr(pe, "html_to_pdf", lambda html: b"MOCK_PDF")

    # 1. Bytes return
    pdf = pe.markdown_to_pdf("content")
    assert pdf == b"MOCK_PDF"

    # 2. File write
    out_file = tmp_path / "report.pdf"
    pdf2 = pe.markdown_to_pdf("content", output_path=out_file)
    assert pdf2 == b"MOCK_PDF"
    assert out_file.read_bytes() == b"MOCK_PDF"


def test_export_evidence_report_to_pdf(monkeypatch):
    class FakeReport:
        pass

    monkeypatch.setattr(pe, "format_evidence_report", lambda r, include_metadata=True: "# Evidence")
    monkeypatch.setattr(pe, "markdown_to_pdf", lambda md, output_path=None: b"EVIDENCE_PDF")

    pdf = pe.export_evidence_report_to_pdf(FakeReport())
    assert pdf == b"EVIDENCE_PDF"
