"""File extraction utilities (DOCX/PPTX/CSV/XLSX parsers)."""

from amprenta_rag.extraction.parsers import (
    get_parser,
    parse_csv,
    parse_docx,
    parse_excel,
    parse_pptx,
)

__all__ = ["get_parser", "parse_csv", "parse_docx", "parse_excel", "parse_pptx"]


