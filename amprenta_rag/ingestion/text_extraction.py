# amprenta_rag/ingestion/text_extraction.py

from __future__ import annotations

from io import BytesIO
from typing import List, Optional

import re

from pypdf import PdfReader

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_text_from_pdf_bytes(data: bytes) -> str:
    """
    Extract text from a PDF using pypdf.
    """
    reader = PdfReader(BytesIO(data))
    text_parts: List[str] = []
    for page in reader.pages:
        t = page.extract_text() or ""
        if t:
            text_parts.append(t)
    return "\n".join(text_parts)


def extract_text_from_attachment_bytes(
    content_type: Optional[str],
    data: bytes,
    filename: str = "",
) -> str:
    """
    Convert different attachment types to text.

    - PDFs: use PdfReader
    - text/*: decode as UTF-8
    - everything else: skip
    """
    ct = (content_type or "").lower()

    if ct == "application/pdf":
        return extract_text_from_pdf_bytes(data)

    if ct.startswith("text/"):
        try:
            return data.decode("utf-8", errors="ignore")
        except Exception as e:
            logger.warning("⚠️ Failed to decode text attachment %s: %s", filename, e)
            return ""

    logger.info("ℹ️ Skipping non-text attachment %r (%s)", filename, content_type)
    return ""


def html_to_text(html: str) -> str:
    """
    Very simple HTML → text conversion for Zotero notes or HTML chunks.
    """
    # Strip comments
    html = re.sub(r"<!--.*?-->", " ", html, flags=re.DOTALL)
    # <br> → newline
    html = re.sub(r"<br\s*/?>", "\n", html, flags=re.IGNORECASE)
    # Remove tags
    text = re.sub(r"<[^>]+>", " ", html)
    # Normalize whitespace
    text = re.sub(r"\s+", " ", text)
    return text.strip()


_BOILERPLATE_PATTERNS = [
    "all rights reserved",
    "no responsibility or liability",
    "creativecommons.org/licenses",
    "this is an open access article",
    "distributed under the terms of the creative commons",
    "medrxiv preprint",
    "biorxiv preprint",
    "author/funder",
    "supplementary table",
    "supplementary figure",
    "view this article online",
]


def is_boilerplate(text: str) -> bool:
    """
    Heuristic to skip obvious license/footer boilerplate.
    """
    t = text.lower()
    return any(pat in t for pat in _BOILERPLATE_PATTERNS)