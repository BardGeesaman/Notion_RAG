"""
Detect omics type from file name and CSV headers.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def detect_omics_type(file_path: str) -> Optional[str]:
    """
    Detect omics type from filename and CSV headers.

    Returns one of: "lipidomics", "metabolomics", "proteomics", "transcriptomics", or None.
    """
    path = Path(file_path)
    name_lower = path.name.lower()

    # Filename heuristics
    if "lipid" in name_lower:
        return "lipidomics"
    if "metabol" in name_lower:
        return "metabolomics"
    if "prote" in name_lower or "protein" in name_lower:
        return "proteomics"
    if "transcript" in name_lower or "gene" in name_lower or "rna" in name_lower:
        return "transcriptomics"

    # Header heuristics
    try:
        with path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            headers = next(reader, [])
            headers_lower = [h.lower() for h in headers]
    except Exception as exc:
        logger.debug("[OMICS-DETECT] Could not read headers for %s: %r", file_path, exc)
        return None

    header_str = ",".join(headers_lower)

    if "lipid" in header_str:
        return "lipidomics"
    if "metabolite" in header_str or "compound" in header_str:
        return "metabolomics"
    if "protein" in header_str or "uniprot" in header_str:
        return "proteomics"
    if (
        "gene" in header_str
        or "log2foldchange" in header_str
        or "ensembl" in header_str
    ):
        return "transcriptomics"

    return None

