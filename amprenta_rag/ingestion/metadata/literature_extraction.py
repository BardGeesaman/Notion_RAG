"""
Literature metadata extraction.

DEPRECATED: Notion support has been removed. Postgres is now the source of truth.
"""

from __future__ import annotations

import re
from typing import Any, Dict

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ingestion.zotero_api import ZoteroItem

logger = get_logger(__name__)


def get_literature_semantic_metadata(
    parent_page_id: str, item: ZoteroItem
) -> Dict[str, Any]:
    """
    DEPRECATED: Notion support has been removed.
    
    Returns minimal metadata structure from Zotero item only.
    """
    logger.debug(
        "[METADATA] get_literature_semantic_metadata() - Notion support removed, using Zotero data only"
    )
    
    # --- Simple doc_type inference from Zotero item_type and tags ---
    doc_type = "JournalArticle"
    if item.item_type.lower() == "book":
        doc_type = "Book"
    elif item.item_type.lower() == "thesis":
        doc_type = "Thesis"
    elif any("review" in t.lower() for t in item.tags):
        doc_type = "Review"
    elif any("trial" in t.lower() for t in item.tags):
        doc_type = "ClinicalTrialDesign"

    # Extract year from date string if available
    year_val = None
    if item.date:
        m = re.search(r"(19|20)\d{2}", item.date)
        if m:
            year_val = int(m.group(0))

    return {
        "doc_id": f"ZOTERO:{item.key}",
        "doc_source": "Zotero",
        "doc_source_subtype": "Attachment",
        "doc_type": doc_type,
        "diseases": [],
        "targets": [],
        "modality": [],
        "stage": None,
        "model_systems": [],
        "biomarker_role": [],
        "year": year_val,
        "journal": item.journal,
        "importance": None,
        "manual_tags": item.tags,
        "lipids_raw": [],
        "lipids": [],
        "lipid_classes": [],
        "lipid_signatures": [],
        "lipid_signature_role": [],
        "phenotype_axes": [],
        "signature_ownership": [],
        "matrix": [],
        "treatment_arms": [],
    }
