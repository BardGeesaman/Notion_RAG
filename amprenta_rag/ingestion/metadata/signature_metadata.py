"""
Signature metadata collection and reverse linking.

Handles fetching signature metadata from Notion and creating reverse links
from signature pages to source pages (literature, email, experiment, dataset).
"""

from __future__ import annotations

import json
from typing import Any, Dict, List

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ingestion.metadata.helpers import (
    fetch_notion_page,
    get_multi_names,
    get_select_name,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def enforce_signature_reverse_link(
    signature_page_id: str, source_page_id: str
) -> None:
    """
    DEPRECATED: Notion support has been removed.
    
    This function previously created reverse links from signature pages to source pages.
    Postgres is now the source of truth.
    
    TODO: Phase 3 - Implement Postgres-based reverse linking if needed.
    """
    logger.debug(
        "[METADATA] enforce_signature_reverse_link() deprecated - Notion support removed"
    )


def collect_signature_metadata(signature_ids: List[str]) -> Dict[str, Any]:
    """
    Fetch Lipid Signature pages and collect their metadata.

    For each signature page ID, fetches the page and extracts:
    - Short ID (rich_text or title)
    - Biomarker Role (multi_select)
    - Phenotype Axes (multi_select)
    - Data Ownership (select)

    Args:
        signature_ids: List of Notion page IDs (without dashes) for signature pages

    Returns:
        Dictionary with:
        - sig_short_ids: List[str] - Short IDs from signature pages
        - sig_roles: List[str] - Biomarker roles (aggregated)
        - sig_axes: List[str] - Phenotype axes (aggregated)
        - sig_ownership: List[str] - Data ownership values (aggregated)
    """
    sig_short_ids: List[str] = []
    sig_roles: List[str] = []
    sig_axes: List[str] = []
    sig_ownership: List[str] = []

    if not signature_ids:
        return {
            "sig_short_ids": [],
            "sig_roles": [],
            "sig_axes": [],
            "sig_ownership": [],
        }

    for sig_id in signature_ids:
        try:
            sig_page = fetch_notion_page(sig_id)
            sig_props = sig_page.get("properties", {}) or {}

            # Extract Short ID from rich_text or title, fallback to Name property
            short_id = None
            short_id_prop = sig_props.get("Short ID", {})
            if "rich_text" in short_id_prop:
                rich_text = short_id_prop.get("rich_text", []) or []
                if rich_text:
                    short_id = rich_text[0].get("plain_text", "").strip()
            elif "title" in short_id_prop:
                title = short_id_prop.get("title", []) or []
                if title:
                    short_id = title[0].get("plain_text", "").strip()

            # Fallback to Name property if Short ID is missing/empty
            if not short_id:
                name_prop = sig_props.get("Name", {})
                if "title" in name_prop:
                    title = name_prop.get("title", []) or []
                    if title:
                        short_id = title[0].get("plain_text", "").strip()
                elif "rich_text" in name_prop:
                    rich_text = name_prop.get("rich_text", []) or []
                    if rich_text:
                        short_id = rich_text[0].get("plain_text", "").strip()

            if short_id:
                sig_short_ids.append(short_id)

            # Extract Biomarker Role (multi-select)
            roles = get_multi_names(sig_props, "Biomarker Role")
            sig_roles.extend(roles)

            # Extract Phenotype Axes (multi-select)
            axes = get_multi_names(sig_props, "Phenotype Axes")
            sig_axes.extend(axes)

            # Extract Data Ownership (select)
            ownership = get_select_name(sig_props, "Data Ownership")
            if ownership:
                sig_ownership.append(ownership)

        except Exception as e:
            logger.warning(
                "[METADATA] Error fetching signature page %s: %r",
                sig_id,
                e,
            )
            continue

    return {
        "sig_short_ids": sig_short_ids,
        "sig_roles": list(set(sig_roles)),  # Deduplicate
        "sig_axes": list(set(sig_axes)),  # Deduplicate
        "sig_ownership": list(set(sig_ownership)),  # Deduplicate
    }

