"""
Signature metadata collection and reverse linking.

Handles fetching signature metadata from Notion and creating reverse links
from signature pages to source pages (literature, email, experiment, dataset).
"""

from __future__ import annotations

import json
from typing import Any, Dict, List

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
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
    Create reverse link from signature page to source page (literature/email/experiment/dataset).

    Updates the signature page's "Source Papers" relation to include the source page.
    """
    cfg = get_config().notion
    try:
        # Fetch current signature page
        sig_page = fetch_notion_page(signature_page_id)
        props = sig_page.get("properties", {}) or {}
        current_rel = props.get("Source Papers", {}).get("relation", []) or []
        existing_ids = {
            r.get("id", "").replace("-", "") for r in current_rel if r.get("id")
        }

        # Check if source page is already linked
        source_page_id_clean = source_page_id.replace("-", "")
        if source_page_id_clean in existing_ids:
            return  # Already linked, skip

        # Add source page to relation
        updated_rel = [{"id": source_page_id_clean}] + current_rel

        # Update signature page
        payload = {
            "properties": {
                "Source Papers": {"relation": updated_rel},
            },
        }

        url = f"{cfg.base_url}/pages/{signature_page_id}"
        resp = requests.patch(
            url,
            headers=notion_headers(),
            data=json.dumps(payload),
            timeout=30,
        )

        if resp.status_code >= 300:
            logger.error(
                "[METADATA] Failed to update signature page %s reverse link: %s",
                signature_page_id,
                resp.text,
            )
        else:
            logger.info(
                "[METADATA] Added reverse link: signature %s -> source %s",
                signature_page_id,
                source_page_id,
            )

    except Exception as e:
        logger.error(
            "[METADATA] Error creating reverse link for signature %s -> source %s: %r",
            signature_page_id,
            source_page_id,
            e,
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

