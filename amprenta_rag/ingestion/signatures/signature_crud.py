"""
Signature page CRUD operations.

Handles creating, finding, and updating signature pages in Notion.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature
from amprenta_rag.ingestion.signatures.short_id import generate_signature_short_id

logger = get_logger(__name__)

def notion_headers() -> Dict[str, str]:
    """DEPRECATED: Notion support removed. Returns empty headers dict."""
    logger.debug("[SIGNATURES][SIGNATURE-CRUD] notion_headers() deprecated - Notion support removed")
    return {}


def find_or_create_signature_page(
    signature: Signature,
    signature_type: str = "Literature-derived",
    data_ownership: str = "Public",
    version: Optional[str] = None,
    description: Optional[str] = None,
) -> Optional[str]:
    """
    Find or create a Lipid Signature page in Notion.

    Args:
        signature: Signature object from loader
        signature_type: "Consortium", "Literature-derived", "Open Dataset", or "Other"
        data_ownership: "Public" or appropriate value
        version: Optional version string
        description: Optional description text

    Returns:
        Notion page ID (with dashes) or None if creation failed
    """

    cfg = get_config()

    # Check if signature DB ID is configured
    if not hasattr(cfg.notion, "signature_db_id") or not cfg.notion.signature_db_id:
        logger.warning(
            "[INGEST][SIGNATURES] Lipid Signatures database ID not configured. "
            "Set NOTION_SIGNATURE_DB_ID in config."
        )
        return None

    db_id = cfg.notion.signature_db_id
    short_id = generate_signature_short_id(signature.name, version)

    # First, try to find existing page by Short ID or Name
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"

        # Try Short ID first
        payload = {
            "filter": {
                "or": [
                    {
                        "property": "Short ID",
                        "rich_text": {"equals": short_id},
                    },
                    {
                        "property": "Name",
                        "title": {"equals": signature.name},
                    },
                ],
            },
            "page_size": 1,
        }

        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()

        results = resp.json().get("results", [])
        if results:
            page_id = results[0].get("id", "")
            logger.debug(
                "[INGEST][SIGNATURES] Found existing signature page for %s: %s",
                signature.name,
                page_id,
            )
            # Update existing page with any missing fields
            update_signature_page_if_needed(
                page_id,
                signature,
                signature_type,
                data_ownership,
                version,
                description,
                short_id,
            )
            return page_id
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error querying for existing signature %s: %r",
            signature.name,
            e,
        )

    # Create new page
    try:
        url = f"{cfg.notion.base_url}/pages"

        properties: Dict[str, Any] = {
            "Name": {
                "title": [{"text": {"content": signature.name}}],
            },
            "Short ID": {
                "rich_text": [{"text": {"content": short_id}}],
            },
            "Status": {
                "select": {"name": "Active"},
            },
            "Signature Type": {
                "select": {"name": signature_type},
            },
            "Data Ownership": {
                "select": {"name": data_ownership},
            },
        }

        if version:
            properties["Version"] = {
                "rich_text": [{"text": {"content": version}}],
            }

        if description or signature.description:
            desc = description or signature.description or ""
            properties["Description"] = {
                "rich_text": [{"text": {"content": desc[:2000]}}],  # Limit length
            }

        payload = {
            "parent": {"database_id": db_id},
            "properties": properties,
        }

        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()

        new_page = resp.json()
        page_id = new_page.get("id", "")
        logger.info(
            "[INGEST][SIGNATURES] Created new signature page for %s: %s (Short ID: %s)",
            signature.name,
            page_id,
            short_id,
        )
        return page_id
    except Exception as e:
        logger.error(
            "[INGEST][SIGNATURES] Error creating signature page for %s: %r",
            signature.name,
            e,
        )
        return None


def update_signature_page_if_needed(
    page_id: str,
    signature: Signature,
    signature_type: str,
    data_ownership: str,
    version: Optional[str],
    description: Optional[str],
    short_id: str,
) -> None:
    """
    Update an existing signature page with missing fields.

    Only updates fields that are empty/missing, preserving existing data.
    """
    cfg = get_config()

    try:
        # Fetch current page to check what's missing
        url = f"{cfg.notion.base_url}/pages/{page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()

        page = resp.json()
        props = page.get("properties", {}) or {}

        updates: Dict[str, Any] = {}

        # Update Short ID if missing
        short_id_prop = props.get("Short ID", {})
        short_id_text = ""
        if "rich_text" in short_id_prop:
            rich_text = short_id_prop.get("rich_text", []) or []
            short_id_text = "".join(rt.get("plain_text", "") for rt in rich_text)
        if not short_id_text.strip():
            updates["Short ID"] = {"rich_text": [{"text": {"content": short_id}}]}

        # Update Version if missing
        if version:
            version_prop = props.get("Version", {})
            version_text = ""
            if "rich_text" in version_prop:
                rich_text = version_prop.get("rich_text", []) or []
                version_text = "".join(rt.get("plain_text", "") for rt in rich_text)
            if not version_text.strip():
                updates["Version"] = {"rich_text": [{"text": {"content": version}}]}

        # Update Description if missing
        desc = description or signature.description
        if desc:
            desc_prop = props.get("Description", {})
            desc_text = ""
            if "rich_text" in desc_prop:
                rich_text = desc_prop.get("rich_text", []) or []
                desc_text = "".join(rt.get("plain_text", "") for rt in rich_text)
            if not desc_text.strip():
                updates["Description"] = {
                    "rich_text": [{"text": {"content": desc[:2000]}}],
                }

        # Update Status if missing
        status_prop = props.get("Status", {})
        if not status_prop.get("select"):
            updates["Status"] = {"select": {"name": "Active"}}

        # Update Signature Type if missing
        sig_type_prop = props.get("Signature Type", {})
        if not sig_type_prop.get("select"):
            updates["Signature Type"] = {"select": {"name": signature_type}}

        # Update Data Ownership if missing
        ownership_prop = props.get("Data Ownership", {})
        if not ownership_prop.get("select"):
            updates["Data Ownership"] = {"select": {"name": data_ownership}}

        # Apply updates if any
        if updates:
            url = f"{cfg.notion.base_url}/pages/{page_id}"
            payload = {"properties": updates}
            resp = requests.patch(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            logger.debug(
                "[INGEST][SIGNATURES] Updated signature page %s with missing fields",
                page_id,
            )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error updating signature page %s: %r",
            page_id,
            e,
        )


def update_signature_modalities(
    page_id: str,
    modalities: List[str],
) -> None:
    """
    Update the Modalities field on a signature page.

    Args:
        page_id: Notion page ID of signature (with dashes)
        modalities: List of modality strings (e.g., ["Gene", "Protein", "Metabolite", "Lipid"])
    """
    cfg = get_config()

    # Map lowercase modality names to Notion select values
    modality_map = {
        "gene": "Gene",
        "protein": "Protein",
        "metabolite": "Metabolite",
        "lipid": "Lipid",
    }

    # Normalize modalities to Notion select values
    normalized_modalities = []
    for mod in modalities:
        mod_lower = mod.lower().strip()
        normalized = modality_map.get(mod_lower, mod.title())  # Default to title case
        normalized_modalities.append(normalized)

    try:
        # Fetch current page
        url = f"{cfg.notion.base_url}/pages/{page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()

        page = resp.json()
        props = page.get("properties", {}) or {}
        
        # Check if Modalities property exists
        if "Modalities" not in props:
            logger.debug(
                "[INGEST][SIGNATURES] Modalities property not found on signature page %s, "
                "skipping update",
                page_id,
            )
            return

        # Update Modalities field (merge with existing if any)
        existing_modalities = []
        modalities_prop = props.get("Modalities", {})
        if modalities_prop.get("type") == "multi_select":
            existing_modalities = [
                item.get("name", "") 
                for item in modalities_prop.get("multi_select", [])
            ]

        # Merge existing and new modalities
        all_modalities = list(set(existing_modalities + normalized_modalities))

        # Update page
        url = f"{cfg.notion.base_url}/pages/{page_id}"
        payload = {
            "properties": {
                "Modalities": {
                    "multi_select": [{"name": mod} for mod in all_modalities],
                },
            },
        }

        resp = requests.patch(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()

        logger.debug(
            "[INGEST][SIGNATURES] Updated Modalities for signature page %s: %s",
            page_id,
            ", ".join(all_modalities),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error updating Modalities for signature page %s: %r",
            page_id,
            e,
        )
        # Non-blocking - modalities update is optional

