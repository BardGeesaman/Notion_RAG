"""
Lipid species page CRUD operations.

Handles creating, finding, and updating lipid species pages in Notion.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.species_matching import (
    classify_lipid_class,
    normalize_species_name,
)

logger = get_logger(__name__)


def find_or_create_lipid_species_page(
    lipid_name: str,
) -> Optional[str]:
    """
    Find or create a Lipid Species page in Notion.

    Args:
        lipid_name: Normalized or original lipid species name

    Returns:
        Notion page ID (with dashes) or None if creation failed
    """

    cfg = get_config()

    # Check if lipid species DB ID is configured
    if (
        not hasattr(cfg.notion, "lipid_species_db_id")
        or not cfg.notion.lipid_species_db_id
    ):
        logger.warning(
            "[INGEST][SIGNATURES] Lipid Species database ID not configured. "
            "Set NOTION_LIPID_SPECIES_DB_ID in config."
        )
        return None

    db_id = cfg.notion.lipid_species_db_id

    # Normalize name for matching
    normalized = normalize_species_name(lipid_name)

    # Classify lipid class
    lipid_class = classify_lipid_class(lipid_name)

    # Try to find existing page by name (normalized matching)
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "or": [
                    {
                        "property": "Name",
                        "title": {"equals": lipid_name},
                    },
                    {
                        "property": "Synonyms",
                        "rich_text": {"contains": normalized},
                    },
                ],
            },
            "page_size": 10,  # Get multiple to check normalization
        }

        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()

        results = resp.json().get("results", [])

        # Check if any result matches normalized name
        for result in results:
            props = result.get("properties", {}) or {}
            name_prop = props.get("Name", {})
            name_title = name_prop.get("title", []) or []
            if name_title:
                existing_name = name_title[0].get("plain_text", "")
                if normalize_species_name(existing_name) == normalized:
                    page_id = result.get("id", "")
                    logger.debug(
                        "[INGEST][SIGNATURES] Found existing lipid species page for %s: %s",
                        lipid_name,
                        page_id,
                    )
                    # Update synonyms if needed
                    update_lipid_species_synonyms(page_id, lipid_name)
                    return page_id
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error querying for existing lipid species %s: %r",
            lipid_name,
            e,
        )

    # Create new lipid species page
    try:
        url = f"{cfg.notion.base_url}/pages"

        properties: Dict[str, Any] = {
            "Name": {
                "title": [{"text": {"content": lipid_name}}],
            },
        }

        if lipid_class:
            properties["Class"] = {
                "select": {"name": lipid_class},
            }

        # Add original name as synonym
        synonyms_text = f"{lipid_name}\n{normalized}"
        properties["Synonyms"] = {
            "rich_text": [{"text": {"content": synonyms_text}}],
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
            "[INGEST][SIGNATURES] Created new lipid species page for %s: %s (Class: %s)",
            lipid_name,
            page_id,
            lipid_class or "Unknown",
        )
        return page_id
    except Exception as e:
        logger.error(
            "[INGEST][SIGNATURES] Error creating lipid species page for %s: %r",
            lipid_name,
            e,
        )
        return None


def update_lipid_species_synonyms(page_id: str, new_name: str) -> None:
    """Update synonyms list for an existing lipid species page."""
    cfg = get_config()

    try:
        # Fetch current page
        url = f"{cfg.notion.base_url}/pages/{page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()

        page = resp.json()
        props = page.get("properties", {}) or {}

        synonyms_prop = props.get("Synonyms", {})
        synonyms_text = ""
        if "rich_text" in synonyms_prop:
            rich_text = synonyms_prop.get("rich_text", []) or []
            synonyms_text = "".join(rt.get("plain_text", "") for rt in rich_text)

        # Check if new name is already in synonyms
        if new_name.lower() not in synonyms_text.lower():
            # Add new name to synonyms
            updated_synonyms = f"{synonyms_text}\n{new_name}".strip()

            url = f"{cfg.notion.base_url}/pages/{page_id}"
            payload = {
                "properties": {
                    "Synonyms": {
                        "rich_text": [{"text": {"content": updated_synonyms}}],
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
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error updating synonyms for lipid species %s: %r",
            page_id,
            e,
        )

