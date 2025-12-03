"""
Notion CRUD operations for lipid signatures, components, and species.

This module handles all Notion database operations for:
- Lipid Signatures database
- Lipid Signature Components database
- Lipid Species database

All functions are idempotent - they find existing pages or create new ones.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import (Signature,
                                                      SignatureComponent)
from amprenta_rag.signatures.species_matching import (classify_lipid_class,
                                                     normalize_species_name)

logger = get_logger(__name__)

__all__ = [
    "generate_signature_short_id",
    "find_or_create_signature_page",
    "update_signature_page_if_needed",
    "find_or_create_component_page",
    "find_or_create_lipid_species_page",
    "update_lipid_species_synonyms",
]


def generate_signature_short_id(signature_name: str, version: Optional[str] = None) -> str:
    """
    Generate a deterministic Short ID for a signature.

    Args:
        signature_name: Signature name
        version: Optional version string

    Returns:
        Deterministic Short ID (e.g., "ALS-CSF-Core-6Cer-v1")
    """
    # Normalize name for ID generation
    normalized = re.sub(r"[^a-zA-Z0-9]+", "-", signature_name).strip("-")
    normalized = normalized[:50]  # Limit length

    if version:
        return f"{normalized}-v{version}"

    # Generate version-less ID
    return normalized


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


def find_or_create_component_page(
    component: SignatureComponent,
    signature_page_id: str,
    disease_context: Optional[List[str]] = None,
    matrix: Optional[List[str]] = None,
) -> Optional[str]:
    """
    Find or create a Lipid Signature Component page in Notion.

    Args:
        component: SignatureComponent object
        signature_page_id: Notion page ID of parent signature (with dashes)
        disease_context: Optional list of disease context strings
        matrix: Optional list of matrix strings

    Returns:
        Notion page ID (with dashes) or None if creation failed
    """
    cfg = get_config()

    # Check if component DB ID is configured
    if (
        not hasattr(cfg.notion, "signature_component_db_id")
        or not cfg.notion.signature_component_db_id
    ):
        logger.warning(
            "[INGEST][SIGNATURES] Lipid Signature Components database ID not configured. "
            "Set NOTION_SIGNATURE_COMPONENT_DB_ID in config."
        )
        return None

    db_id = cfg.notion.signature_component_db_id

    # Map direction to Notion select values
    direction_map = {
        "↑": "Up",
        "↓": "Down",
        "neutral": "NoChange",
        "complex": "Complex",
        None: "Unknown",
    }
    direction_value = direction_map.get(component.direction, "Unknown")

    # Try to find existing component by Component Name and Signature relation
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "and": [
                    {
                        "property": "Component Name",
                        "title": {"equals": component.species},
                    },
                    {
                        "property": "Signature",
                        "relation": {"contains": signature_page_id},
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
                "[INGEST][SIGNATURES] Found existing component page for %s: %s",
                component.species,
                page_id,
            )
            return page_id
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error querying for existing component %s: %r",
            component.species,
            e,
        )

    # Create new component page
    try:
        url = f"{cfg.notion.base_url}/pages"

        properties: Dict[str, Any] = {
            "Component Name": {
                "title": [{"text": {"content": component.species}}],
            },
            "Raw Name": {
                "rich_text": [{"text": {"content": component.species}}],
            },
            "Direction": {
                "select": {"name": direction_value},
            },
            "Signature": {
                "relation": [{"id": signature_page_id}],
            },
        }

        if component.weight is not None:
            properties["Weight"] = {
                "number": component.weight,
            }

        if disease_context:
            properties["Disease Context"] = {
                "multi_select": [{"name": d} for d in disease_context],
            }

        if matrix:
            properties["Matrix"] = {
                "multi_select": [{"name": m} for m in matrix],
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
            "[INGEST][SIGNATURES] Created new component page for %s: %s",
            component.species,
            page_id,
        )
        return page_id
    except Exception as e:
        logger.error(
            "[INGEST][SIGNATURES] Error creating component page for %s: %r",
            component.species,
            e,
        )
        return None


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

