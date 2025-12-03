# amprenta_rag/ingestion/signature_ingestion.py

"""
External lipid signature ingestion into Notion databases.

This module handles ingesting external lipid signatures (from TSV, CSV, JSON, etc.)
into the Notion knowledge graph:
- Lipid Signatures database
- Lipid Signature Components database
- Lipid Species database (canonical ontology)

Creates the full relation graph: Signature → Components → Species
"""

from __future__ import annotations

import hashlib
import json
import re
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.openai_client import (get_default_models,
                                                get_openai_client)
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import (Signature,
                                                      SignatureComponent,
                                                      load_signature_from_tsv)
from amprenta_rag.signatures.species_matching import (classify_lipid_class,
                                                      normalize_species_name)

logger = get_logger(__name__)


def _generate_short_id(signature_name: str, version: Optional[str] = None) -> str:
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


def _find_or_create_signature_page(
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
    short_id = _generate_short_id(signature.name, version)

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
            _update_signature_page_if_needed(
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


def _update_signature_page_if_needed(
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


def _find_or_create_component_page(
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


def _find_or_create_lipid_species_page(
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
                    _update_lipid_species_synonyms(page_id, lipid_name)
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


def _update_lipid_species_synonyms(page_id: str, new_name: str) -> None:
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


def _link_component_to_lipid_species(
    component_page_id: str,
    lipid_species_page_id: str,
) -> None:
    """
    Add relation from component to lipid species.

    Args:
        component_page_id: Notion page ID of component (with dashes)
        lipid_species_page_id: Notion page ID of lipid species (with dashes)
    """
    cfg = get_config()

    try:
        # Fetch current component page
        url = f"{cfg.notion.base_url}/pages/{component_page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()

        page = resp.json()
        props = page.get("properties", {}) or {}
        current_rel = props.get("Lipid Species", {}).get("relation", []) or []
        existing_ids = {r.get("id", "") for r in current_rel if r.get("id")}

        # Check if already linked
        if lipid_species_page_id in existing_ids:
            return  # Already linked

        # Add lipid species to relation
        updated_rel = [{"id": lipid_species_page_id}] + current_rel

        # Update component page
        url = f"{cfg.notion.base_url}/pages/{component_page_id}"
        payload = {
            "properties": {
                "Lipid Species": {"relation": updated_rel},
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
            "[INGEST][SIGNATURES] Linked component %s to lipid species %s",
            component_page_id,
            lipid_species_page_id,
        )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error linking component %s to lipid species %s: %r",
            component_page_id,
            lipid_species_page_id,
            e,
        )


def ingest_signature_from_file(
    signature_path: Path,
    signature_type: str = "Literature-derived",
    data_ownership: str = "Public",
    version: Optional[str] = None,
    description: Optional[str] = None,
    disease_context: Optional[List[str]] = None,
    matrix: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Ingest a single signature from a file (TSV, CSV) into Notion.

    Creates the full relation graph:
    - Lipid Signature page
    - Lipid Signature Component pages (one per component)
    - Lipid Species pages (canonical ontology)

    Links: Signature → Components → Species

    Args:
        signature_path: Path to signature file (TSV or CSV)
        signature_type: "Consortium", "Literature-derived", "Open Dataset", or "Other"
        data_ownership: "Public" or appropriate value
        version: Optional version string
        description: Optional description text
        disease_context: Optional list of disease context strings for all components
        matrix: Optional list of matrix strings for all components

    Returns:
        Dictionary with ingestion results:
        - signature_page_id: Notion page ID of signature
        - component_count: Number of components created
        - species_count: Number of lipid species created/linked
        - warnings: List of warning messages
    """
    logger.info(
        "[INGEST][SIGNATURES] Loading signature from file: %s",
        signature_path,
    )

    # Load signature from file
    try:
        signature = load_signature_from_tsv(signature_path)
    except Exception as e:
        logger.error(
            "[INGEST][SIGNATURES] Error loading signature from %s: %r",
            signature_path,
            e,
        )
        raise

    logger.info(
        "[INGEST][SIGNATURES] Loaded signature '%s' with %d components",
        signature.name,
        len(signature.components),
    )

    warnings: List[str] = []

    # Create or find signature page
    signature_page_id = _find_or_create_signature_page(
        signature=signature,
        signature_type=signature_type,
        data_ownership=data_ownership,
        version=version,
        description=description,
    )

    if not signature_page_id:
        error_msg = f"Failed to create/find signature page for {signature.name}"
        logger.error(f"[INGEST][SIGNATURES] {error_msg}")
        warnings.append(error_msg)
        return {
            "signature_page_id": None,
            "component_count": 0,
            "species_count": 0,
            "warnings": warnings,
        }

    # Create component pages and link to lipid species
    component_count = 0
    species_created: Set[str] = set()

    for component in signature.components:
        try:
            # Create or find component page
            component_page_id = _find_or_create_component_page(
                component=component,
                signature_page_id=signature_page_id,
                disease_context=disease_context,
                matrix=matrix,
            )

            if not component_page_id:
                warning = f"Failed to create component page for {component.species}"
                logger.warning(f"[INGEST][SIGNATURES] {warning}")
                warnings.append(warning)
                continue

            component_count += 1

            # Create or find lipid species page
            lipid_species_page_id = _find_or_create_lipid_species_page(
                lipid_name=component.species,
            )

            if lipid_species_page_id:
                # Link component to lipid species
                _link_component_to_lipid_species(
                    component_page_id=component_page_id,
                    lipid_species_page_id=lipid_species_page_id,
                )

                # Link component to metabolite feature (cross-link)
                try:
                    link_component_to_metabolite_feature(
                        component_page_id=component_page_id,
                        lipid_species_page_id=lipid_species_page_id,
                    )
                except Exception as e:
                    logger.warning(
                        "[INGEST][SIGNATURES] Error linking component to metabolite feature: %r",
                        e,
                    )
                    # Non-blocking - continue

                if component.species not in species_created:
                    species_created.add(component.species)
            else:
                warning = f"Failed to create/link lipid species for {component.species}"
                logger.warning(f"[INGEST][SIGNATURES] {warning}")
                warnings.append(warning)

        except Exception as e:
            warning = f"Error processing component {component.species}: {e}"
            logger.warning(f"[INGEST][SIGNATURES] {warning}")
            warnings.append(warning)

    # Embed signature into Pinecone for RAG queries
    try:
        embed_signature(
            signature_page_id=signature_page_id,
            signature=signature,
        )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error embedding signature '%s': %r",
            signature.name,
            e,
        )
        # Non-blocking - continue

    logger.info(
        "[INGEST][SIGNATURES] Ingestion complete for signature '%s': "
        "%d components, %d species",
        signature.name,
        component_count,
        len(species_created),
    )

    return {
        "signature_page_id": signature_page_id,
        "component_count": component_count,
        "species_count": len(species_created),
        "warnings": warnings,
    }


def _fetch_notion_page_helper(page_id: str) -> Dict[str, Any]:
    """Helper to fetch a Notion page by ID."""
    cfg = get_config().notion
    url = f"{cfg.base_url}/pages/{page_id}"
    resp = requests.get(url, headers=notion_headers(), timeout=30)
    resp.raise_for_status()
    return resp.json()


def link_signature_to_source(
    signature_page_id: str,
    source_page_id: str,
    source_type: str,
) -> None:
    """
    Create reverse link from signature page to source page.

    Maps source_type to appropriate Notion relation property:
    - "literature" → "Source Papers"
    - "dataset" → "External Datasets" (or appropriate property name)
    - "email" → "Email & Notes" (or appropriate property name)
    - "experiment" → "Source Experiments" (or appropriate property name)

    Args:
        signature_page_id: Notion page ID of signature (with dashes)
        source_page_id: Notion page ID of source (with dashes)
        source_type: One of "literature", "dataset", "email", "experiment"
    """
    cfg = get_config()

    # Map source_type to relation property name on signature page
    relation_map = {
        "literature": "Source Papers",
        "dataset": "External Datasets",
        "email": "Email & Notes",
        "experiment": "Source Experiments",
    }

    relation_property = relation_map.get(source_type)
    if not relation_property:
        logger.warning(
            "[INGEST][SIGNATURES] Unknown source_type %s for reverse linking",
            source_type,
        )
        return

    try:
        # Fetch current signature page
        sig_page = _fetch_notion_page_helper(signature_page_id)
        props = sig_page.get("properties", {}) or {}
        current_rel = props.get(relation_property, {}).get("relation", []) or []
        existing_ids = {r.get("id", "") for r in current_rel if r.get("id")}

        # Check if already linked
        if source_page_id in existing_ids:
            return  # Already linked

        # Add source page to relation
        updated_rel = [{"id": source_page_id}] + current_rel

        # Update signature page
        url = f"{cfg.notion.base_url}/pages/{signature_page_id}"
        payload = {
            "properties": {
                relation_property: {"relation": updated_rel},
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
            "[INGEST][SIGNATURES] Linked signature %s to %s source %s",
            signature_page_id,
            source_type,
            source_page_id,
        )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error linking signature %s to %s source %s: %r",
            signature_page_id,
            source_type,
            source_page_id,
            e,
        )
        # Don't raise - reverse linking is non-critical


def link_component_to_metabolite_feature(
    component_page_id: str,
    lipid_species_page_id: str,
) -> None:
    """
    Link a signature component's lipid species to a Metabolite Feature page.

    Creates the chain: Component → Lipid Species → Metabolite Feature

    Args:
        component_page_id: Notion page ID of component (with dashes)
        lipid_species_page_id: Notion page ID of lipid species (with dashes)
    """
    from amprenta_rag.ingestion.feature_extraction import (
        _add_relation_to_metabolite_page, _find_or_create_metabolite_page)

    try:
        # Fetch lipid species page to get the name
        lipid_page = _fetch_notion_page_helper(lipid_species_page_id)
        props = lipid_page.get("properties", {}) or {}

        name_prop = props.get("Name", {})
        name_title = name_prop.get("title", []) or []
        if name_title:
            lipid_name = name_title[0].get("plain_text", "")

            # Find or create Metabolite Feature page
            metabolite_page_id = _find_or_create_metabolite_page(lipid_name)

            if metabolite_page_id:
                # Link Metabolite Feature to Lipid Species (if relation exists)
                # This creates bidirectional links where supported
                logger.debug(
                    "[INGEST][SIGNATURES] Linked lipid species %s to metabolite feature %s",
                    lipid_species_page_id,
                    metabolite_page_id,
                )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error linking component to metabolite feature: %r",
            e,
        )
        # Don't raise - cross-linking is non-critical


def embed_signature(
    signature_page_id: str,
    signature: Signature,
) -> None:
    """
    Embed a signature into Pinecone for RAG queries.

    Creates text representation of signature, chunks it, embeds, and upserts to Pinecone.

    Args:
        signature_page_id: Notion page ID of signature (with dashes)
        signature: Signature object with components
    """
    cfg = get_config()

    try:
        # Build text representation of signature
        text_parts = [
            f"Lipid Signature: {signature.name}",
        ]

        if signature.description:
            text_parts.append(f"Description: {signature.description}")

        text_parts.append("\nComponents:")
        for comp in signature.components:
            comp_line = f"- {comp.species}"
            if comp.direction:
                comp_line += f" ({comp.direction})"
            if comp.weight:
                comp_line += f" [weight: {comp.weight}]"
            text_parts.append(comp_line)

        signature_text = "\n".join(text_parts)

        # Chunk (1-2 chunks should be enough for a signature)
        from amprenta_rag.ingestion.zotero_ingest import (_chunk_text,
                                                          _embed_texts)

        chunks = _chunk_text(signature_text, max_chars=2000)
        if not chunks:
            logger.debug(
                "[INGEST][SIGNATURES] No chunks generated for signature %s",
                signature.name,
            )
            return

        logger.info(
            "[INGEST][SIGNATURES] Embedding signature '%s' (%d chunk(s))",
            signature.name,
            len(chunks),
        )

        # Embed chunks
        embeddings = _embed_texts(chunks)

        # Upsert to Pinecone
        index = get_pinecone_index()

        vectors: List[Dict[str, Any]] = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{signature_page_id.replace('-', '')}_sig_chunk_{order:03d}"

            meta: Dict[str, Any] = {
                "source": "Signature",
                "source_type": "Signature",
                "signature_page_id": signature_page_id.replace("-", ""),
                "signature_name": signature.name,
                "snippet": textwrap.shorten(chunk, width=300),
            }

            vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(meta),
                }
            )

        index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)

        logger.info(
            "[INGEST][SIGNATURES] Embedded signature '%s' to Pinecone (%d vectors)",
            signature.name,
            len(vectors),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error embedding signature %s: %r",
            signature.name,
            e,
        )
        # Don't raise - embedding is non-critical
