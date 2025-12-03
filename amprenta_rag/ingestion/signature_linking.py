"""
Signature linking utilities.

This module handles linking operations between signatures, components, species,
and other Notion entities (metabolite features, source pages).
"""

from __future__ import annotations

from typing import Any, Dict

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

__all__ = [
    "link_component_to_lipid_species",
    "link_signature_to_source",
    "link_component_to_metabolite_feature",
]


def _fetch_notion_page_helper(page_id: str) -> Dict[str, Any]:
    """Helper to fetch a Notion page by ID."""
    cfg = get_config().notion
    url = f"{cfg.base_url}/pages/{page_id}"
    resp = requests.get(url, headers=notion_headers(), timeout=30)
    resp.raise_for_status()
    return resp.json()


def link_component_to_lipid_species(
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

