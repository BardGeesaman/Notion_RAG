# amprenta_rag/ingestion/feature_extraction.py

"""
Metabolite/Feature extraction and linking for ingestion pipelines.

This module provides functions to:
- Normalize metabolite names
- Extract metabolite features from mwTab data and text
- Link features to Notion items (datasets, literature, experiments, emails)
"""

from __future__ import annotations

import json
import re
from typing import Any, Dict, List, Optional, Set

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Common metabolite synonyms mapping (expandable)
METABOLITE_SYNONYMS: Dict[str, str] = {
    "l-glutamate": "glutamate",
    "l-glutamic acid": "glutamate",
    "d-glutamate": "glutamate",
    "glutamic acid": "glutamate",
    "l-glutamine": "glutamine",
    "d-glutamine": "glutamine",
    "gln": "glutamine",
    "glu": "glutamate",
}

# Common amino acids (for text scanning)
AMINO_ACIDS = [
    "alanine",
    "arginine",
    "asparagine",
    "aspartic acid",
    "cysteine",
    "glutamine",
    "glutamate",
    "glutamic acid",
    "glycine",
    "histidine",
    "isoleucine",
    "leucine",
    "lysine",
    "methionine",
    "phenylalanine",
    "proline",
    "serine",
    "threonine",
    "tryptophan",
    "tyrosine",
    "valine",
]

# Common nucleotides (for text scanning)
NUCLEOTIDES = [
    "adenosine",
    "adenine",
    "guanosine",
    "guanine",
    "cytidine",
    "cytosine",
    "thymidine",
    "thymine",
    "uracil",
    "uridine",
    "atp",
    "adp",
    "amp",
    "gtp",
    "gdp",
    "gmp",
    "ctp",
    "cdp",
    "cmp",
    "utp",
    "udp",
    "ump",
]


def normalize_metabolite_name(raw: str) -> str:
    """
    Normalize metabolite names from mwTab and CSV formats.

    Args:
        raw: Raw metabolite name (e.g., "HMDB:12345 Glutamine", "L-glutamate")

    Returns:
        Canonical normalized name (e.g., "glutamine", "glutamate")
    """
    if not raw:
        return ""

    # Lowercase and strip
    normalized = raw.lower().strip()

    # Remove prefixes like "HMDB:", "KEGG:", "CHEBI:", etc.
    normalized = re.sub(
        r"^(hmdb|kegg|chebi|pubchem|cas)[:\s]+", "", normalized, flags=re.IGNORECASE
    )

    # Remove common prefixes/suffixes
    normalized = re.sub(r"^\d+\s*[-:]?\s*", "", normalized)  # Leading numbers
    normalized = re.sub(r"\s*\(.*?\)\s*$", "", normalized)  # Trailing parentheses

    # Strip again after cleanup
    normalized = normalized.strip()

    # Apply synonym mapping
    normalized_lower = normalized.lower()
    if normalized_lower in METABOLITE_SYNONYMS:
        normalized = METABOLITE_SYNONYMS[normalized_lower]

    # Capitalize first letter for canonical form
    if normalized:
        normalized = (
            normalized[0].upper() + normalized[1:]
            if len(normalized) > 1
            else normalized.upper()
        )

    return normalized


def extract_features_from_mwtab(mwtab_json: Dict[str, Any]) -> List[str]:
    """
    Extract metabolite names from mwTab JSON data.

    Args:
        mwtab_json: Parsed mwTab JSON dictionary

    Returns:
        List of normalized metabolite names
    """
    metabolite_names: Set[str] = set()

    # Check for MS_METABOLITE_DATA section
    metabolite_sections = [
        "MS_METABOLITE_DATA",
        "GC_METABOLITE_DATA",
        "LC_METABOLITE_DATA",
        "METABOLITE_DATA",
    ]

    for section_key in metabolite_sections:
        if section_key not in mwtab_json:
            continue

        section = mwtab_json[section_key]
        if not isinstance(section, dict):
            continue

        # Extract data array
        data_array = section.get("Data", [])
        if not isinstance(data_array, list):
            continue

        for row in data_array:
            if isinstance(row, dict):
                # Look for Metabolite key (case-insensitive)
                metabolite_key = None
                for key in row.keys():
                    if key.lower() in ["metabolite", "metabolite_name", "compound"]:
                        metabolite_key = key
                        break

                if metabolite_key:
                    metabolite_raw = row.get(metabolite_key)
                    if metabolite_raw and isinstance(metabolite_raw, str):
                        normalized = normalize_metabolite_name(metabolite_raw)
                        if normalized:
                            metabolite_names.add(normalized)

            elif isinstance(row, list) and len(row) > 0:
                # First element might be metabolite name
                if isinstance(row[0], str):
                    normalized = normalize_metabolite_name(row[0])
                    if normalized:
                        metabolite_names.add(normalized)

    return sorted(list(metabolite_names))


def extract_features_from_text(text: str) -> List[str]:
    """
    Extract metabolite names from text content using pattern matching.

    This is a lightweight scanner for literature chunks, email/note chunks,
    and experiment descriptions.

    Args:
        text: Text content to scan

    Returns:
        List of normalized metabolite names found in text
    """
    if not text:
        return []

    text_lower = text.lower()
    found_metabolites: Set[str] = set()

    # Scan for common amino acids
    for aa in AMINO_ACIDS:
        # Look for whole word matches (not substring)
        pattern = r"\b" + re.escape(aa) + r"\b"
        if re.search(pattern, text_lower, re.IGNORECASE):
            found_metabolites.add(normalize_metabolite_name(aa))

    # Scan for nucleotides
    for nt in NUCLEOTIDES:
        pattern = r"\b" + re.escape(nt) + r"\b"
        if re.search(pattern, text_lower, re.IGNORECASE):
            found_metabolites.add(normalize_metabolite_name(nt))

    # Scan for ceramide patterns (basic)
    ceramide_patterns = [
        r"\bcer\s*\([^)]+\)",
        r"\bceramide\s*\([^)]+\)",
        r"\bcer\([^)]+\)",
    ]
    for pattern in ceramide_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            found_metabolites.add(normalize_metabolite_name(match.group(0)))

    # Scan for common small molecules mentioned in signatures
    common_metabolites = [
        "glucose",
        "lactate",
        "pyruvate",
        "citrate",
        "succinate",
        "fumarate",
        "malate",
        "oxaloacetate",
        "alpha-ketoglutarate",
        "acetylate",
        "carnitine",
        "creatine",
        "choline",
    ]
    for metabolite in common_metabolites:
        pattern = r"\b" + re.escape(metabolite) + r"\b"
        if re.search(pattern, text_lower, re.IGNORECASE):
            found_metabolites.add(normalize_metabolite_name(metabolite))

    return sorted(list(found_metabolites))


def _find_or_create_metabolite_page(metabolite_name: str) -> Optional[str]:
    """
    Find or create a Metabolite Features page in Notion.

    Args:
        metabolite_name: Canonical metabolite name

    Returns:
        Notion page ID (with dashes) or None if creation failed
    """
    cfg = get_config()

    # Check if metabolite_features_db_id is configured
    if (
        not hasattr(cfg.notion, "metabolite_features_db_id")
        or not cfg.notion.metabolite_features_db_id
    ):
        logger.warning(
            "[INGEST][FEATURES] Metabolite Features database ID not configured. "
            "Set NOTION_METABOLITE_FEATURES_DB_ID in config."
        )
        return None

    db_id = cfg.notion.metabolite_features_db_id

    # First, try to find existing page
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "property": "Name",
                "title": {
                    "equals": metabolite_name,
                },
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
                "[INGEST][FEATURES] Found existing metabolite page for %s: %s",
                metabolite_name,
                page_id,
            )
            return page_id
    except Exception as e:
        logger.warning(
            "[INGEST][FEATURES] Error querying for existing metabolite %s: %r",
            metabolite_name,
            e,
        )

    # Create new page
    try:
        url = f"{cfg.notion.base_url}/pages"
        payload = {
            "parent": {"database_id": db_id},
            "properties": {
                "Name": {
                    "title": [
                        {
                            "text": {"content": metabolite_name},
                        }
                    ]
                },
            },
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
            "[INGEST][FEATURES] Created new metabolite page for %s: %s",
            metabolite_name,
            page_id,
        )
        return page_id
    except Exception as e:
        logger.error(
            "[INGEST][FEATURES] Error creating metabolite page for %s: %r",
            metabolite_name,
            e,
        )
        return None


def _add_relation_to_metabolite_page(
    metabolite_page_id: str,
    target_page_id: str,
    relation_property: str,
) -> None:
    """
    Add a relation from a metabolite page to a target page.

    Args:
        metabolite_page_id: Notion page ID of metabolite (with dashes)
        target_page_id: Notion page ID to link to (with dashes)
        relation_property: Name of relation property on metabolite page
    """
    cfg = get_config()

    try:
        # Fetch current metabolite page to get existing relations
        url = f"{cfg.notion.base_url}/pages/{metabolite_page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()

        page = resp.json()
        props = page.get("properties", {}) or {}
        current_rel = props.get(relation_property, {}).get("relation", []) or []
        existing_ids = {r.get("id", "") for r in current_rel if r.get("id")}

        # Check if already linked
        if target_page_id in existing_ids:
            return  # Already linked

        # Add target page to relation
        updated_rel = [{"id": target_page_id}] + current_rel

        # Update metabolite page
        url = f"{cfg.notion.base_url}/pages/{metabolite_page_id}"
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
            "[INGEST][FEATURES] Added relation %s from metabolite %s to target %s",
            relation_property,
            metabolite_page_id,
            target_page_id,
        )
    except Exception as e:
        logger.warning(
            "[INGEST][FEATURES] Error adding relation %s from metabolite %s to target %s: %r",
            relation_property,
            metabolite_page_id,
            target_page_id,
            e,
        )
        # Don't raise - relation updates are non-critical


def link_features_to_notion_items(
    feature_names: List[str],
    item_page_id: str,
    item_type: str,
) -> None:
    """
    Link metabolite features to a Notion item (dataset, literature, email, experiment).

    For each feature:
    1. Creates or retrieves the corresponding Metabolite Features page
    2. Adds the correct relation to the target item

    Args:
        feature_names: List of normalized metabolite names
        item_page_id: Notion page ID of the target item (with dashes)
        item_type: One of "dataset", "literature", "email", "experiment"
    """
    if not feature_names:
        logger.debug(
            "[INGEST][FEATURES] No features to link for %s %s",
            item_type,
            item_page_id,
        )
        return

    # Map item_type to relation property name on Metabolite Features page
    relation_map = {
        "dataset": "Datasets",
        "literature": "Literature Mentions",
        "email": "Emails / Notes",
        "experiment": "Experiments",
    }

    relation_property = relation_map.get(item_type)
    if not relation_property:
        logger.warning(
            "[INGEST][FEATURES] Unknown item_type %s for page %s",
            item_type,
            item_page_id,
        )
        return

    logger.info(
        "[INGEST][FEATURES] Linking %d feature(s) to %s page %s",
        len(feature_names),
        item_type,
        item_page_id,
    )

    linked_count = 0
    for feature_name in feature_names:
        try:
            metabolite_page_id = _find_or_create_metabolite_page(feature_name)
            if metabolite_page_id:
                _add_relation_to_metabolite_page(
                    metabolite_page_id,
                    item_page_id,
                    relation_property,
                )
                linked_count += 1
        except Exception as e:
            logger.warning(
                "[INGEST][FEATURES] Error linking feature %s to %s %s: %r",
                feature_name,
                item_type,
                item_page_id,
                e,
            )

    logger.info(
        "[INGEST][FEATURES] Successfully linked %d/%d features to %s page %s",
        linked_count,
        len(feature_names),
        item_type,
        item_page_id,
    )
