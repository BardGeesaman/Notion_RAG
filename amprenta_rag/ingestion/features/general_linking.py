"""
General feature linking functions.

This module provides general-purpose functions for creating/finding feature pages
(lipids, metabolites, proteins, genes) and linking them to datasets.
"""

from __future__ import annotations

from typing import Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _find_or_create_feature_page(feature_type: str, feature_name: str) -> Optional[str]:
    """
    Find or create a feature page in the correct Notion DB.

    Feature Type → Notion DB Mapping:
    - lipid → Lipid Species DB
    - metabolite → Metabolite Features DB
    - protein → Protein Features DB
    - gene → Gene Features DB

    Args:
        feature_type: One of "lipid", "metabolite", "protein", "gene"
        feature_name: Normalized feature name

    Returns:
        Notion page ID (with dashes) or None if creation failed
    """
    cfg = get_config()

    # Map feature type to config attribute and DB name
    db_map = {
        "lipid": ("lipid_species_db_id", "Lipid Species"),
        "metabolite": ("metabolite_features_db_id", "Metabolite Features"),
        "protein": ("protein_features_db_id", "Protein Features"),
        "gene": ("gene_features_db_id", "Gene Features"),
    }

    if feature_type not in db_map:
        logger.warning(
            "[INGEST][FEATURE] Unknown feature_type '%s' for feature '%s'",
            feature_type,
            feature_name,
        )
        return None

    db_attr, db_name = db_map[feature_type]

    # Check if DB is configured
    if not hasattr(cfg.notion, db_attr) or not getattr(cfg.notion, db_attr):
        logger.debug(
            "[INGEST][FEATURE][WARN] Feature DB for %s not configured; skipping linking",
            feature_type,
        )
        return None

    db_id = getattr(cfg.notion, db_attr)

    # Special handling for lipid species (uses existing function)
    if feature_type == "lipid":
        from amprenta_rag.ingestion.signature_notion_crud import (
            find_or_create_lipid_species_page,
        )

        page_id = find_or_create_lipid_species_page(feature_name)
        if page_id:
            logger.debug(
                "[INGEST][FEATURE] Found/created lipid species page '%s' (id: %s)",
                feature_name,
                page_id,
            )
        return page_id

    # For other feature types, find or create page
    # First, try to find existing page
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "property": "Name",
                "title": {
                    "equals": feature_name,
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
                "[INGEST][FEATURE] Found existing %s feature page '%s' (id: %s)",
                feature_type,
                feature_name,
                page_id,
            )
            return page_id
    except Exception as e:
        logger.warning(
            "[INGEST][FEATURE] Error querying for existing %s feature '%s': %r",
            feature_type,
            feature_name,
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
                            "text": {"content": feature_name},
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
            "[INGEST][FEATURE] Created new %s feature page '%s' (id: %s)",
            feature_type,
            feature_name,
            page_id,
        )
        return page_id
    except Exception as e:
        logger.error(
            "[INGEST][FEATURE] Error creating %s feature page for '%s': %r",
            feature_type,
            feature_name,
            e,
        )
        return None


def _add_dataset_relation(
    feature_page_id: str, dataset_page_id: str, feature_type: str
) -> None:
    """
    Add dataset_page_id to the correct Feature → Dataset relation.

    Relation Type → Notion Relation Property:
    - Gene Features → "Transcriptomics Datasets" or "Datasets"
    - Protein Features → "Proteomics Datasets" or "Datasets"
    - Metabolite Features → "Metabolomics Datasets" or "Datasets"
    - Lipid Species → Try multiple candidates or skip if no suitable property exists

    Args:
        feature_page_id: Notion page ID of feature (with dashes)
        dataset_page_id: Notion page ID of dataset (with dashes)
        feature_type: One of "lipid", "metabolite", "protein", "gene"
    """
    cfg = get_config()

    # Map feature type to candidate relation property names (in order of preference)
    relation_property_candidates = {
        "gene": ["Transcriptomics Datasets", "Datasets", "Related Datasets"],
        "protein": ["Proteomics Datasets", "Datasets", "Related Datasets"],
        "metabolite": ["Metabolomics Datasets", "Datasets", "Related Datasets"],
        "lipid": [
            "Lipidomics Datasets",
            "Datasets",
            "Related Datasets",
            "Experimental Data Assets",
        ],
    }

    candidates = relation_property_candidates.get(feature_type, ["Datasets"])

    try:
        # Fetch current feature page to get existing relations and schema
        url = f"{cfg.notion.base_url}/pages/{feature_page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()

        page = resp.json()
        props = page.get("properties", {}) or {}

        # Find the first available relation property that matches our candidates
        relation_property = None
        current_rel = []

        # First, try candidate names
        for candidate in candidates:
            if candidate in props:
                prop_data = props.get(candidate, {})
                if prop_data.get("type") == "relation":
                    current_rel = prop_data.get("relation", []) or []
                    relation_property = candidate
                    break

        # If no candidate found, search for any relation property containing "dataset" (case-insensitive)
        if not relation_property:
            for prop_name, prop_data in props.items():
                if prop_data.get("type") == "relation":
                    prop_lower = prop_name.lower()
                    if "dataset" in prop_lower or "data asset" in prop_lower:
                        current_rel = prop_data.get("relation", []) or []
                        relation_property = prop_name
                        logger.info(
                            "[INGEST][FEATURE] Found dataset relation property '%s' for %s feature %s",
                            relation_property,
                            feature_type,
                            feature_page_id,
                        )
                        break

        # If still no suitable property found, log and return
        if not relation_property:
            available_relations = [
                name
                for name, data in props.items()
                if data.get("type") == "relation"
            ]
            logger.debug(
                "[INGEST][FEATURE] No suitable dataset relation property found for %s feature %s. "
                "Available relations: %s. Skipping relation update.",
                feature_type,
                feature_page_id,
                available_relations,
            )
            return

        existing_ids = {r.get("id", "") for r in current_rel if r.get("id")}

        # Check if already linked
        if dataset_page_id in existing_ids:
            logger.debug(
                "[INGEST][FEATURE] Dataset %s already linked to %s feature %s",
                dataset_page_id,
                feature_type,
                feature_page_id,
            )
            return  # Already linked

        # Add target page to relation
        updated_rel = [{"id": dataset_page_id}] + current_rel

        # Update feature page
        url = f"{cfg.notion.base_url}/pages/{feature_page_id}"
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

        logger.info(
            "[INGEST][FEATURE] Added dataset %s to %s feature %s (via property '%s')",
            dataset_page_id,
            feature_type,
            feature_page_id,
            relation_property,
        )
    except Exception as e:
        logger.warning(
            "[INGEST][FEATURE] Error adding dataset relation from %s feature %s to dataset %s: %r",
            feature_type,
            feature_page_id,
            dataset_page_id,
            e,
        )
        # Don't raise - relation updates are non-critical


def link_feature(
    feature_type: str,
    feature_name: str,
    dataset_page_id: str,
) -> None:
    """
    General-purpose feature linking helper.
    Creates or finds a feature page in the correct DB, then links dataset.

    Args:
        feature_type: One of "lipid", "metabolite", "protein", "gene"
        feature_name: Normalized feature name
        dataset_page_id: Notion page ID of dataset (with dashes)
    """
    try:
        feature_page_id = _find_or_create_feature_page(feature_type, feature_name)
        if feature_page_id:
            _add_dataset_relation(feature_page_id, dataset_page_id, feature_type)
            logger.debug(
                "[INGEST][FEATURE] Linked dataset %s to %s '%s'",
                dataset_page_id,
                feature_type,
                feature_name,
            )
    except Exception as e:
        logger.warning(
            "[INGEST][FEATURE] Error linking %s feature '%s' to dataset %s: %r",
            feature_type,
            feature_name,
            dataset_page_id,
            e,
        )
        # Non-blocking - continue with other features

