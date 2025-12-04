"""
Feature linking functions.

Handles creating/finding feature pages in Notion and linking them to
datasets, literature, experiments, and other items.
"""

from __future__ import annotations

from typing import List, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


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

