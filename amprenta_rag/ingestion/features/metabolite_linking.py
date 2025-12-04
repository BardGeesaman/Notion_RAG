"""
Metabolite-specific linking functions.

This module provides functions for creating/finding metabolite feature pages
and linking them to Notion items (datasets, literature, emails, experiments).
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

