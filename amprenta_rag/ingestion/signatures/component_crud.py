"""
Component page CRUD operations.

Handles creating and finding signature component pages in Notion.
Supports multi-omics components (genes, proteins, metabolites, lipids).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import SignatureComponent

logger = get_logger(__name__)


def find_or_create_component_page(
    component: SignatureComponent,
    signature_page_id: str,
    disease_context: Optional[List[str]] = None,
    matrix: Optional[List[str]] = None,
) -> Optional[str]:
    """
    Find or create a multi-omics Signature Component page in Notion.

    Supports components for genes, proteins, metabolites, and lipids.

    Args:
        component: SignatureComponent object (with feature_name and feature_type)
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
            "[INGEST][SIGNATURES] Signature Components database ID not configured. "
            "Set NOTION_SIGNATURE_COMPONENT_DB_ID in config."
        )
        return None

    db_id = cfg.notion.signature_component_db_id

    # Get feature name (use feature_name, fallback to species for backward compatibility)
    feature_name = getattr(component, "feature_name", component.species)
    feature_type = getattr(component, "feature_type", "lipid")  # Default to lipid for backward compat

    # Map direction to Notion select values
    direction_map = {
        "↑": "Up",
        "↓": "Down",
        "neutral": "NoChange",
        "complex": "Complex",
        None: "Unknown",
    }
    direction_value = direction_map.get(component.direction, "Unknown")

    # Map feature_type to Notion select values
    feature_type_map = {
        "gene": "Gene",
        "protein": "Protein",
        "metabolite": "Metabolite",
        "lipid": "Lipid",
    }
    feature_type_value = feature_type_map.get(feature_type, "Lipid")  # Default to Lipid

    # Try to find existing component by Component Name and Signature relation
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "and": [
                    {
                        "property": "Component Name",
                        "title": {"equals": feature_name},
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
                "[INGEST][SIGNATURES] Found existing component page for %s (%s): %s",
                feature_name,
                feature_type,
                page_id,
            )
            return page_id
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error querying for existing component %s: %r",
            feature_name,
            e,
        )

    # Create new component page
    try:
        url = f"{cfg.notion.base_url}/pages"

        properties: Dict[str, Any] = {
            "Component Name": {
                "title": [{"text": {"content": feature_name}}],
            },
            "Raw Name": {
                "rich_text": [{"text": {"content": feature_name}}],
            },
            "Direction": {
                "select": {"name": direction_value},
            },
            "Signature": {
                "relation": [{"id": signature_page_id}],
            },
        }

        # Add Feature Type if property exists (gracefully skip if not)
        try:
            properties["Feature Type"] = {
                "select": {"name": feature_type_value},
            }
        except Exception:
            # Feature Type property may not exist yet - log and continue
            logger.debug(
                "[INGEST][SIGNATURES] Feature Type property not available for component %s",
                feature_name,
            )

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
            "[INGEST][SIGNATURES] Created new component page for %s (%s): %s",
            feature_name,
            feature_type,
            page_id,
        )
        return page_id
    except Exception as e:
        logger.error(
            "[INGEST][SIGNATURES] Error creating component page for %s: %r",
            feature_name,
            e,
        )
        # If Feature Type property caused the error, retry without it
        if "Feature Type" in str(e) or "feature type" in str(e).lower():
            logger.debug(
                "[INGEST][SIGNATURES] Retrying component creation without Feature Type property"
            )
            try:
                # Remove Feature Type from properties
                if "Feature Type" in properties:
                    del properties["Feature Type"]
                
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
                    "[INGEST][SIGNATURES] Created component page (without Feature Type) for %s: %s",
                    feature_name,
                    page_id,
                )
                return page_id
            except Exception as e2:
                logger.error(
                    "[INGEST][SIGNATURES] Error creating component page (retry): %r",
                    e2,
                )
        return None

