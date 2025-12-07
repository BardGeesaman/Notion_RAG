"""
Signature loading from Notion.

Handles fetching signatures from Notion and loading them into Signature objects.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import (
    Signature,
    SignatureComponent,
)

logger = get_logger(__name__)

def notion_headers() -> Dict[str, str]:
    """DEPRECATED: Notion support removed. Returns empty headers dict."""
    logger.debug("[SIGNATURE-MATCHING][LOADER] notion_headers() deprecated - Notion support removed")
    return {}


def fetch_all_signatures_from_notion() -> List[Dict[str, Any]]:
    """
    Fetch all signature pages from Notion Lipid Signatures database.

    Returns:
        List of signature page dictionaries
    """
    cfg = get_config()

    if not cfg.notion.signature_db_id:
        logger.warning("[INGEST][SIGNATURE-MATCH] Signature database ID not configured")
        return []

    try:
        url = f"{cfg.notion.base_url}/databases/{cfg.notion.signature_db_id}/query"
        all_pages = []
        has_more = True
        start_cursor = None

        while has_more:
            payload = {
                "page_size": 100,
            }
            if start_cursor:
                payload["start_cursor"] = start_cursor

            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()

            data = resp.json()
            all_pages.extend(data.get("results", []))
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")

        logger.debug(
            "[INGEST][SIGNATURE-MATCH] Fetched %d signatures from Notion",
            len(all_pages),
        )
        return all_pages

    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Error fetching signatures from Notion: %r",
            e,
        )
        return []


def load_signature_from_notion_page(
    signature_page: Dict[str, Any],
) -> Optional[Signature]:
    """
    Load a Signature object from a Notion page.

    Fetches components from the signature page and builds a Signature object.

    Args:
        signature_page: Notion page dictionary for a signature

    Returns:
        Signature object or None if loading failed
    """
    cfg = get_config()
    signature_page_id = signature_page.get("id", "")

    if not signature_page_id:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Signature page has no ID",
        )
        return None

    # Get signature name
    props = signature_page.get("properties", {}) or {}
    name_prop = props.get("Name", {}).get("title", []) or []
    signature_name = name_prop[0].get("plain_text", "") if name_prop else "Unknown"

    # Get description if available
    desc_prop = props.get("Description", {}) or {}
    if desc_prop.get("rich_text"):
        description = desc_prop["rich_text"][0].get("plain_text", "")
    else:
        description = None

    # Fetch components from Signature Components DB
    if (
        not hasattr(cfg.notion, "signature_component_db_id")
        or not cfg.notion.signature_component_db_id
    ):
        logger.debug(
            "[INGEST][SIGNATURE-MATCH] Component DB ID not configured",
        )
        return None

    components: List[SignatureComponent] = []

    try:
        # Query for components linked to this signature
        url = f"{cfg.notion.base_url}/databases/{cfg.notion.signature_component_db_id}/query"
        all_components = []
        has_more = True
        start_cursor = None

        while has_more:
            payload = {
                "filter": {
                    "property": "Signature",
                    "relation": {"contains": signature_page_id},
                },
                "page_size": 100,
            }
            if start_cursor:
                payload["start_cursor"] = start_cursor

            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()

            data = resp.json()
            all_components.extend(data.get("results", []))
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")

        # Parse components
        for comp_page in all_components:
            comp_props = comp_page.get("properties", {}) or {}

            # Get species name (Component Name)
            name_prop = comp_props.get("Component Name", {}).get("title", []) or []
            species = name_prop[0].get("plain_text", "") if name_prop else ""

            if not species:
                continue

            # Get direction
            direction_prop = comp_props.get("Direction", {}).get("select")
            direction = None
            if direction_prop:
                direction_str = direction_prop.get("name", "")
                # Map Notion values to signature format
                direction_map = {
                    "Up": "↑",
                    "Down": "↓",
                    "NoChange": "neutral",
                    "Complex": "complex",
                    "Unknown": None,
                }
                direction = direction_map.get(direction_str, None)

            # Get weight
            weight_prop = comp_props.get("Weight", {}).get("number")
            weight = float(weight_prop) if weight_prop is not None else None

            # Get feature_type (for multi-omics support)
            feature_type_prop = comp_props.get("Feature Type", {}).get("select")
            feature_type = None
            if feature_type_prop:
                feature_type_str = feature_type_prop.get("name", "").lower()
                # Normalize to lowercase (gene, protein, metabolite, lipid)
                feature_type = feature_type_str
            else:
                # Default to lipid for backward compatibility
                feature_type = "lipid"

            components.append(
                SignatureComponent(
                    feature_name=species,  # Use feature_name for multi-omics
                    feature_type=feature_type,
                    direction=direction,
                    weight=weight,
                )
            )

        if not components:
            logger.debug(
                "[INGEST][SIGNATURE-MATCH] No components found for signature %s",
                signature_name,
            )
            return None

        signature = Signature(
            name=signature_name,
            components=components,
            description=description,
        )

        logger.debug(
            "[INGEST][SIGNATURE-MATCH] Loaded signature '%s' with %d components",
            signature_name,
            len(components),
        )

        return signature

    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Error loading signature %s: %r",
            signature_name,
            e,
        )
        return None

