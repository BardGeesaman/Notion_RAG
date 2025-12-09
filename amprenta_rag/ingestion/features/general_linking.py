"""
General feature linking functions.

Notion support has been removed - Postgres is now the source of truth.
These functions are stubs that log and return without action.
"""

from __future__ import annotations

from typing import Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _find_or_create_feature_page(feature_type: str, feature_name: str) -> Optional[str]:
    """
    Stub: Notion support removed. Returns None.
    
    Previously found or created a feature page in Notion.
    """
    logger.debug(
        "[FEATURES][GENERAL-LINKING] _find_or_create_feature_page() is a no-op (Notion removed)"
    )
    return None


def link_feature_to_dataset(
    feature_page_id: str,
    dataset_page_id: str,
    feature_type: str,
) -> bool:
    """
    Stub: Notion support removed. Returns False.
    
    Previously linked a feature page to a dataset page in Notion.
    """
    logger.debug(
        "[FEATURES][GENERAL-LINKING] link_feature_to_dataset() is a no-op (Notion removed)"
    )
    return False
