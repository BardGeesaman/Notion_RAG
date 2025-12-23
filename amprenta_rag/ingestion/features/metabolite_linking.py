"""
Metabolite-specific linking functions.

Notion support has been removed - Postgres is now the source of truth.
These functions are stubs that log and return without action.
"""

from __future__ import annotations

from typing import List, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _find_or_create_metabolite_page(metabolite_name: str) -> Optional[str]:
    """
    Stub: Notion support removed. Returns None.

    Previously found or created a Metabolite Features page in Notion.
    """
    logger.debug(
        "[FEATURES][METABOLITE-LINKING] _find_or_create_metabolite_page() is a no-op (Notion removed)"
    )
    return None


def link_metabolites_to_item(
    metabolite_names: List[str],
    item_page_id: str,
    item_type: str,
) -> int:
    """
    Stub: Notion support removed. Returns 0.

    Previously linked metabolites to a Notion item.
    """
    del item_page_id
    logger.debug(
        "[FEATURES][METABOLITE-LINKING] link_metabolites_to_item() is a no-op (Notion removed)"
    )
    return 0
