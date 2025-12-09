"""
Signature loading from Notion.

Notion support has been removed - Postgres is now the source of truth.
These functions are stubs that log and return without action.
"""

from __future__ import annotations

from typing import Any, Dict, List

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature

logger = get_logger(__name__)


def fetch_all_signatures_from_notion() -> List[Dict[str, Any]]:
    """
    Stub: Notion support removed. Returns empty list.
    
    Previously fetched all signature pages from Notion.
    """
    logger.debug(
        "[SIGNATURE-MATCHING][LOADER] fetch_all_signatures_from_notion() is a no-op (Notion removed)"
    )
    return []


def load_signatures_from_notion() -> List[Signature]:
    """
    Stub: Notion support removed. Returns empty list.
    
    Previously loaded signatures from Notion into Signature objects.
    """
    logger.debug(
        "[SIGNATURE-MATCHING][LOADER] load_signatures_from_notion() is a no-op (Notion removed)"
    )
    return []


def load_signature_from_notion_page(page: Dict[str, Any]) -> Signature | None:
    """
    Stub: Notion support removed. Returns None.
    
    Previously loaded a single signature from a Notion page.
    """
    logger.debug(
        "[SIGNATURE-MATCHING][LOADER] load_signature_from_notion_page() is a no-op (Notion removed)"
    )
    return None

