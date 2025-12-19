"""
Signature metadata collection and reverse linking.

DEPRECATED: Notion support has been removed. Postgres is now the source of truth.
These functions are stubs for backward compatibility.
"""

from __future__ import annotations

from typing import Any, Dict, List

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def enforce_signature_reverse_link(
    signature_page_id: str, source_page_id: str
) -> None:
    """
    DEPRECATED: Notion support has been removed.
    """
    logger.debug(
        "[METADATA] enforce_signature_reverse_link() deprecated - Notion support removed"
    )


def collect_signature_metadata(signature_ids: List[str]) -> Dict[str, Any]:
    """
    DEPRECATED: Notion support has been removed.

    Returns empty metadata structure for backward compatibility.
    """
    logger.debug(
        "[METADATA] collect_signature_metadata() deprecated - Notion support removed"
    )
    return {
        "sig_short_ids": [],
        "sig_roles": [],
        "sig_axes": [],
        "sig_ownership": [],
    }
