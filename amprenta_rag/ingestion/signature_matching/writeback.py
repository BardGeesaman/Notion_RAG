"""
Notion writeback for signature matches.

Notion support has been removed - Postgres is now the source of truth.
These functions are stubs that log and return without action.
"""

from __future__ import annotations

from typing import List

from amprenta_rag.ingestion.signature_matching.models import SignatureMatchResult
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def update_dataset_with_signature_matches(
    dataset_page_id: str,
    matches: List[SignatureMatchResult],
) -> None:
    """
    Stub: Notion support removed. No-op.
    
    Previously updated a dataset's Notion page with signature match information.
    """
    logger.debug(
        "[SIGNATURE-MATCHING][WRITEBACK] update_dataset_with_signature_matches() is a no-op (Notion removed)"
    )
