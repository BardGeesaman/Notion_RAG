"""
Signature page CRUD operations.

Notion support has been removed - Postgres is now the source of truth.
These functions are stubs that log and return without action.
"""

from __future__ import annotations

from typing import Optional

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature

logger = get_logger(__name__)


def find_or_create_signature_page(
    signature: Signature,
    signature_type: str = "Literature-derived",
    data_ownership: str = "Public",
    version: Optional[str] = None,
    description: Optional[str] = None,
) -> Optional[str]:
    """
    Stub: Notion support removed. Returns None.
    
    Previously found or created a Lipid Signature page in Notion.
    """
    logger.debug(
        "[SIGNATURES][SIGNATURE-CRUD] find_or_create_signature_page() is a no-op (Notion removed)"
    )
    return None


def update_signature_validation_metrics(
    signature_page_id: str,
    validation_results: dict,
) -> bool:
    """
    Stub: Notion support removed. Returns False.
    
    Previously updated validation metrics on a Signature page in Notion.
    """
    logger.debug(
        "[SIGNATURES][SIGNATURE-CRUD] update_signature_validation_metrics() is a no-op (Notion removed)"
    )
    return False


def update_signature_modalities(
    signature_page_id: str,
    modalities: list[str],
) -> bool:
    """
    Stub: Notion support removed. Returns False.
    
    Previously updated modalities on a Signature page in Notion.
    """
    logger.debug(
        "[SIGNATURES][SIGNATURE-CRUD] update_signature_modalities() is a no-op (Notion removed)"
    )
    return False


def update_signature_page_if_needed(
    signature_page_id: str,
    updates: dict,
) -> bool:
    """
    Stub: Notion support removed. Returns False.
    
    Previously updated a Signature page in Notion if changes were detected.
    """
    logger.debug(
        "[SIGNATURES][SIGNATURE-CRUD] update_signature_page_if_needed() is a no-op (Notion removed)"
    )
    return False
