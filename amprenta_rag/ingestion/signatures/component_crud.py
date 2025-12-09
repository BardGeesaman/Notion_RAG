"""
Component page CRUD operations.

Notion support has been removed - Postgres is now the source of truth.
These functions are stubs that log and return without action.
"""

from __future__ import annotations

from typing import List, Optional

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
    Stub: Notion support removed. Returns None.
    
    Previously found or created a Signature Component page in Notion.
    """
    logger.debug(
        "[SIGNATURES][COMPONENT-CRUD] find_or_create_component_page() is a no-op (Notion removed)"
    )
    return None
