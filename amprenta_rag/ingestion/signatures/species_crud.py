"""
Lipid species page CRUD operations.

Notion support has been removed - Postgres is now the source of truth.
These functions are stubs that log and return without action.
"""

from __future__ import annotations

from typing import Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def find_or_create_lipid_species_page(
    lipid_name: str,
) -> Optional[str]:
    """
    Stub: Notion support removed. Returns None.
    
    Previously found or created a Lipid Species page in Notion.
    """
    logger.debug(
        "[SIGNATURES][SPECIES-CRUD] find_or_create_lipid_species_page() is a no-op (Notion removed)"
    )
    return None


def update_lipid_species_synonyms(
    species_page_id: str,
    synonyms: list[str],
) -> bool:
    """
    Stub: Notion support removed. Returns False.
    
    Previously updated synonyms on a Lipid Species page in Notion.
    """
    logger.debug(
        "[SIGNATURES][SPECIES-CRUD] update_lipid_species_synonyms() is a no-op (Notion removed)"
    )
    return False
