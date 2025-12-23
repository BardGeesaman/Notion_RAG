"""
Notion CRUD operations for lipid signatures, components, and species.

DEPRECATED: Notion support has been removed. Postgres is now the source of truth.
These functions are stubs that return None/False for backward compatibility.
"""

from __future__ import annotations

from typing import List, Optional

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ingestion.signatures.short_id import generate_signature_short_id
from amprenta_rag.signatures.signature_loader import Signature, SignatureComponent

logger = get_logger(__name__)


def find_or_create_component_page(
    component: SignatureComponent,
    signature_page_id: str,
    disease_context: Optional[List[str]] = None,
    matrix: Optional[List[str]] = None,
) -> Optional[str]:
    """DEPRECATED: Notion support removed. Returns None."""
    logger.debug("[SIGNATURE-CRUD] find_or_create_component_page() is a no-op (Notion removed)")
    return None


def find_or_create_lipid_species_page(
    lipid_name: str,
) -> Optional[str]:
    """DEPRECATED: Notion support removed. Returns None."""
    logger.debug("[SIGNATURE-CRUD] find_or_create_lipid_species_page() is a no-op (Notion removed)")
    return None


def find_or_create_signature_page(
    signature: Signature,
    signature_type: str = "Literature-derived",
    data_ownership: str = "Public",
    version: Optional[str] = None,
    description: Optional[str] = None,
) -> Optional[str]:
    """DEPRECATED: Notion support removed. Returns None."""
    logger.debug("[SIGNATURE-CRUD] find_or_create_signature_page() is a no-op (Notion removed)")
    return None


def update_lipid_species_synonyms(
    species_page_id: str,
    synonyms: List[str],
) -> bool:
    """DEPRECATED: Notion support removed. Returns False."""
    del species_page_id, synonyms
    logger.debug("[SIGNATURE-CRUD] update_lipid_species_synonyms() is a no-op (Notion removed)")
    return False


def update_signature_modalities(
    signature_page_id: str,
    modalities: List[str],
) -> bool:
    """DEPRECATED: Notion support removed. Returns False."""
    logger.debug("[SIGNATURE-CRUD] update_signature_modalities() is a no-op (Notion removed)")
    return False


def update_signature_page_if_needed(
    signature_page_id: str,
    updates: dict,
) -> bool:
    """DEPRECATED: Notion support removed. Returns False."""
    logger.debug("[SIGNATURE-CRUD] update_signature_page_if_needed() is a no-op (Notion removed)")
    return False


__all__ = [
    "generate_signature_short_id",
    "find_or_create_signature_page",
    "update_signature_page_if_needed",
    "update_signature_modalities",
    "find_or_create_component_page",
    "find_or_create_lipid_species_page",
    "update_lipid_species_synonyms",
]
