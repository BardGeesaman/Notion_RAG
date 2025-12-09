"""
Signature linking utilities.

This module previously handled linking operations between signatures, components,
species, and Notion entities. Notion support has been removed - Postgres is now
the source of truth. These functions are stubs that log and return without action.
"""

from __future__ import annotations

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

__all__ = [
    "link_component_to_lipid_species",
    "link_signature_to_source",
    "link_component_to_metabolite_feature",
    "link_component_to_feature",
]


def link_component_to_lipid_species(
    component_page_id: str,
    lipid_species_page_id: str,
) -> None:
    """
    Stub: Notion linking removed. No-op.
    
    Previously linked component to lipid species in Notion.
    """
    logger.debug(
        "[SIGNATURE-LINKING] link_component_to_lipid_species() is a no-op (Notion removed)"
    )


def link_signature_to_source(
    signature_page_id: str,
    source_page_id: str,
    source_type: str,
) -> None:
    """
    Stub: Notion linking removed. No-op.
    
    Previously created reverse link from signature page to source page in Notion.
    """
    logger.debug(
        "[SIGNATURE-LINKING] link_signature_to_source() is a no-op (Notion removed)"
    )


def link_component_to_metabolite_feature(
    component_page_id: str,
    lipid_species_page_id: str,
) -> None:
    """
    Stub: Notion linking removed. No-op.
    
    Previously linked component's lipid species to Metabolite Feature page in Notion.
    """
    logger.debug(
        "[SIGNATURE-LINKING] link_component_to_metabolite_feature() is a no-op (Notion removed)"
    )


def link_component_to_feature(
    component_page_id: str,
    feature_type: str,
    feature_name: str,
) -> None:
    """
    Stub: Notion linking removed. No-op.
    
    Previously linked signature component to feature page in Notion.
    """
    logger.debug(
        "[SIGNATURE-LINKING] link_component_to_feature() is a no-op (Notion removed)"
    )
