"""
Notion CRUD operations for signatures, components, and species.

This package handles all Notion database operations for:
- Lipid Signatures database
- Lipid Signature Components database
- Lipid Species database

All functions are idempotent - they find existing pages or create new ones.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.signatures.component_crud import find_or_create_component_page
from amprenta_rag.ingestion.signatures.short_id import generate_signature_short_id
from amprenta_rag.ingestion.signatures.signature_crud import (
    find_or_create_signature_page,
    update_signature_modalities,
    update_signature_page_if_needed,
)
from amprenta_rag.ingestion.signatures.species_crud import (
    find_or_create_lipid_species_page,
    update_lipid_species_synonyms,
)

__all__ = [
    "generate_signature_short_id",
    "find_or_create_signature_page",
    "update_signature_page_if_needed",
    "update_signature_modalities",
    "find_or_create_component_page",
    "find_or_create_lipid_species_page",
    "update_lipid_species_synonyms",
]

