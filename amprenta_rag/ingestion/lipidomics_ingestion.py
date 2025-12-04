"""
Compatibility wrapper for the refactored lipidomics ingestion pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from lipidomics_ingestion.py will continue to work
without changes.
"""

from __future__ import annotations

from amprenta_rag.ingestion.lipidomics import (
    create_lipidomics_dataset_page,
    embed_lipidomics_dataset,
    extract_species_from_file,
    ingest_lipidomics_file,
    normalize_lipid_species,
)

__all__ = [
    "ingest_lipidomics_file",
    "normalize_lipid_species",
    "extract_species_from_file",
    "create_lipidomics_dataset_page",
    "embed_lipidomics_dataset",
]
