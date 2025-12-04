"""
Compatibility wrapper for the refactored proteomics ingestion pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from proteomics_ingestion.py will continue to work
without changes.
"""

from __future__ import annotations

from amprenta_rag.ingestion.proteomics import (
    create_proteomics_dataset_page,
    embed_proteomics_dataset,
    extract_protein_set_from_file,
    ingest_proteomics_file,
    normalize_protein_identifier,
)

__all__ = [
    "ingest_proteomics_file",
    "normalize_protein_identifier",
    "extract_protein_set_from_file",
    "create_proteomics_dataset_page",
    "embed_proteomics_dataset",
]
