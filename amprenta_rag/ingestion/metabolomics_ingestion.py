"""
Compatibility wrapper for the refactored metabolomics ingestion pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from metabolomics_ingestion.py will continue to work
without changes.
"""

from __future__ import annotations

from amprenta_rag.ingestion.metabolomics import (
    create_metabolomics_dataset_page,
    embed_metabolomics_dataset,
    extract_metabolite_set_from_file,
    ingest_metabolomics_file,
    normalize_metabolite_name,
)

__all__ = [
    "ingest_metabolomics_file",
    "normalize_metabolite_name",
    "extract_metabolite_set_from_file",
    "create_metabolomics_dataset_page",
    "embed_metabolomics_dataset",
]
