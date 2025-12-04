"""
Compatibility wrapper for the refactored transcriptomics ingestion pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from transcriptomics_ingestion.py will continue to work
without changes.
"""

from __future__ import annotations

from amprenta_rag.ingestion.transcriptomics import (
    build_dge_text_representation,
    create_transcriptomics_dataset_page,
    embed_transcriptomics_dataset,
    extract_gene_set_from_file,
    ingest_transcriptomics_file,
    normalize_gene_identifier,
)

__all__ = [
    "ingest_transcriptomics_file",
    "normalize_gene_identifier",
    "extract_gene_set_from_file",
    "create_transcriptomics_dataset_page",
    "embed_transcriptomics_dataset",
    "build_dge_text_representation",
]
