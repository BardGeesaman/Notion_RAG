"""
Lipidomics ingestion pipeline.

This package provides functions for ingesting internal lipidomics datasets,
including species normalization, file parsing, Notion integration, and RAG embedding.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.lipidomics.embedding import embed_lipidomics_dataset
from amprenta_rag.ingestion.lipidomics.file_parsing import extract_species_from_file
from amprenta_rag.ingestion.lipidomics.ingestion import (
    create_lipidomics_dataset_page,
    ingest_lipidomics_file,
)
from amprenta_rag.ingestion.lipidomics.normalization import normalize_lipid_species

__all__ = [
    "ingest_lipidomics_file",
    "normalize_lipid_species",
    "extract_species_from_file",
    "create_lipidomics_dataset_page",
    "embed_lipidomics_dataset",
]

