"""
Proteomics ingestion pipeline.

This package provides functions for ingesting internal proteomics datasets,
including protein identifier normalization, file parsing, Notion integration, and RAG embedding.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.proteomics.embedding import embed_proteomics_dataset
from amprenta_rag.ingestion.proteomics.file_parsing import extract_protein_set_from_file
from amprenta_rag.ingestion.proteomics.ingestion import (
    create_proteomics_dataset_page,
    ingest_proteomics_file,
)
from amprenta_rag.ingestion.proteomics.normalization import normalize_protein_identifier

__all__ = [
    "ingest_proteomics_file",
    "normalize_protein_identifier",
    "extract_protein_set_from_file",
    "create_proteomics_dataset_page",
    "embed_proteomics_dataset",
]

