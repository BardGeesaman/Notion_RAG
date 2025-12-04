"""
Metabolomics ingestion pipeline.

This package provides functions for ingesting internal metabolomics datasets,
including metabolite name normalization, file parsing, Notion integration, and RAG embedding.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.metabolomics.embedding import embed_metabolomics_dataset
from amprenta_rag.ingestion.metabolomics.file_parsing import extract_metabolite_set_from_file
from amprenta_rag.ingestion.metabolomics.ingestion import (
    create_metabolomics_dataset_page,
    ingest_metabolomics_file,
)
from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name

__all__ = [
    "ingest_metabolomics_file",
    "normalize_metabolite_name",
    "extract_metabolite_set_from_file",
    "create_metabolomics_dataset_page",
    "embed_metabolomics_dataset",
]

