"""
Transcriptomics ingestion pipeline.

This package provides functions for ingesting internal transcriptomics datasets,
including gene identifier normalization, file parsing, Notion integration, and RAG embedding.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.transcriptomics.embedding import embed_transcriptomics_dataset
from amprenta_rag.ingestion.transcriptomics.file_parsing import extract_gene_set_from_file
from amprenta_rag.ingestion.transcriptomics.ingestion import (
    create_transcriptomics_dataset_page,
    ingest_transcriptomics_file,
)
from amprenta_rag.ingestion.transcriptomics.normalization import normalize_gene_identifier
from amprenta_rag.ingestion.transcriptomics.text_building import build_dge_text_representation

__all__ = [
    "ingest_transcriptomics_file",
    "normalize_gene_identifier",
    "extract_gene_set_from_file",
    "create_transcriptomics_dataset_page",
    "embed_transcriptomics_dataset",
    "build_dge_text_representation",
]

