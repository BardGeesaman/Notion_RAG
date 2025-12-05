"""
RAG builders for Postgres + Notion hybrid architecture.

This package provides utilities for building RAG chunks and metadata
using Postgres as the source of truth, with Notion providing narrative context.
"""

from amprenta_rag.rag.postgres_builder import (
    build_dataset_rag_metadata,
    build_dataset_rag_text,
    build_experiment_rag_metadata,
    build_feature_rag_metadata,
    build_program_rag_metadata,
    build_program_rag_text,
    build_signature_rag_metadata,
)

__all__ = [
    "build_dataset_rag_metadata",
    "build_dataset_rag_text",
    "build_experiment_rag_metadata",
    "build_feature_rag_metadata",
    "build_program_rag_metadata",
    "build_program_rag_text",
    "build_signature_rag_metadata",
]

