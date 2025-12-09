"""
Cross-omics RAG reasoning capabilities.

This package provides high-level reasoning functions that generate multi-omics summaries
across Programs, Signatures, Features, and Datasets.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.query.cross_omics.dataset_summary_postgres import cross_omics_dataset_summary_postgres
from amprenta_rag.query.cross_omics.feature_summary_postgres import cross_omics_feature_summary_postgres
from amprenta_rag.query.cross_omics.program_summary_postgres import cross_omics_program_summary_postgres
from amprenta_rag.query.cross_omics.signature_summary_postgres import cross_omics_signature_summary_postgres

__all__ = [
    "cross_omics_program_summary_postgres",
    "cross_omics_signature_summary_postgres",
    "cross_omics_feature_summary_postgres",
    "cross_omics_dataset_summary_postgres",
]
