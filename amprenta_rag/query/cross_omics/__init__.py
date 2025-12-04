"""
Cross-omics RAG reasoning capabilities.

This package provides high-level reasoning functions that generate multi-omics summaries
across Programs, Signatures, Features, and Datasets.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.query.cross_omics.dataset_summary import cross_omics_dataset_summary
from amprenta_rag.query.cross_omics.feature_summary import cross_omics_feature_summary
from amprenta_rag.query.cross_omics.program_summary import cross_omics_program_summary
from amprenta_rag.query.cross_omics.signature_summary import cross_omics_signature_summary

__all__ = [
    "cross_omics_program_summary",
    "cross_omics_signature_summary",
    "cross_omics_feature_summary",
    "cross_omics_dataset_summary",
]

