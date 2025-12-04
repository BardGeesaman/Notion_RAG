"""
Cross-omics RAG reasoning capabilities.

High-level reasoning functions that generate multi-omics summaries across:
- Programs
- Signatures
- Features (genes, proteins, metabolites, lipids)
- Datasets

**Note**: This module has been refactored. All functions are now in the
`amprenta_rag.query.cross_omics` package. This file maintains backward
compatibility by re-exporting all public functions.
"""

from __future__ import annotations

# Re-export all public functions from the new modular structure
from amprenta_rag.query.cross_omics import (
    cross_omics_dataset_summary,
    cross_omics_feature_summary,
    cross_omics_program_summary,
    cross_omics_signature_summary,
)

__all__ = [
    "cross_omics_program_summary",
    "cross_omics_signature_summary",
    "cross_omics_feature_summary",
    "cross_omics_dataset_summary",
]
