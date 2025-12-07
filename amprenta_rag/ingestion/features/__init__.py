"""
Feature extraction and linking utilities.

This package provides functions for extracting features from various sources
and linking them to Notion items.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.features.constants import AMINO_ACIDS, METABOLITE_SYNONYMS, NUCLEOTIDES
from amprenta_rag.ingestion.features.extraction import extract_features_from_mwtab, extract_features_from_text

# Conditional imports for Notion-dependent linking modules (legacy)
try:
    from amprenta_rag.ingestion.features.general_linking import link_feature
except (ImportError, RuntimeError, ModuleNotFoundError):
    # Notion sync disabled - general_linking is not available
    link_feature = None

try:
    from amprenta_rag.ingestion.features.metabolite_linking import link_features_to_notion_items
except (ImportError, RuntimeError, ModuleNotFoundError):
    # Notion sync disabled - metabolite_linking is not available
    link_features_to_notion_items = None
from amprenta_rag.ingestion.features.normalization import normalize_metabolite_name

__all__ = [
    # Constants
    "METABOLITE_SYNONYMS",
    "AMINO_ACIDS",
    "NUCLEOTIDES",
    # Normalization
    "normalize_metabolite_name",
    # Extraction
    "extract_features_from_mwtab",
    "extract_features_from_text",
    # Linking
    "link_feature",  # May be None if Notion disabled
    "link_features_to_notion_items",
]
