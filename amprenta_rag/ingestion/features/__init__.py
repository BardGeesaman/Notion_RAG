"""
Feature extraction and linking utilities.

This package provides functions for extracting features from various sources
and linking them to Notion items.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from typing import Any, Callable, Optional, cast

from amprenta_rag.ingestion.features.constants import AMINO_ACIDS, METABOLITE_SYNONYMS, NUCLEOTIDES
from amprenta_rag.ingestion.features.extraction import extract_features_from_mwtab, extract_features_from_text

link_feature: Any = None
link_features_to_notion_items: Any = None
try:
    from amprenta_rag.ingestion.features import general_linking as _gl  # type: ignore[import]
    link_feature = getattr(_gl, "link_feature", None)
except Exception:
    link_feature = None
try:
    from amprenta_rag.ingestion.features import metabolite_linking as _ml  # type: ignore[import]
    link_features_to_notion_items = getattr(_ml, "link_features_to_notion_items", None)
except Exception:
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
