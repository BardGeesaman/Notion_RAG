"""
Feature extraction and linking for ingestion pipelines.

This package provides functions to:
- Normalize feature names (metabolites, etc.)
- Extract features from mwTab data and text
- Link features to Notion items (datasets, literature, experiments, emails)

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.features.constants import (
    AMINO_ACIDS,
    METABOLITE_SYNONYMS,
    NUCLEOTIDES,
)
from amprenta_rag.ingestion.features.extraction import (
    extract_features_from_mwtab,
    extract_features_from_text,
)
from amprenta_rag.ingestion.features.linking import (
    link_feature,
    link_features_to_notion_items,
)
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
    "link_feature",
    "link_features_to_notion_items",
]

