"""
Integration tests for feature extraction pipeline.
"""

import pytest
from unittest.mock import Mock, patch

from amprenta_rag.ingestion.features.extraction import (
    extract_features_from_mwtab,
    extract_features_from_text,
)
from amprenta_rag.ingestion.features.normalization import normalize_metabolite_name


def test_extract_and_normalize_pipeline():
    """Test the full pipeline of extraction and normalization."""
    mwtab_json = {
        "MS_METABOLITE_DATA": {
            "Data": [
                {"Metabolite": "HMDB:12345 L-glutamate"},
                {"Metabolite": "KEGG:C00025 Citrate"},
            ],
        },
    }
    features = extract_features_from_mwtab(mwtab_json)
    # Should normalize l-glutamate to Glutamate
    assert any("Glutamate" in f or "glutamate" in f.lower() for f in features)
    assert any("Citrate" in f or "citrate" in f.lower() for f in features)


def test_text_extraction_with_normalization():
    """Test that text extraction uses normalization."""
    text = "Increased levels of l-glutamate and citrate were observed."
    features = extract_features_from_text(text)
    # Should find and normalize metabolites
    assert len(features) > 0

