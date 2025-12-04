"""
Tests for feature extraction functions.
"""

from amprenta_rag.ingestion.features.extraction import (
    extract_features_from_mwtab,
    extract_features_from_text,
)


def test_extract_features_from_mwtab():
    """Test extracting features from mwTab JSON."""
    mwtab_json = {
        "MS_METABOLITE_DATA": {
            "Data": [
                {"Metabolite": "HMDB:12345 Glutamine"},
                {"Metabolite": "Citrate"},
            ],
        },
    }
    result = extract_features_from_mwtab(mwtab_json)
    assert len(result) >= 2
    assert "Glutamine" in result
    assert "Citrate" in result


def test_extract_features_from_text():
    """Test extracting features from text content."""
    text = "The study found increased levels of glutamate and citrate in the samples."
    result = extract_features_from_text(text)
    # Should find at least glutamate (amino acid)
    assert len(result) > 0


def test_extract_features_from_text_amino_acids():
    """Test extracting amino acids from text."""
    text = "Alanine, arginine, and serine were measured."
    result = extract_features_from_text(text)
    assert "Alanine" in result or "alanine" in result.lower()
    assert "Serine" in result or "serine" in result.lower()

