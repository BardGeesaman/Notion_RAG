"""
Tests for lipid species mapping functions.
"""

from amprenta_rag.ingestion.signature_matching.species_mapping import (
    map_raw_lipid_to_canonical_species,
)


def test_map_raw_lipid_cer():
    """Test mapping CER vendor format to canonical."""
    result = map_raw_lipid_to_canonical_species("CER 16:0")
    assert result == "Cer(d18:1/16:0)"

    result = map_raw_lipid_to_canonical_species("cer 24:1")
    assert result == "Cer(d18:1/24:1)"


def test_map_raw_lipid_sm():
    """Test mapping SM vendor format to canonical."""
    result = map_raw_lipid_to_canonical_species("SM 16:0")
    assert result == "SM(d18:1/16:0)"

    result = map_raw_lipid_to_canonical_species("sm 24:1")
    assert result == "SM(d18:1/24:1)"


def test_map_raw_lipid_already_canonical():
    """Test that already canonical names are returned as-is."""
    result = map_raw_lipid_to_canonical_species("Cer(d18:1/16:0)")
    assert result == "Cer(d18:1/16:0)"


def test_map_raw_lipid_invalid():
    """Test that invalid names return None."""
    result = map_raw_lipid_to_canonical_species("Invalid Lipid")
    assert result is None

    result = map_raw_lipid_to_canonical_species("")
    assert result is None

    result = map_raw_lipid_to_canonical_species(None)
    assert result is None

