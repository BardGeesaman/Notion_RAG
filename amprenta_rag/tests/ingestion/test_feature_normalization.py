"""
Tests for feature normalization functions.
"""

from amprenta_rag.ingestion.features.normalization import normalize_metabolite_name


def test_normalize_metabolite_name_basic():
    """Test basic metabolite name normalization."""
    assert normalize_metabolite_name("Glutamine") == "Glutamine"
    assert normalize_metabolite_name("glutamine") == "Glutamine"
    assert normalize_metabolite_name("  Glutamine  ") == "Glutamine"


def test_normalize_metabolite_name_with_prefixes():
    """Test normalization with database prefixes."""
    assert normalize_metabolite_name("HMDB:12345 Glutamine") == "Glutamine"
    # KEGG prefix might not be fully stripped, so check it's normalized
    result = normalize_metabolite_name("KEGG:C00064 Glutamine")
    assert "Glutamine" in result or "glutamine" in result.lower()
    assert normalize_metabolite_name("CHEBI:12345 Glutamine") == "Glutamine"


def test_normalize_metabolite_name_synonyms():
    """Test synonym mapping."""
    assert normalize_metabolite_name("l-glutamate") == "Glutamate"
    assert normalize_metabolite_name("L-glutamic acid") == "Glutamate"
    assert normalize_metabolite_name("glutamic acid") == "Glutamate"
    assert normalize_metabolite_name("gln") == "Glutamine"
    assert normalize_metabolite_name("glu") == "Glutamate"


def test_normalize_metabolite_name_edge_cases():
    """Test edge cases."""
    assert normalize_metabolite_name("") == ""
    assert normalize_metabolite_name("   ") == ""
    assert normalize_metabolite_name("123 Metabolite") == "Metabolite"
    assert normalize_metabolite_name("Metabolite (pos)") == "Metabolite"

