"""
Tests for signature short ID generation.
"""

from amprenta_rag.ingestion.signatures.short_id import generate_signature_short_id


def test_generate_signature_short_id_basic():
    """Test basic short ID generation."""
    result = generate_signature_short_id("ALS-CSF-Core-6Ceramides")
    assert result == "ALS-CSF-Core-6Ceramides"
    assert len(result) <= 50


def test_generate_signature_short_id_with_version():
    """Test short ID generation with version."""
    result = generate_signature_short_id("Test Signature", version="1.0")
    assert result == "Test-Signature-v1.0"


def test_generate_signature_short_id_special_chars():
    """Test short ID generation with special characters."""
    result = generate_signature_short_id("Test/Signature@2024")
    assert "/" not in result
    assert "@" not in result
    assert "-" in result or result.isalnum()


def test_generate_signature_short_id_long_name():
    """Test short ID generation with very long name."""
    long_name = "A" * 100
    result = generate_signature_short_id(long_name)
    assert len(result) <= 50

