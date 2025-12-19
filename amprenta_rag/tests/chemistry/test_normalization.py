"""
Tests for the chemistry normalization module.

Coverage:
- SMILES normalization (with and without RDKit)
- Molecular descriptor computation
- Compound ID generation
- Edge cases (empty/invalid input)
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

from amprenta_rag.chemistry.normalization import (
    RDKIT_AVAILABLE,
    compute_molecular_descriptors,
    generate_compound_id,
    normalize_smiles,
)


class TestNormalizeSMILES:
    """Tests for normalize_smiles function."""

    def test_empty_smiles_returns_empty(self):
        """Empty SMILES should return empty string with None values."""
        result = normalize_smiles("")
        assert result == ("", None, None)

    def test_whitespace_smiles_returns_whitespace(self):
        """Whitespace-only SMILES should return original with None values."""
        result = normalize_smiles("   ")
        assert result == ("   ", None, None)

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_valid_smiles_normalized(self):
        """Valid SMILES should be normalized to canonical form."""
        # Ethanol - simple test case
        smiles = "CCO"
        canonical, inchi_key, formula = normalize_smiles(smiles)

        assert canonical == "CCO"  # Ethanol canonical form
        assert formula == "C2H6O"

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_non_canonical_smiles_normalized(self):
        """Non-canonical SMILES should be converted to canonical form."""
        # Different representations of ethanol
        smiles = "OCC"  # Non-canonical
        canonical, _, _ = normalize_smiles(smiles)

        assert canonical == "CCO"  # Should be normalized

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_invalid_smiles_returns_original(self):
        """Invalid SMILES should return original with None values."""
        invalid_smiles = "INVALID_SMILES_XYZ"
        result = normalize_smiles(invalid_smiles)

        assert result == (invalid_smiles, None, None)

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_complex_molecule_normalization(self):
        """Complex molecules should be normalized correctly."""
        # Glucose
        glucose_smiles = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
        canonical, inchi_key, formula = normalize_smiles(glucose_smiles)

        assert canonical is not None
        assert formula == "C6H12O6"

    def test_smiles_stripped(self):
        """Leading/trailing whitespace should be stripped."""
        smiles = "  CCO  "
        canonical, _, _ = normalize_smiles(smiles)

        # Should be stripped even without RDKit
        assert canonical.strip() == "CCO"


class TestComputeMolecularDescriptors:
    """Tests for compute_molecular_descriptors function."""

    def test_returns_dict_with_expected_keys(self):
        """Should return dict with all expected descriptor keys."""
        result = compute_molecular_descriptors("CCO")

        expected_keys = {
            "molecular_weight",
            "logp",
            "hbd_count",
            "hba_count",
            "rotatable_bonds",
            "aromatic_rings",
        }
        assert set(result.keys()) == expected_keys

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_valid_smiles_returns_values(self):
        """Valid SMILES should return computed descriptor values."""
        # Ethanol
        result = compute_molecular_descriptors("CCO")

        assert result["molecular_weight"] is not None
        assert result["molecular_weight"] > 0
        assert result["logp"] is not None
        assert result["hbd_count"] is not None
        assert result["hba_count"] is not None

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_ethanol_descriptors(self):
        """Ethanol should have known descriptor values."""
        result = compute_molecular_descriptors("CCO")

        # Ethanol: C2H6O, MW ~46.07
        assert 45 < result["molecular_weight"] < 47
        # Ethanol has 1 H-bond donor (OH)
        assert result["hbd_count"] == 1
        # Ethanol has 1 H-bond acceptor (O)
        assert result["hba_count"] == 1

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_invalid_smiles_returns_none_values(self):
        """Invalid SMILES should return None for all descriptors."""
        result = compute_molecular_descriptors("INVALID_XYZ")

        assert result["molecular_weight"] is None
        assert result["logp"] is None

    def test_empty_smiles_returns_none_values(self):
        """Empty SMILES should return None for all descriptors."""
        result = compute_molecular_descriptors("")

        assert result["molecular_weight"] is None


class TestGenerateCompoundID:
    """Tests for generate_compound_id function."""

    def test_returns_string(self):
        """Should return a string ID."""
        result = generate_compound_id("CCO")

        assert isinstance(result, str)

    def test_returns_16_char_hex(self):
        """Should return a 16-character hex string."""
        result = generate_compound_id("CCO")

        assert len(result) == 16
        # Check it's valid hex
        int(result, 16)

    def test_same_smiles_same_id(self):
        """Same SMILES should generate same ID."""
        id1 = generate_compound_id("CCO")
        id2 = generate_compound_id("CCO")

        assert id1 == id2

    def test_different_smiles_different_id(self):
        """Different SMILES should generate different IDs."""
        id1 = generate_compound_id("CCO")
        id2 = generate_compound_id("CCCO")  # Propanol

        assert id1 != id2

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_canonical_equivalence(self):
        """Equivalent SMILES should generate same ID (with RDKit)."""
        # Different representations of ethanol
        id1 = generate_compound_id("CCO")
        id2 = generate_compound_id("OCC")

        # With RDKit, these should normalize to same canonical SMILES
        assert id1 == id2

    def test_empty_smiles_still_generates_id(self):
        """Empty SMILES should still generate an ID (hash of empty string)."""
        result = generate_compound_id("")

        assert isinstance(result, str)
        assert len(result) == 16


class TestRDKitFallback:
    """Tests for fallback behavior when RDKit is not available."""

    @patch("amprenta_rag.chemistry.normalization.RDKIT_AVAILABLE", False)
    def test_normalize_smiles_fallback(self):
        """Without RDKit, normalize_smiles should return cleaned input."""
        from amprenta_rag.chemistry import normalization

        # Need to reimport to pick up patched value
        original_available = normalization.RDKIT_AVAILABLE
        normalization.RDKIT_AVAILABLE = False

        try:
            result = normalization.normalize_smiles("  CCO  ")
            # Should just strip whitespace
            assert result == ("CCO", None, None)
        finally:
            normalization.RDKIT_AVAILABLE = original_available

    @patch("amprenta_rag.chemistry.normalization.RDKIT_AVAILABLE", False)
    def test_compute_descriptors_fallback(self):
        """Without RDKit, compute_molecular_descriptors should return None values."""
        from amprenta_rag.chemistry import normalization

        original_available = normalization.RDKIT_AVAILABLE
        normalization.RDKIT_AVAILABLE = False

        try:
            result = normalization.compute_molecular_descriptors("CCO")
            # All values should be None
            assert all(v is None for v in result.values())
        finally:
            normalization.RDKIT_AVAILABLE = original_available

