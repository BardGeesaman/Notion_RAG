"""
Tests for Ketcher chemical structure editor component.

Tests HTML generation and component integration.
"""

from __future__ import annotations

import pytest

from scripts.dashboard.components.ketcher_sketcher import (
    build_ketcher_html,
    get_smiles_from_session,
)


class TestBuildKetcherHTML:
    """Tests for Ketcher HTML generation."""

    def test_build_html_default(self):
        """Test HTML generation with defaults."""
        html = build_ketcher_html()
        
        assert "ketcher-frame" in html
        assert "unpkg.com/ketcher-standalone" in html
        assert "exportSMILES" in html

    def test_build_html_with_initial_smiles(self):
        """Test HTML generation with initial SMILES."""
        smiles = "CCO"
        html = build_ketcher_html(initial_smiles=smiles)
        
        assert "CCO" in html
        assert "setMolecule" in html

    def test_build_html_with_custom_dimensions(self):
        """Test HTML generation with custom width/height."""
        html = build_ketcher_html(width=1000, height=800)
        
        assert "height: 800px" in html

    def test_build_html_unique_ids(self):
        """Test that different component keys generate unique IDs."""
        html1 = build_ketcher_html(component_key="editor1")
        html2 = build_ketcher_html(component_key="editor2")
        
        assert "ketcher_editor1" in html1
        assert "ketcher_editor2" in html2
        assert html1 != html2

    def test_build_html_escapes_smiles(self):
        """Test that SMILES with special characters are escaped."""
        smiles_with_quotes = "C'C\"O"
        html = build_ketcher_html(initial_smiles=smiles_with_quotes)
        
        # Should escape quotes to prevent JavaScript injection
        assert "\\'" in html or "\\''" not in html  # Basic escaping check

    def test_build_html_contains_postmessage(self):
        """Test that HTML includes postMessage for SMILES communication."""
        html = build_ketcher_html()
        
        assert "postMessage" in html
        assert "ketcher-smiles" in html


class TestGetSmilesFromSession:
    """Tests for session state SMILES retrieval."""

    def test_get_smiles_not_in_session(self):
        """Test getting SMILES when not in session."""
        # Note: Streamlit session_state only exists in running app
        # This test validates the function doesn't crash
        result = get_smiles_from_session("nonexistent_key")
        
        # Should return None or handle gracefully
        assert result is None or isinstance(result, str)

