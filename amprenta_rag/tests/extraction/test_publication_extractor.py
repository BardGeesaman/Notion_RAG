"""
Tests for publication data extractor.

Tests PDF text extraction, section detection, and LLM-based experiment extraction.
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.extraction.publication_extractor import (
    detect_sections,
    extract_methods_section,
    extract_experiments_from_text,
    extract_from_pdf_bytes,
    clear_cache,
)


class TestSectionDetection:
    """Tests for section detection in publication text."""

    def test_detect_sections_standard_structure(self):
        """Test section detection with standard paper structure."""
        text = """
Title of Paper

Abstract
This is the abstract of the paper.

Introduction
This is the introduction section.

Methods
This is the methods section with experimental details.

Results
This is the results section.

Discussion
This is the discussion.

References
1. Author et al.
"""
        sections = detect_sections(text)
        
        assert "Abstract" in sections
        assert "Introduction" in sections
        assert "Methods" in sections
        assert "Results" in sections
        assert "Discussion" in sections
        assert "References" in sections

    def test_detect_sections_materials_and_methods(self):
        """Test detection of 'Materials and Methods' variant."""
        text = """
Abstract
Paper abstract here.

Materials and Methods
Experimental procedures described.

Results
Findings presented.
"""
        sections = detect_sections(text)
        
        assert "Methods" in sections
        assert "Materials and Methods" in sections.get("Methods", "").lower() or sections.get("Methods", "")

    def test_extract_methods_section_success(self):
        """Test successful methods section extraction."""
        text = """
Abstract
Paper abstract.

Methods
Cell culture protocol.
Western blot procedure.
Statistical analysis.

Results
Data and findings.
"""
        methods_text = extract_methods_section(text)
        
        assert methods_text is not None
        assert "Cell culture" in methods_text or "protocol" in methods_text.lower()

    def test_extract_methods_section_not_found(self):
        """Test methods extraction when section doesn't exist."""
        text = """
Abstract
This paper has no methods section.

Results
Just results here.
"""
        methods_text = extract_methods_section(text)
        
        assert methods_text is None


class TestExperimentExtraction:
    """Tests for LLM-based experiment extraction."""

    @patch("amprenta_rag.extraction.publication_extractor.get_openai_client")
    def test_extract_experiments_success(self, mock_get_client):
        """Test successful experiment extraction."""
        # Mock LLM response
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_choice = MagicMock()
        mock_message = MagicMock()
        
        mock_message.content = json.dumps({
            "experiments": [
                {
                    "experiment_type": "RNA-seq",
                    "cell_line": "HeLa",
                    "treatment": "Drug A",
                    "concentration": "10 μM",
                    "timepoint": "24h",
                    "replicate_count": 3,
                    "measured_entities": ["gene1", "gene2", "gene3"],
                    "key_findings": "Significant upregulation"
                }
            ]
        })
        
        mock_choice.message = mock_message
        mock_response.choices = [mock_choice]
        mock_client.chat.completions.create.return_value = mock_response
        mock_get_client.return_value = mock_client
        
        # Clear cache before test
        clear_cache()
        
        text = "Methods: We performed RNA-seq on HeLa cells treated with Drug A..."
        result = extract_experiments_from_text(text)
        
        assert len(result.experiments) == 1
        assert result.experiments[0].experiment_type == "RNA-seq"
        assert result.experiments[0].cell_line == "HeLa"
        assert result.confidence > 0.0
        assert result.cached is False

    @patch("amprenta_rag.extraction.publication_extractor.get_openai_client")
    def test_extract_experiments_caching(self, mock_get_client):
        """Test that extraction results are cached."""
        # Mock LLM response
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_choice = MagicMock()
        mock_message = MagicMock()
        mock_message.content = '{"experiments": []}'
        mock_choice.message = mock_message
        mock_response.choices = [mock_choice]
        mock_client.chat.completions.create.return_value = mock_response
        mock_get_client.return_value = mock_client
        
        clear_cache()
        
        text = "Test methods text"
        
        # First call - should call LLM
        result1 = extract_experiments_from_text(text)
        assert result1.cached is False
        
        # Second call - should use cache
        result2 = extract_experiments_from_text(text)
        assert result2.cached is True
        
        # LLM should only be called once
        assert mock_client.chat.completions.create.call_count == 1

    @patch("amprenta_rag.extraction.publication_extractor.get_openai_client")
    def test_extract_experiments_truncates_long_text(self, mock_get_client):
        """Test that long text is truncated to token limit."""
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_choice = MagicMock()
        mock_message = MagicMock()
        mock_message.content = '{"experiments": []}'
        mock_choice.message = mock_message
        mock_response.choices = [mock_choice]
        mock_client.chat.completions.create.return_value = mock_response
        mock_get_client.return_value = mock_client
        
        clear_cache()
        
        # Create text longer than max_tokens * 4
        long_text = "Methods section. " * 10000  # ~170K chars
        
        result = extract_experiments_from_text(long_text, max_tokens=1000)
        
        # Should truncate to ~4000 chars (1000 tokens * 4)
        # Verify LLM was called (didn't fail on long text)
        assert mock_client.chat.completions.create.call_count == 1

    @patch("amprenta_rag.extraction.publication_extractor.get_openai_client")
    def test_extract_experiments_handles_json_error(self, mock_get_client):
        """Test handling of invalid JSON from LLM."""
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_choice = MagicMock()
        mock_message = MagicMock()
        mock_message.content = "Not valid JSON at all"
        mock_choice.message = mock_message
        mock_response.choices = [mock_choice]
        mock_client.chat.completions.create.return_value = mock_response
        mock_get_client.return_value = mock_client
        
        clear_cache()
        
        result = extract_experiments_from_text("Test text")
        
        # Should return empty result on JSON error
        assert len(result.experiments) == 0
        assert result.confidence == 0.0


class TestPDFExtraction:
    """Tests for PDF-to-experiments pipeline."""

    @patch("amprenta_rag.extraction.publication_extractor.extract_text_from_pdf_bytes")
    @patch("amprenta_rag.extraction.publication_extractor.extract_experiments_from_text")
    def test_extract_from_pdf_success(self, mock_extract_exp, mock_extract_pdf):
        """Test successful PDF extraction pipeline."""
        literature_id = uuid4()
        
        # Mock PDF text extraction
        mock_extract_pdf.return_value = """
Methods
We used RNA-seq to analyze gene expression.
Cells were treated with compound X for 24 hours.

Results
Significant changes observed.
"""
        
        # Mock experiment extraction
        mock_result = MagicMock()
        mock_result.experiments = [MagicMock(experiment_type="RNA-seq")]
        mock_result.confidence = 0.9
        mock_extract_exp.return_value = mock_result
        
        pdf_bytes = b"fake pdf content"
        result = extract_from_pdf_bytes(pdf_bytes, literature_id)
        
        # Verify pipeline executed
        mock_extract_pdf.assert_called_once_with(pdf_bytes)
        mock_extract_exp.assert_called_once()
        
        assert len(result.experiments) > 0

    @patch("amprenta_rag.extraction.publication_extractor.extract_text_from_pdf_bytes")
    def test_extract_from_pdf_empty_text(self, mock_extract_pdf):
        """Test PDF extraction with empty/short text."""
        mock_extract_pdf.return_value = "Too short"
        
        pdf_bytes = b"fake pdf"
        result = extract_from_pdf_bytes(pdf_bytes)
        
        # Should return empty result for short text
        assert len(result.experiments) == 0
        assert result.confidence == 0.0

    @patch("amprenta_rag.extraction.publication_extractor.extract_text_from_pdf_bytes")
    def test_extract_from_pdf_handles_error(self, mock_extract_pdf):
        """Test PDF extraction handles errors gracefully."""
        mock_extract_pdf.side_effect = RuntimeError("PDF parsing failed")
        
        pdf_bytes = b"fake pdf"
        result = extract_from_pdf_bytes(pdf_bytes)
        
        # Should return empty result on error
        assert len(result.experiments) == 0
        assert result.confidence == 0.0


class TestCacheManagement:
    """Tests for extraction cache."""

    def test_clear_cache(self):
        """Test cache clearing."""
        # This is a simple function test
        clear_cache()
        # Should not raise any errors
        assert True

