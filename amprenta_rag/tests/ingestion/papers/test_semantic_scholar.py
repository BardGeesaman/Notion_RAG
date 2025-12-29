"""
Unit and functional tests for Semantic Scholar client.

Includes mocked unit tests and @pytest.mark.slow functional tests with real API.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from amprenta_rag.ingestion.papers.semantic_scholar import SemanticScholarRepository


class TestSemanticScholarUnit:
    """Unit tests for Semantic Scholar client with mocked responses."""

    @patch("amprenta_rag.ingestion.papers.semantic_scholar.requests.get")
    def test_search_papers_success(self, mock_get):
        """Test successful paper search."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "data": [
                {"paperId": "abc123", "title": "Test Paper 1"},
                {"paperId": "def456", "title": "Test Paper 2"},
            ]
        }
        mock_get.return_value = mock_response
        
        repo = SemanticScholarRepository()
        paper_ids = repo.search_papers("test query", max_results=10)
        
        assert len(paper_ids) == 2
        assert "abc123" in paper_ids
        assert "def456" in paper_ids

    @patch("amprenta_rag.ingestion.papers.semantic_scholar.requests.get")
    def test_fetch_metadata_success(self, mock_get):
        """Test successful metadata fetch."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "paperId": "abc123",
            "title": "Test Paper",
            "abstract": "This is a test abstract.",
            "authors": [{"name": "Author One"}, {"name": "Author Two"}],
            "year": 2023,
            "venue": "Test Journal",
            "externalIds": {"DOI": "10.1234/test", "PubMed": "12345678"},
            "url": "https://www.semanticscholar.org/paper/abc123",
        }
        mock_get.return_value = mock_response
        
        repo = SemanticScholarRepository()
        metadata = repo.fetch_metadata("abc123")
        
        assert metadata is not None
        assert metadata.title == "Test Paper"
        assert len(metadata.authors) == 2
        assert metadata.doi == "10.1234/test"
        assert metadata.pmid == "12345678"

    @patch("amprenta_rag.ingestion.papers.semantic_scholar.requests.get")
    def test_fetch_metadata_not_found(self, mock_get):
        """Test metadata fetch for non-existent paper."""
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response
        
        repo = SemanticScholarRepository()
        metadata = repo.fetch_metadata("nonexistent")
        
        assert metadata is None

    @patch("amprenta_rag.ingestion.papers.semantic_scholar.requests.get")
    def test_get_citations_success(self, mock_get):
        """Test successful citations retrieval."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "data": [
                {"citingPaper": {"paperId": "cite1", "title": "Citing Paper 1"}},
                {"citingPaper": {"paperId": "cite2", "title": "Citing Paper 2"}},
            ]
        }
        mock_get.return_value = mock_response
        
        repo = SemanticScholarRepository()
        citations = repo.get_citations("abc123", limit=10)
        
        assert len(citations) == 2
        assert citations[0]["paperId"] == "cite1"

    @patch("amprenta_rag.ingestion.papers.semantic_scholar.requests.get")
    def test_get_references_success(self, mock_get):
        """Test successful references retrieval."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "data": [
                {"citedPaper": {"paperId": "ref1", "title": "Reference 1"}},
                {"citedPaper": {"paperId": "ref2", "title": "Reference 2"}},
            ]
        }
        mock_get.return_value = mock_response
        
        repo = SemanticScholarRepository()
        references = repo.get_references("abc123", limit=10)
        
        assert len(references) == 2
        assert references[0]["paperId"] == "ref1"

    @patch("amprenta_rag.ingestion.papers.semantic_scholar.requests.get")
    def test_get_recommendations_success(self, mock_get):
        """Test successful recommendations retrieval."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "recommendedPapers": [
                {"paperId": "rec1", "title": "Recommended 1"},
                {"paperId": "rec2", "title": "Recommended 2"},
            ]
        }
        mock_get.return_value = mock_response
        
        repo = SemanticScholarRepository()
        recommendations = repo.get_recommendations("abc123", limit=10)
        
        assert len(recommendations) == 2
        assert "rec1" in recommendations
        assert "rec2" in recommendations


@pytest.mark.slow
class TestSemanticScholarFunctional:
    """Functional tests with real Semantic Scholar API calls."""

    def test_search_papers_real(self):
        """Test real paper search."""
        repo = SemanticScholarRepository()
        
        # Search for well-known topic
        paper_ids = repo.search_papers("machine learning", max_results=5)
        
        # Should get some results
        assert isinstance(paper_ids, list)
        assert len(paper_ids) > 0
        assert len(paper_ids) <= 5

    def test_fetch_metadata_real(self):
        """Test real metadata fetch for a known paper."""
        repo = SemanticScholarRepository()
        
        # Use a well-known paper ID (AlexNet paper)
        paper_id = "abd1c342495432171beb7ca8fd9551ef13cbd0ff"
        
        metadata = repo.fetch_metadata(paper_id)
        
        assert metadata is not None
        assert metadata.title is not None
        assert len(metadata.title) > 10
        assert metadata.source == "semantic_scholar"

    def test_get_citations_real(self):
        """Test real citations fetch."""
        repo = SemanticScholarRepository()
        
        # Use a well-known paper ID
        paper_id = "abd1c342495432171beb7ca8fd9551ef13cbd0ff"
        
        citations = repo.get_citations(paper_id, limit=10)
        
        # This highly cited paper should have citations
        assert isinstance(citations, list)
        assert len(citations) > 0

