"""
Unit and functional tests for OpenAlex client.

Includes mocked unit tests and @pytest.mark.slow functional tests with real API.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from amprenta_rag.ingestion.papers.openalex import OpenAlexRepository


class TestOpenAlexUnit:
    """Unit tests for OpenAlex client with mocked responses."""

    @patch("amprenta_rag.ingestion.papers.openalex.requests.get")
    def test_search_works_success(self, mock_get):
        """Test successful works search."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "results": [
                {"id": "https://openalex.org/W123", "title": "Work 1"},
                {"id": "https://openalex.org/W456", "title": "Work 2"},
            ]
        }
        mock_get.return_value = mock_response
        
        repo = OpenAlexRepository(email="test@example.com")
        work_ids = repo.search_works("test query", max_results=10)
        
        assert len(work_ids) == 2
        assert "https://openalex.org/W123" in work_ids
        assert "https://openalex.org/W456" in work_ids

    @patch("amprenta_rag.ingestion.papers.openalex.requests.get")
    def test_get_work_success(self, mock_get):
        """Test successful work metadata fetch."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "id": "https://openalex.org/W123",
            "title": "Test Work",
            "abstract_inverted_index": {"This": [0], "is": [1], "test": [2]},
            "authorships": [
                {"author": {"display_name": "Author One"}},
                {"author": {"display_name": "Author Two"}},
            ],
            "publication_year": 2023,
            "primary_location": {"source": {"display_name": "Test Journal"}},
            "doi": "https://doi.org/10.1234/test",
            "ids": {"pmid": "https://pubmed.ncbi.nlm.nih.gov/12345678/"},
        }
        mock_get.return_value = mock_response
        
        repo = OpenAlexRepository()
        metadata = repo.get_work("W123")
        
        assert metadata is not None
        assert metadata.title == "Test Work"
        assert len(metadata.authors) == 2
        assert metadata.doi == "10.1234/test"
        assert metadata.pmid == "12345678"

    @patch("amprenta_rag.ingestion.papers.openalex.requests.get")
    def test_get_work_not_found(self, mock_get):
        """Test work fetch for non-existent work."""
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response
        
        repo = OpenAlexRepository()
        metadata = repo.get_work("W999999")
        
        assert metadata is None

    @patch("amprenta_rag.ingestion.papers.openalex.requests.get")
    def test_get_author_success(self, mock_get):
        """Test successful author fetch."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "id": "https://openalex.org/A123",
            "display_name": "Test Author",
            "works_count": 50,
        }
        mock_get.return_value = mock_response
        
        repo = OpenAlexRepository()
        author = repo.get_author("A123")
        
        assert author is not None
        assert author["display_name"] == "Test Author"

    @patch("amprenta_rag.ingestion.papers.openalex.requests.get")
    def test_get_institution_success(self, mock_get):
        """Test successful institution fetch."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "id": "https://openalex.org/I123",
            "display_name": "Test University",
            "country_code": "US",
        }
        mock_get.return_value = mock_response
        
        repo = OpenAlexRepository()
        institution = repo.get_institution("I123")
        
        assert institution is not None
        assert institution["display_name"] == "Test University"

    def test_reconstruct_abstract(self):
        """Test abstract reconstruction from inverted index."""
        repo = OpenAlexRepository()
        
        inverted_index = {
            "This": [0],
            "is": [1],
            "a": [2],
            "test": [3],
            "abstract": [4],
        }
        
        abstract = repo._reconstruct_abstract(inverted_index)
        
        assert abstract == "This is a test abstract"


@pytest.mark.slow
class TestOpenAlexFunctional:
    """Functional tests with real OpenAlex API calls."""

    def test_search_works_real(self):
        """Test real works search."""
        repo = OpenAlexRepository(email="test@amprenta.ai")
        
        # Search for well-known topic
        work_ids = repo.search_works("deep learning", max_results=5)
        
        # Should get some results
        assert isinstance(work_ids, list)
        assert len(work_ids) > 0
        assert len(work_ids) <= 5
        # OpenAlex IDs start with W
        assert all("openalex.org/W" in wid for wid in work_ids)

    def test_get_work_real(self):
        """Test real work metadata fetch."""
        repo = OpenAlexRepository(email="test@amprenta.ai")
        
        # Use a known OpenAlex work ID (AlexNet paper)
        work_id = "W2100837269"
        
        metadata = repo.get_work(work_id)
        
        assert metadata is not None
        assert metadata.title is not None
        assert len(metadata.title) > 10
        assert metadata.source == "openalex"

    def test_get_work_not_found_real(self):
        """Test work fetch for non-existent work."""
        repo = OpenAlexRepository()
        
        # Use an invalid work ID
        metadata = repo.get_work("W99999999999")
        
        assert metadata is None

