"""
Unit tests for papers API endpoints.

Tests search, ingestion, and retrieval with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.ingestion.papers.base import PaperMetadata

client = TestClient(app)


class TestPapersSearch:
    """Tests for POST /api/v1/papers/search endpoint."""

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    def test_search_papers_success(self, mock_repo_class):
        """Test successful paper search."""
        # Mock repository
        mock_repo = MagicMock()
        mock_repo.search_papers.return_value = ["12345", "67890"]
        mock_repo.fetch_metadata.side_effect = [
            PaperMetadata(
                paper_id="12345",
                title="Test Paper 1",
                abstract="Abstract 1",
                authors=["Author A"],
                pmid="12345",
            ),
            PaperMetadata(
                paper_id="67890",
                title="Test Paper 2",
                abstract="Abstract 2",
                authors=["Author B"],
                pmid="67890",
            ),
        ]
        mock_repo_class.return_value = mock_repo

        # Make request
        response = client.post(
            "/api/v1/papers/search",
            json={"query": "sphingolipid", "source": "pubmed", "limit": 10, "offset": 0},
        )

        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 2
        assert len(data["results"]) == 2
        assert data["results"][0]["title"] == "Test Paper 1"
        assert data["results"][1]["title"] == "Test Paper 2"

    def test_search_papers_unsupported_source(self):
        """Test search with unsupported source."""
        response = client.post(
            "/api/v1/papers/search",
            json={"query": "test", "source": "biorxiv"},
        )

        assert response.status_code == 400
        assert "Unsupported source" in response.json()["detail"]


class TestPapersIngest:
    """Tests for POST /api/v1/papers/ingest endpoint."""

    def test_ingest_paper_missing_identifiers(self):
        """Test ingestion without PMID or DOI."""
        response = client.post(
            "/api/v1/papers/ingest",
            json={},
        )

        assert response.status_code == 400
        assert "pmid or doi" in response.json()["detail"].lower()

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    def test_ingest_paper_already_exists(self, mock_repo_class):
        """Test ingestion of existing paper."""
        # Mock database query to return existing literature
        mock_lit = MagicMock()
        mock_lit.id = uuid4()
        mock_lit.pmid = "12345"
        mock_lit.doi = "10.1234/test"
        mock_lit.title = "Existing Paper"

        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_lit

        def mock_get_db():
            yield mock_session

        # Use dependency override
        from amprenta_rag.api.dependencies import get_database_session

        app.dependency_overrides[get_database_session] = mock_get_db

        try:
            response = client.post(
                "/api/v1/papers/ingest",
                json={"pmid": "12345"},
            )

            assert response.status_code == 201
            data = response.json()
            assert data["already_exists"] is True
            assert data["pmid"] == "12345"
            assert data["title"] == "Existing Paper"
        finally:
            app.dependency_overrides.clear()

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    def test_ingest_paper_new(self, mock_repo_class):
        """Test ingestion of new paper."""
        # Mock repository
        mock_repo = MagicMock()
        mock_repo.fetch_metadata.return_value = PaperMetadata(
            paper_id="99999",
            title="New Paper",
            abstract="New abstract",
            authors=["New Author"],
            journal="New Journal",
            year=2025,
            doi="10.9999/new",
            pmid="99999",
            pmc_id="PMC99999",
            mesh_terms=["Term1", "Term2"],
            source="pubmed",
            url="https://pubmed.ncbi.nlm.nih.gov/99999/",
        )
        mock_repo_class.return_value = mock_repo

        # Mock database session
        mock_lit_id = uuid4()
        mock_session = MagicMock()
        # First call: check for existing (returns None)
        mock_session.query.return_value.filter.return_value.first.return_value = None

        # Mock add and flush
        def mock_flush():
            # Set ID after flush
            if mock_session.add.call_args:
                mock_session.add.call_args[0][0].id = mock_lit_id

        mock_session.flush = mock_flush

        def mock_get_db():
            yield mock_session

        # Use dependency override
        from amprenta_rag.api.dependencies import get_database_session

        app.dependency_overrides[get_database_session] = mock_get_db

        try:
            response = client.post(
                "/api/v1/papers/ingest",
                json={"pmid": "99999", "fetch_fulltext": False},
            )

            assert response.status_code == 201
            data = response.json()
            assert data["already_exists"] is False
            assert data["pmid"] == "99999"
            assert data["title"] == "New Paper"
        finally:
            app.dependency_overrides.clear()


class TestPapersGet:
    """Tests for GET /api/v1/papers/{paper_id} endpoint."""

    def test_get_paper_not_found(self):
        """Test retrieving non-existent paper."""
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None

        def mock_get_db():
            yield mock_session

        # Use dependency override
        from amprenta_rag.api.dependencies import get_database_session

        app.dependency_overrides[get_database_session] = mock_get_db

        try:
            paper_id = uuid4()
            response = client.get(f"/api/v1/papers/{paper_id}")

            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()

    def test_get_paper_success(self):
        """Test retrieving existing paper."""
        paper_id = uuid4()

        # Mock literature record
        mock_lit = MagicMock()
        mock_lit.id = paper_id
        mock_lit.pmid = "12345"
        mock_lit.pmc_id = "PMC12345"
        mock_lit.doi = "10.1234/test"
        mock_lit.title = "Test Paper"
        mock_lit.abstract = "Test abstract"
        mock_lit.journal = "Test Journal"
        mock_lit.year = 2024
        mock_lit.mesh_terms = ["MeSH1", "MeSH2"]

        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_lit

        def mock_get_db():
            yield mock_session

        # Use dependency override
        from amprenta_rag.api.dependencies import get_database_session

        app.dependency_overrides[get_database_session] = mock_get_db

        try:
            response = client.get(f"/api/v1/papers/{paper_id}")

            assert response.status_code == 200
            data = response.json()
            assert data["pmid"] == "12345"
            assert data["title"] == "Test Paper"
            assert data["journal"] == "Test Journal"
            assert len(data["mesh_terms"]) == 2
        finally:
            app.dependency_overrides.clear()

