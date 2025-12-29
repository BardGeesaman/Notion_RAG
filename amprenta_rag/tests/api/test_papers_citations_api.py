"""
Unit tests for papers citations API endpoints.

Tests citation graph and enrichment endpoints with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestGetPaperCitations:
    """Tests for GET /api/v1/papers/{paper_id}/citations endpoint."""

    def test_get_citations_success(self):
        """Test successful citations retrieval."""
        paper_id = uuid4()
        citing_paper_id = uuid4()
        
        # Mock paper and citations
        mock_literature = MagicMock()
        mock_literature.id = paper_id
        
        mock_citing_paper = MagicMock()
        mock_citing_paper.id = citing_paper_id
        mock_citing_paper.doi = "10.1234/citing"
        mock_citing_paper.title = "Citing Paper"
        
        mock_citation = MagicMock()
        mock_citation.citing_paper = mock_citing_paper
        mock_citation.is_influential = True
        mock_citation.citation_context = "As shown in ref..."
        
        # Mock database session
        mock_session = MagicMock()
        
        # First query: verify paper exists
        mock_session.query.return_value.filter.return_value.first.return_value = mock_literature
        
        # Second query: get citations
        mock_citations_query = MagicMock()
        mock_citations_query.filter.return_value = mock_citations_query
        mock_citations_query.order_by.return_value = mock_citations_query
        mock_citations_query.offset.return_value = mock_citations_query
        mock_citations_query.limit.return_value = mock_citations_query
        mock_citations_query.all.return_value = [mock_citation]
        
        call_count = [0]
        def mock_query_side_effect(*args):
            call_count[0] += 1
            if call_count[0] == 1:
                # First call: Literature query
                return mock_session.query.return_value
            else:
                # Second call: PaperCitation query
                return mock_citations_query
        
        mock_session.query.side_effect = mock_query_side_effect
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get(f"/api/v1/papers/{paper_id}/citations")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["title"] == "Citing Paper"
            assert data[0]["is_influential"] is True
        finally:
            app.dependency_overrides.clear()

    def test_get_citations_paper_not_found(self):
        """Test citations for non-existent paper."""
        paper_id = uuid4()
        
        # Mock database session - paper not found
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get(f"/api/v1/papers/{paper_id}/citations")
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()


class TestGetPaperReferences:
    """Tests for GET /api/v1/papers/{paper_id}/references endpoint."""

    def test_get_references_success(self):
        """Test successful references retrieval."""
        paper_id = uuid4()
        cited_paper_id = uuid4()
        
        # Mock paper and references
        mock_literature = MagicMock()
        mock_literature.id = paper_id
        
        mock_cited_paper = MagicMock()
        mock_cited_paper.id = cited_paper_id
        mock_cited_paper.doi = "10.1234/cited"
        mock_cited_paper.title = "Cited Paper"
        
        mock_reference = MagicMock()
        mock_reference.cited_paper = mock_cited_paper
        mock_reference.cited_doi = None
        mock_reference.cited_title = None
        mock_reference.is_influential = False
        
        # Mock database session
        mock_session = MagicMock()
        
        # First query: verify paper exists
        mock_session.query.return_value.filter.return_value.first.return_value = mock_literature
        
        # Second query: get references
        mock_refs_query = MagicMock()
        mock_refs_query.filter.return_value = mock_refs_query
        mock_refs_query.order_by.return_value = mock_refs_query
        mock_refs_query.offset.return_value = mock_refs_query
        mock_refs_query.limit.return_value = mock_refs_query
        mock_refs_query.all.return_value = [mock_reference]
        
        call_count = [0]
        def mock_query_side_effect(*args):
            call_count[0] += 1
            if call_count[0] == 1:
                return mock_session.query.return_value
            else:
                return mock_refs_query
        
        mock_session.query.side_effect = mock_query_side_effect
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get(f"/api/v1/papers/{paper_id}/references")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["title"] == "Cited Paper"
        finally:
            app.dependency_overrides.clear()


class TestEnrichPaper:
    """Tests for POST /api/v1/papers/{paper_id}/enrich endpoint."""

    @patch("amprenta_rag.api.routers.papers.SemanticScholarRepository")
    @patch("amprenta_rag.api.routers.papers.OpenAlexRepository")
    def test_enrich_paper_success(self, mock_oa_class, mock_s2_class):
        """Test successful paper enrichment."""
        paper_id = uuid4()
        
        # Mock paper
        mock_literature = MagicMock()
        mock_literature.id = paper_id
        mock_literature.doi = "10.1234/test"
        mock_literature.pmid = "12345678"
        mock_literature.semantic_scholar_id = None
        mock_literature.openalex_id = None
        
        # Mock S2 repository
        mock_s2 = MagicMock()
        mock_s2_metadata = MagicMock()
        mock_s2_metadata.paper_id = "s2_abc123"
        mock_s2.fetch_metadata.return_value = mock_s2_metadata
        mock_s2.get_citations.return_value = []
        mock_s2.get_references.return_value = []
        mock_s2_class.return_value = mock_s2
        
        # Mock OA repository
        mock_oa = MagicMock()
        mock_oa_metadata = MagicMock()
        mock_oa_metadata.paper_id = "https://openalex.org/W123"
        mock_oa.get_work.return_value = mock_oa_metadata
        mock_oa_class.return_value = mock_oa
        
        # Mock database session
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_literature
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(
                f"/api/v1/papers/{paper_id}/enrich",
                json={"fetch_citations": False, "fetch_references": False},
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["enriched"] is True
            assert data["semantic_scholar_id"] == "s2_abc123"
            assert data["openalex_id"] == "https://openalex.org/W123"
        finally:
            app.dependency_overrides.clear()

    def test_enrich_paper_not_found(self):
        """Test enrichment for non-existent paper."""
        paper_id = uuid4()
        
        # Mock database session - paper not found
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(
                f"/api/v1/papers/{paper_id}/enrich",
                json={"fetch_citations": True, "fetch_references": True},
            )
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()

