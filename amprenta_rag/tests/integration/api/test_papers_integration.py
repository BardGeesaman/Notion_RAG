"""Integration tests for papers API with real database."""

import pytest
from datetime import datetime
from unittest.mock import patch, MagicMock
from uuid import uuid4

from amprenta_rag.models.content import Literature, PaperCitation
from amprenta_rag.ingestion.papers.base import PaperMetadata


@pytest.mark.integration
class TestPapersAPIIntegration:
    """Integration tests for papers API endpoints."""

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    def test_search_papers_success(self, mock_repo_class, integration_client, 
                                 db_session, timed_request):
        """Test successful paper search with mocked external API."""
        # Mock PubMed repository (external service)
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
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/papers/search",
            "test_search_papers_success",
            json={"query": "sphingolipid", "source": "pubmed", "limit": 10, "offset": 0}
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["total"] == 2
        assert len(data["results"]) == 2
        assert data["results"][0]["title"] == "Test Paper 1"
        assert data["results"][1]["title"] == "Test Paper 2"

    def test_search_papers_unsupported_source(self, integration_client, timed_request):
        """Test search with unsupported source."""
        response, benchmark = timed_request(
            "POST",
            "/api/v1/papers/search",
            "test_search_unsupported_source",
            json={"query": "test", "source": "biorxiv"}
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "Unsupported source" in response.json()["detail"]

    def test_ingest_paper_missing_identifiers(self, integration_client, timed_request):
        """Test ingestion without PMID or DOI."""
        response, benchmark = timed_request(
            "POST",
            "/api/v1/papers/ingest",
            "test_ingest_missing_identifiers",
            json={}
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "pmid or doi" in response.json()["detail"].lower()

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    def test_ingest_paper_already_exists(self, mock_repo_class, integration_client, 
                                       db_session, timed_request):
        """Test ingestion of existing paper in real database."""
        # Create existing literature in database
        existing_lit = Literature(
            id=uuid4(),
            title="Existing Paper",
            pmid="12345",
            doi="10.1234/test",
            abstract="Existing abstract",
            source="pubmed",
            created_at=datetime.utcnow()
        )
        db_session.add(existing_lit)
        db_session.commit()
        db_session.refresh(existing_lit)
        
        # Mock repository (won't be called due to existing record)
        mock_repo = MagicMock()
        mock_repo_class.return_value = mock_repo
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/papers/ingest",
            "test_ingest_already_exists",
            json={"pmid": "12345"}
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["paper_id"] == str(existing_lit.id)
        assert data["status"] == "existing"
        assert data["title"] == "Existing Paper"

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    def test_ingest_paper_new(self, mock_repo_class, integration_client, 
                            db_session, timed_request):
        """Test ingestion of new paper with real database storage."""
        # Mock repository to return new paper metadata
        mock_repo = MagicMock()
        mock_metadata = PaperMetadata(
            paper_id="67890",
            title="New Research Paper",
            abstract="New abstract content",
            authors=["Dr. Smith", "Dr. Jones"],
            pmid="67890",
            doi="10.1234/new",
            journal="Nature",
            year=2024
        )
        mock_repo.fetch_metadata.return_value = mock_metadata
        mock_repo_class.return_value = mock_repo
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/papers/ingest",
            "test_ingest_new_paper",
            json={"pmid": "67890"}
        )
        
        assert response.status_code == 201
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["status"] == "created"
        assert data["title"] == "New Research Paper"
        
        # Verify paper was stored in database
        paper_id = data["paper_id"]
        stored_lit = db_session.query(Literature).filter_by(id=paper_id).first()
        assert stored_lit is not None
        assert stored_lit.title == "New Research Paper"
        assert stored_lit.pmid == "67890"
        assert stored_lit.doi == "10.1234/new"
        assert stored_lit.journal == "Nature"

    def test_get_paper_not_found(self, integration_client, timed_request):
        """Test retrieving non-existent paper."""
        fake_paper_id = uuid4()
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/papers/{fake_paper_id}",
            "test_get_paper_not_found"
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    def test_get_paper_success(self, integration_client, db_session, timed_request):
        """Test successful paper retrieval from real database."""
        # Create paper in database
        paper = Literature(
            id=uuid4(),
            title="Retrieved Paper",
            abstract="Paper abstract for retrieval",
            pmid="54321",
            doi="10.5678/retrieved",
            journal="Science",
            year=2023,
            citation_count=25,
            source="pubmed",
            created_at=datetime.utcnow()
        )
        db_session.add(paper)
        db_session.commit()
        db_session.refresh(paper)
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/papers/{paper.id}",
            "test_get_paper_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["id"] == str(paper.id)
        assert data["title"] == "Retrieved Paper"
        assert data["pmid"] == "54321"
        assert data["doi"] == "10.5678/retrieved"
        assert data["citation_count"] == 25

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    async def test_search_pubmed_async(self, mock_repo_class, integration_client, timed_request):
        """Test async PubMed search with mocked repository."""
        # Mock repository for async operation
        mock_repo = MagicMock()
        mock_repo.search_papers.return_value = ["11111"]
        mock_repo.fetch_metadata.return_value = PaperMetadata(
            paper_id="11111",
            title="Async Paper",
            abstract="Async abstract",
            authors=["Async Author"],
            pmid="11111"
        )
        mock_repo_class.return_value = mock_repo
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/papers/search",
            "test_search_pubmed_async",
            json={"query": "async test", "source": "pubmed", "limit": 5}
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["total"] == 1
        assert data["results"][0]["title"] == "Async Paper"

    @patch("amprenta_rag.api.routers.papers.OpenAlexRepository")
    async def test_search_openalex_async(self, mock_repo_class, integration_client, timed_request):
        """Test async OpenAlex search with mocked repository."""
        # Mock OpenAlex repository
        mock_repo = MagicMock()
        mock_repo.search_papers.return_value = ["W22222"]
        mock_repo.fetch_metadata.return_value = PaperMetadata(
            paper_id="W22222",
            title="OpenAlex Paper",
            abstract="OpenAlex abstract",
            authors=["OpenAlex Author"],
            openalex_id="W22222"
        )
        mock_repo_class.return_value = mock_repo
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/papers/search",
            "test_search_openalex_async",
            json={"query": "openalex test", "source": "openalex", "limit": 5}
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["total"] == 1
        assert data["results"][0]["title"] == "OpenAlex Paper"

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    async def test_ingest_async(self, mock_repo_class, integration_client, 
                              db_session, timed_request):
        """Test async paper ingestion with real database."""
        # Mock repository for async ingestion
        mock_repo = MagicMock()
        mock_metadata = PaperMetadata(
            paper_id="33333",
            title="Async Ingested Paper",
            abstract="Async ingested abstract",
            authors=["Async Author"],
            pmid="33333",
            doi="10.1111/async"
        )
        mock_repo.fetch_metadata.return_value = mock_metadata
        mock_repo_class.return_value = mock_repo
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/papers/ingest",
            "test_ingest_async",
            json={"pmid": "33333"}
        )
        
        assert response.status_code == 201
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["status"] == "created"
        assert data["title"] == "Async Ingested Paper"
        
        # Verify in database
        stored_paper = db_session.query(Literature).filter_by(pmid="33333").first()
        assert stored_paper is not None
        assert stored_paper.title == "Async Ingested Paper"

    @patch("amprenta_rag.api.routers.papers.enrich_paper_with_semantic_scholar")
    async def test_enrich_async(self, mock_enrich, integration_client, 
                              db_session, timed_request):
        """Test async paper enrichment with real database."""
        # Create paper in database
        paper = Literature(
            id=uuid4(),
            title="Paper to Enrich",
            pmid="44444",
            source="pubmed",
            created_at=datetime.utcnow()
        )
        db_session.add(paper)
        db_session.commit()
        db_session.refresh(paper)
        
        # Mock enrichment service (external API)
        mock_enrich.return_value = {
            "semantic_scholar_id": "S44444",
            "citation_count": 42,
            "tldr_summary": "This paper is about enrichment"
        }
        
        response, benchmark = timed_request(
            "POST",
            f"/api/v1/papers/{paper.id}/enrich",
            "test_enrich_async",
            json={"source": "semantic_scholar"}
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["status"] == "enriched"
        
        # Verify enrichment was stored in database
        db_session.refresh(paper)
        assert paper.semantic_scholar_id == "S44444"
        assert paper.citation_count == 42
        assert paper.tldr_summary == "This paper is about enrichment"

    @patch("amprenta_rag.api.routers.papers.PubMedRepository")
    async def test_concurrent_search(self, mock_repo_class, integration_client, 
                                   benchmark_tracker, timed_request):
        """Test concurrent paper searches with performance tracking."""
        # Mock repository for concurrent operations
        mock_repo = MagicMock()
        mock_repo.search_papers.return_value = ["55555"]
        mock_repo.fetch_metadata.return_value = PaperMetadata(
            paper_id="55555",
            title="Concurrent Paper",
            abstract="Concurrent abstract",
            authors=["Concurrent Author"],
            pmid="55555"
        )
        mock_repo_class.return_value = mock_repo
        
        # Perform concurrent searches
        search_queries = [
            {"query": f"concurrent test {i}", "source": "pubmed", "limit": 3}
            for i in range(3)
        ]
        
        responses = []
        for i, query in enumerate(search_queries):
            response, benchmark = timed_request(
                "POST",
                "/api/v1/papers/search",
                f"test_concurrent_search_{i}",
                json=query
            )
            responses.append((response, benchmark))
        
        # Verify all searches succeeded
        for response, benchmark in responses:
            assert response.status_code == 200
            assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
            
            data = response.json()
            assert data["total"] == 1
            assert data["results"][0]["title"] == "Concurrent Paper"

    def test_paper_citations_relationship(self, integration_client, db_session, timed_request):
        """Test paper citation relationships in real database."""
        # Create citing and cited papers
        cited_paper = Literature(
            id=uuid4(),
            title="Cited Paper",
            pmid="66666",
            source="pubmed",
            created_at=datetime.utcnow()
        )
        
        citing_paper = Literature(
            id=uuid4(),
            title="Citing Paper",
            pmid="77777",
            source="pubmed",
            created_at=datetime.utcnow()
        )
        
        db_session.add_all([cited_paper, citing_paper])
        db_session.commit()
        
        # Create citation relationship
        citation = PaperCitation(
            id=uuid4(),
            citing_paper_id=citing_paper.id,
            cited_paper_id=cited_paper.id,
            citation_context="This paper builds on previous work...",
            created_at=datetime.utcnow()
        )
        
        db_session.add(citation)
        db_session.commit()
        
        # Verify citation relationship exists in database
        stored_citation = db_session.query(PaperCitation).filter_by(
            citing_paper_id=citing_paper.id,
            cited_paper_id=cited_paper.id
        ).first()
        
        assert stored_citation is not None
        assert stored_citation.citation_context == "This paper builds on previous work..."
        
        # Verify papers can be retrieved with citation data
        citing_paper_response, benchmark = timed_request(
            "GET",
            f"/api/v1/papers/{citing_paper.id}",
            "test_citing_paper_retrieval"
        )
        
        assert citing_paper_response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
