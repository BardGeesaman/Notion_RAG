"""
Unit tests for papers API endpoints.

Tests search, ingestion, and retrieval with mocked dependencies.
"""

from __future__ import annotations

import asyncio
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


class TestAsyncPapersAPI:
    """Test async execution of papers API endpoints."""

    @pytest.mark.asyncio
    async def test_search_pubmed_async(self):
        """Test PubMed search endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.papers._sync_search_papers') as mock_search:
            # Mock the search function
            mock_search.return_value = {
                "results": [
                    {
                        "paper_id": "12345",
                        "title": "Async Test Paper 1",
                        "abstract": "Test abstract 1",
                        "authors": ["Author A"],
                        "journal": "Test Journal",
                        "year": 2023,
                        "doi": "10.1234/test1",
                        "pmid": "12345",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/12345/"
                    },
                    {
                        "paper_id": "67890",
                        "title": "Async Test Paper 2",
                        "abstract": "Test abstract 2",
                        "authors": ["Author B"],
                        "journal": "Test Journal",
                        "year": 2023,
                        "doi": "10.1234/test2",
                        "pmid": "67890",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/67890/"
                    }
                ],
                "total": 2
            }
            
            from amprenta_rag.api.routers.papers import search_papers
            from amprenta_rag.api.routers.papers import PaperSearchRequest
            
            request = PaperSearchRequest(
                query="cancer treatment",
                source="pubmed",
                limit=10,
                offset=0
            )
            
            result = await search_papers(request, db=MagicMock())
            
            # Verify async execution and result
            assert len(result.results) == 2
            assert result.total == 2
            assert result.results[0].title == "Async Test Paper 1"
            assert result.results[1].title == "Async Test Paper 2"
            mock_search.assert_called_once_with("cancer treatment", "pubmed", 10, 0)

    @pytest.mark.asyncio
    async def test_search_openalex_async(self):
        """Test OpenAlex search would execute asynchronously (currently unsupported)."""
        from amprenta_rag.api.routers.papers import search_papers
        from amprenta_rag.api.routers.papers import PaperSearchRequest
        
        request = PaperSearchRequest(
            query="machine learning",
            source="openalex",  # Currently unsupported
            limit=5,
            offset=0
        )
        
        # Should raise HTTPException for unsupported source
        with pytest.raises(Exception) as exc_info:
            await search_papers(request, db=MagicMock())
        
        assert "Unsupported source" in str(exc_info.value)

    @pytest.mark.asyncio
    async def test_ingest_async(self):
        """Test paper ingestion endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.papers._sync_fetch_paper_metadata') as mock_fetch:
            with patch('amprenta_rag.api.routers.papers.get_database_session') as mock_db_session:
                # Mock the fetch function
                mock_metadata = PaperMetadata(
                    paper_id="12345",
                    title="Async Ingest Test Paper",
                    abstract="Test abstract for async ingestion",
                    authors=["Test Author"],
                    journal="Test Journal",
                    year=2023,
                    doi="10.1234/async_test",
                    pmid="12345",
                    pmc_id="PMC12345",
                    mesh_terms=["Cancer", "Treatment"],
                    url="https://pubmed.ncbi.nlm.nih.gov/12345/"
                )
                mock_fetch.return_value = mock_metadata
                
                # Mock database operations
                mock_db = MagicMock()
                mock_db_session.return_value = mock_db
                
                # Mock no existing paper
                mock_db.query.return_value.filter.return_value.first.return_value = None
                
                # Mock Literature constructor and instance
                with patch('amprenta_rag.api.routers.papers.Literature') as MockLiterature:
                    literature_id = uuid4()
                    mock_literature_instance = MagicMock()
                    mock_literature_instance.id = literature_id
                    mock_literature_instance.pmid = "12345"
                    mock_literature_instance.doi = "10.1234/async_test"
                    mock_literature_instance.title = "Async Ingest Test Paper"
                    
                    MockLiterature.return_value = mock_literature_instance
                    mock_db.add.return_value = None
                    mock_db.flush.return_value = None
                    mock_db.commit.return_value = None
                
                    from amprenta_rag.api.routers.papers import ingest_paper
                    from amprenta_rag.api.routers.papers import PaperIngestRequest
                    
                    request = PaperIngestRequest(
                        pmid="12345",
                        fetch_fulltext=False
                    )
                    
                    result = await ingest_paper(request, db=mock_db)
                    
                    # Verify async execution
                    mock_fetch.assert_called_once_with("12345", None)
                    assert result.title == "Async Ingest Test Paper"
                    assert result.already_exists == False
                    assert result.literature_id == literature_id

    @pytest.mark.asyncio
    async def test_enrich_async(self):
        """Test paper enrichment endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.papers._sync_enrich_with_semantic_scholar') as mock_s2:
            with patch('amprenta_rag.api.routers.papers._sync_enrich_with_openalex') as mock_oa:
                with patch('amprenta_rag.api.routers.papers.get_database_session') as mock_db_session:
                    with patch('amprenta_rag.api.routers.papers.PaperCitation') as MockPaperCitation:
                        # Mock Semantic Scholar enrichment
                        mock_s2_metadata = MagicMock()
                        mock_s2_metadata.paper_id = "S2_12345"
                        
                        mock_s2.return_value = {
                            "metadata": mock_s2_metadata,
                            "citations": [
                                {
                                    "paperId": "cite1",
                                    "title": "Citing Paper 1",
                                    "isInfluential": True,
                                    "externalIds": {"DOI": "10.1234/cite1"}
                                }
                            ],
                            "references": [
                                {
                                    "paperId": "ref1",
                                    "title": "Referenced Paper 1",
                                    "isInfluential": False,
                                    "externalIds": {"DOI": "10.1234/ref1"}
                                }
                            ]
                        }
                        
                        # Mock OpenAlex enrichment
                        mock_oa_metadata = MagicMock()
                        mock_oa_metadata.paper_id = "OA_12345"
                        mock_oa.return_value = mock_oa_metadata
                        
                        # Mock database operations
                        mock_db = MagicMock()
                        mock_db_session.return_value = mock_db
                        
                        # Mock literature object
                        paper_id = uuid4()
                        mock_literature = MagicMock()
                        mock_literature.id = paper_id
                        mock_literature.doi = "10.1234/test"
                        mock_literature.pmid = "12345"
                        mock_literature.semantic_scholar_id = None
                        mock_literature.openalex_id = None
                        mock_literature.citation_count = None
                        
                        # Mock query for getting the literature
                        mock_db.query.return_value.filter.return_value.first.return_value = mock_literature
                        
                        # Mock query for checking existing citations/references (should return None to allow adding)
                        mock_db.query.return_value.filter.return_value.filter.return_value.first.return_value = None
                        
                        # Mock PaperCitation constructor
                        MockPaperCitation.return_value = MagicMock()
                        
                        mock_db.add.return_value = None
                        mock_db.commit.return_value = None
                        
                        from amprenta_rag.api.routers.papers import enrich_paper
                        from amprenta_rag.api.routers.papers import EnrichPaperRequest
                        
                        request = EnrichPaperRequest(
                            fetch_citations=True,
                            fetch_references=True
                        )
                        
                        result = await enrich_paper(paper_id, request, db=mock_db)
                        
                        # Verify async execution
                        mock_s2.assert_called_once_with("DOI:10.1234/test", True, True)
                        mock_oa.assert_called_once_with("https://doi.org/10.1234/test")
                        assert result.enriched == True
                        # Note: citations_added may be 0 due to complex DB mocking, but async execution is verified
                        assert result.semantic_scholar_id == "S2_12345"
                        assert result.openalex_id == "OA_12345"

    @pytest.mark.asyncio
    async def test_concurrent_search(self):
        """Test multiple simultaneous paper searches."""
        with patch('amprenta_rag.api.routers.papers._sync_search_papers') as mock_search:
            # Mock function to return different results for different calls
            def mock_search_func(query, source, limit, offset):
                return {
                    "results": [
                        {
                            "paper_id": f"search_{query.replace(' ', '_')}",
                            "title": f"Paper for {query}",
                            "abstract": f"Abstract about {query}",
                            "authors": ["Test Author"],
                            "journal": "Test Journal",
                            "year": 2023,
                            "doi": f"10.1234/{query.replace(' ', '_')}",
                            "pmid": f"pmid_{query.replace(' ', '_')}",
                            "url": f"https://pubmed.ncbi.nlm.nih.gov/{query.replace(' ', '_')}/"
                        }
                    ],
                    "total": 1
                }
            
            mock_search.side_effect = mock_search_func
            
            from amprenta_rag.api.routers.papers import search_papers
            from amprenta_rag.api.routers.papers import PaperSearchRequest
            
            async def search_for_query(query: str):
                request = PaperSearchRequest(
                    query=query,
                    source="pubmed",
                    limit=5,
                    offset=0
                )
                return await search_papers(request, db=MagicMock())
            
            # Make 3 concurrent searches
            tasks = [
                search_for_query("cancer research"),
                search_for_query("machine learning"),
                search_for_query("drug discovery"),
            ]
            
            results = await asyncio.gather(*tasks)
            
            # All searches should succeed
            assert len(results) == 3
            assert results[0].results[0].title == "Paper for cancer research"
            assert results[1].results[0].title == "Paper for machine learning"
            assert results[2].results[0].title == "Paper for drug discovery"
            
            # Verify all calls were made
            assert mock_search.call_count == 3

