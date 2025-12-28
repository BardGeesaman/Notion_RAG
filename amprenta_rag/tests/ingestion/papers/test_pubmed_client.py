"""
Unit tests for PubMed client.

Tests search and metadata fetch functionality with mocked API responses.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from amprenta_rag.ingestion.papers.pubmed_client import PubMedRepository


# Sample XML responses for mocking
ESEARCH_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<eSearchResult>
    <IdList>
        <Id>12345678</Id>
        <Id>87654321</Id>
    </IdList>
</eSearchResult>
"""

EFETCH_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<PubmedArticleSet>
    <PubmedArticle>
        <MedlineCitation>
            <Article>
                <ArticleTitle>Sphingolipid metabolism in cancer</ArticleTitle>
                <Abstract>
                    <AbstractText>This is a test abstract about sphingolipids.</AbstractText>
                </Abstract>
                <AuthorList>
                    <Author>
                        <LastName>Smith</LastName>
                        <ForeName>John</ForeName>
                    </Author>
                    <Author>
                        <LastName>Doe</LastName>
                        <ForeName>Jane</ForeName>
                    </Author>
                </AuthorList>
                <Journal>
                    <Title>Journal of Lipid Research</Title>
                </Journal>
            </Article>
            <MeshHeadingList>
                <MeshHeading>
                    <DescriptorName>Sphingolipids</DescriptorName>
                </MeshHeading>
                <MeshHeading>
                    <DescriptorName>Neoplasms</DescriptorName>
                </MeshHeading>
            </MeshHeadingList>
        </MedlineCitation>
        <PubmedData>
            <ArticleIdList>
                <ArticleId IdType="pubmed">12345678</ArticleId>
                <ArticleId IdType="doi">10.1194/jlr.test</ArticleId>
                <ArticleId IdType="pmc">PMC123456</ArticleId>
            </ArticleIdList>
        </PubmedData>
    </PubmedArticle>
</PubmedArticleSet>
"""


class TestPubMedRepository:
    """Tests for PubMedRepository class."""

    def test_initialization(self):
        """Test repository initializes correctly."""
        repo = PubMedRepository()
        assert repo is not None
        assert hasattr(repo, "api_key")

    @patch("amprenta_rag.ingestion.papers.pubmed_client.requests.get")
    def test_search_papers_success(self, mock_get):
        """Test successful paper search returns PMIDs."""
        # Mock API response
        mock_response = MagicMock()
        mock_response.content = ESEARCH_RESPONSE.encode()
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        repo = PubMedRepository()
        pmids = repo.search_papers("sphingolipid cancer", max_results=10)

        assert len(pmids) == 2
        assert "12345678" in pmids
        assert "87654321" in pmids
        mock_get.assert_called_once()

    @patch("amprenta_rag.ingestion.papers.pubmed_client.requests.get")
    def test_search_papers_empty_results(self, mock_get):
        """Test search with no results returns empty list."""
        # Mock empty response
        empty_response = """<?xml version="1.0" encoding="UTF-8"?>
        <eSearchResult>
            <IdList></IdList>
        </eSearchResult>
        """
        mock_response = MagicMock()
        mock_response.content = empty_response.encode()
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        repo = PubMedRepository()
        pmids = repo.search_papers("nonexistent query")

        assert pmids == []

    @patch("amprenta_rag.ingestion.papers.pubmed_client.requests.get")
    def test_fetch_metadata_success(self, mock_get):
        """Test successful metadata fetch returns PaperMetadata."""
        # Mock API response
        mock_response = MagicMock()
        mock_response.content = EFETCH_RESPONSE.encode()
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        repo = PubMedRepository()
        metadata = repo.fetch_metadata("12345678")

        assert metadata is not None
        assert metadata.paper_id == "12345678"
        assert metadata.title == "Sphingolipid metabolism in cancer"
        assert metadata.abstract == "This is a test abstract about sphingolipids."
        assert len(metadata.authors) == 2
        assert "John Smith" in metadata.authors
        assert "Jane Doe" in metadata.authors
        assert metadata.journal == "Journal of Lipid Research"
        assert metadata.doi == "10.1194/jlr.test"
        assert metadata.pmid == "12345678"
        assert metadata.pmc_id == "PMC123456"
        assert len(metadata.mesh_terms) == 2
        assert "Sphingolipids" in metadata.mesh_terms
        assert "Neoplasms" in metadata.mesh_terms
        assert metadata.source == "pubmed"
        assert "pubmed.ncbi.nlm.nih.gov" in metadata.url

    @patch("amprenta_rag.ingestion.papers.pubmed_client.requests.get")
    def test_fetch_metadata_not_found(self, mock_get):
        """Test metadata fetch for non-existent paper returns None."""
        # Mock empty response
        empty_response = """<?xml version="1.0" encoding="UTF-8"?>
        <PubmedArticleSet></PubmedArticleSet>
        """
        mock_response = MagicMock()
        mock_response.content = empty_response.encode()
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        repo = PubMedRepository()
        metadata = repo.fetch_metadata("99999999")

        assert metadata is None

    def test_fetch_fulltext_not_implemented(self):
        """Test fulltext fetch returns None (not yet implemented)."""
        repo = PubMedRepository()
        fulltext = repo.fetch_fulltext("12345678")

        assert fulltext is None

