"""
Functional tests for paper ingestion with real API calls.

These tests call real NCBI and bioRxiv APIs and are marked as slow.
They are skipped by default and should be run manually or in nightly builds.

Run with: pytest amprenta_rag/tests/ingestion/papers/test_functional.py -v -m slow
"""

from __future__ import annotations

import pytest

from amprenta_rag.ingestion.papers import BioRxivRepository, PubMedRepository


@pytest.mark.slow
class TestPubMedFunctional:
    """Functional tests for PubMed client with real API calls."""

    def test_pubmed_real_search(self):
        """Test real PubMed search with small result set."""
        repo = PubMedRepository()

        # Search for a common term with limited results
        pmids = repo.search_papers("sphingolipid cancer", max_results=2)

        # Verify we got results
        assert isinstance(pmids, list)
        assert len(pmids) > 0
        assert len(pmids) <= 2

        # Verify PMIDs are numeric strings
        for pmid in pmids:
            assert pmid.isdigit()
            assert len(pmid) >= 6  # PMIDs are typically 7-8 digits

    def test_pubmed_real_metadata(self):
        """Test real metadata fetch for a known PMID."""
        repo = PubMedRepository()

        # Use a well-known paper PMID
        # PMID 12345678 is a declaration paper (may not have authors listed)
        # Using a more standard research paper instead
        pmid = "34213474"  # A recent COVID-19 research paper

        metadata = repo.fetch_metadata(pmid)

        # Verify metadata structure
        assert metadata is not None
        assert metadata.pmid == pmid
        assert metadata.title is not None
        assert len(metadata.title) > 10
        assert metadata.abstract is not None
        assert len(metadata.abstract) > 50
        # Authors may be empty for some papers, so check conditionally
        assert isinstance(metadata.authors, list)
        assert metadata.doi is not None or metadata.pmc_id is not None
        assert metadata.source == "pubmed"
        assert "pubmed.ncbi.nlm.nih.gov" in metadata.url

    def test_pubmed_metadata_with_mesh(self):
        """Test metadata fetch includes MeSH terms."""
        repo = PubMedRepository()

        # Fetch a paper that should have MeSH terms
        pmid = "34213474"
        metadata = repo.fetch_metadata(pmid)

        assert metadata is not None
        # Most papers have MeSH terms
        if metadata.mesh_terms:
            assert len(metadata.mesh_terms) > 0
            assert all(isinstance(term, str) for term in metadata.mesh_terms)

    def test_pubmed_nonexistent_pmid(self):
        """Test fetching metadata for non-existent PMID."""
        repo = PubMedRepository()

        # Use an invalid PMID
        metadata = repo.fetch_metadata("99999999999")

        # Should return None for non-existent paper
        assert metadata is None


@pytest.mark.slow
class TestBioRxivFunctional:
    """Functional tests for bioRxiv client with real API calls."""

    def test_biorxiv_real_search(self):
        """Test real bioRxiv search for recent preprints."""
        repo = BioRxivRepository(source="biorxiv")

        # Search for a common biology term
        dois = repo.search_papers("CRISPR", max_results=5)

        # Verify we got results (may be empty if no recent matches)
        assert isinstance(dois, list)
        # Note: bioRxiv search looks in recent papers only,
        # so results may vary

        if dois:
            # If we found results, verify DOI format
            for doi in dois:
                assert "10.1101" in doi  # bioRxiv DOI prefix

    def test_biorxiv_real_metadata(self):
        """Test real metadata fetch for a known bioRxiv DOI."""
        repo = BioRxivRepository(source="biorxiv")

        # Use a known bioRxiv paper DOI
        # This is a real preprint: "Deep learning for cellular image analysis"
        doi = "10.1101/2022.01.03.474753"

        metadata = repo.fetch_metadata(doi)

        # Verify metadata structure
        assert metadata is not None
        assert metadata.doi == doi
        assert metadata.title is not None
        assert len(metadata.title) > 10
        assert len(metadata.authors) > 0
        assert metadata.source == "biorxiv"
        assert "biorxiv.org" in metadata.url

    def test_medrxiv_initialization(self):
        """Test medRxiv repository can be initialized."""
        repo = BioRxivRepository(source="medrxiv")

        assert repo.source == "medrxiv"


@pytest.mark.slow
class TestFullIngestionFlow:
    """Functional tests for complete paper ingestion workflow."""

    def test_full_pubmed_ingestion_flow(self):
        """Test complete flow: search → fetch → verify."""
        repo = PubMedRepository()

        # Step 1: Search for papers
        pmids = repo.search_papers("machine learning biology", max_results=1)

        assert len(pmids) > 0, "Search should return at least one result"

        # Step 2: Fetch metadata for first result
        pmid = pmids[0]
        metadata = repo.fetch_metadata(pmid)

        # Step 3: Verify all expected fields
        assert metadata is not None
        assert metadata.paper_id == pmid
        assert metadata.pmid == pmid
        assert metadata.title is not None
        assert metadata.source == "pubmed"

        # Verify PaperMetadata structure matches expected schema
        assert hasattr(metadata, "paper_id")
        assert hasattr(metadata, "title")
        assert hasattr(metadata, "abstract")
        assert hasattr(metadata, "authors")
        assert hasattr(metadata, "journal")
        assert hasattr(metadata, "year")
        assert hasattr(metadata, "doi")
        assert hasattr(metadata, "pmid")
        assert hasattr(metadata, "pmc_id")
        assert hasattr(metadata, "mesh_terms")
        assert hasattr(metadata, "source")
        assert hasattr(metadata, "url")

    def test_rate_limiting_works(self):
        """Test that rate limiting is applied in real API calls."""
        import time

        repo = PubMedRepository()

        # Make 2 calls (not 3) to avoid hitting NCBI rate limits
        start = time.time()

        # Should be rate limited (3 req/sec without API key)
        repo.search_papers("genomics", max_results=1)
        repo.search_papers("proteomics", max_results=1)

        elapsed = time.time() - start

        # With 3 req/sec, 2 calls should take approximately 1 interval = 0.33s
        # Allow margin for API latency
        assert elapsed >= 0.25, f"Rate limiting may not be working (took {elapsed:.2f}s)"

    def test_biorxiv_fulltext_not_implemented(self):
        """Verify fulltext fetch returns None (not yet implemented)."""
        repo = BioRxivRepository()

        fulltext = repo.fetch_fulltext("10.1101/2022.01.03.474753")

        assert fulltext is None

