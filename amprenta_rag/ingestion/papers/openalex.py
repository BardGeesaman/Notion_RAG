"""
OpenAlex API client for searching and fetching scientific papers.

Uses OpenAlex REST API to search works, authors, and institutions.
"""

from __future__ import annotations

import os
import time
from typing import Dict, List, Optional

import requests

from amprenta_rag.ingestion.papers.base import PaperMetadata, PaperRepositoryInterface
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# OpenAlex API base URL
OPENALEX_API_BASE = "https://api.openalex.org"

class OpenAlexRepository(PaperRepositoryInterface):
    """
    OpenAlex repository client.

    Implements search and metadata fetch for scientific works from OpenAlex,
    an open catalog of scholarly papers, authors, and institutions.

    API Reference: https://docs.openalex.org/

    Example:
        >>> repo = OpenAlexRepository(email="your@email.com")
        >>> work_ids = repo.search_works("CRISPR gene editing", max_results=20)
        >>> metadata = repo.get_work(work_ids[0])
    """

    def __init__(self, email: Optional[str] = None):
        """
        Initialize OpenAlex repository client.

        Args:
            email: Email for polite pool access (gets better rate limits)
                   If not provided, uses OPENALEX_EMAIL env var
        """
        self.email = email or os.getenv("OPENALEX_EMAIL")
        self.headers = {"User-Agent": "Amprenta RAG (https://github.com/amprenta/rag)"}
        self._last_request_time = 0.0
        self.min_interval = 0.1  # 10 req/sec for polite pool
        
        if self.email:
            self.headers["mailto"] = self.email
            logger.info("[OPENALEX] Initialized with polite pool (email: %s)", self.email)
        else:
            logger.warning("[OPENALEX] No email configured - using basic rate limits. Set OPENALEX_EMAIL for better performance.")
    
    def _rate_limit(self) -> None:
        """Apply simple rate limiting with sleep."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last_request_time = time.time()
    
    def search_works(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search OpenAlex for works (papers) matching a query.

        Args:
            query: Search query string
            max_results: Maximum number of work IDs to return (max 200 per page)

        Returns:
            List of OpenAlex work IDs (URLs like https://openalex.org/W...)

        Raises:
            requests.RequestException: If API call fails
        """
        # OpenAlex uses pagination with per_page max 200
        per_page = min(max_results, 200)
        
        params = {
            "search": query,
            "per_page": per_page,
        }
        
        try:
            logger.debug("[OPENALEX] Searching works: %s (max_results=%d)", query, max_results)
            url = f"{OPENALEX_API_BASE}/works"
            
            self._rate_limit()
            response = requests.get(url, params=params, headers=self.headers, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            results = data.get("results", [])
            
            work_ids = [work.get("id") for work in results if work.get("id")]
            logger.info("[OPENALEX] Found %d works for query: %s", len(work_ids), query)
            return work_ids
        
        except requests.RequestException as e:
            logger.error("[OPENALEX] Search failed for query '%s': %r", query, e)
            raise
    
    # Alias for PaperRepositoryInterface compatibility
    def search_papers(self, query: str, max_results: int = 100) -> List[str]:
        """Alias for search_works to implement PaperRepositoryInterface."""
        return self.search_works(query, max_results)
    
    def get_work(self, work_id: str) -> Optional[PaperMetadata]:
        """
        Fetch metadata for a specific OpenAlex work.

        Args:
            work_id: OpenAlex work ID (URL like https://openalex.org/W123 or just W123)

        Returns:
            PaperMetadata object, or None if work not found
        """
        # Handle both full URLs and short IDs
        if not work_id.startswith("http"):
            work_id = f"https://openalex.org/{work_id}"
        
        try:
            logger.debug("[OPENALEX] Fetching work: %s", work_id)
            
            self._rate_limit()
            response = requests.get(work_id, headers=self.headers, timeout=30)
            
            if response.status_code == 404:
                logger.warning("[OPENALEX] Work not found: %s", work_id)
                return None
            
            response.raise_for_status()
            data = response.json()
            
            metadata = self._parse_openalex_work(data)
            logger.info("[OPENALEX] Fetched work %s: %s", work_id, metadata.title[:50] if metadata.title else "")
            return metadata
        
        except requests.RequestException as e:
            logger.error("[OPENALEX] Work fetch failed for %s: %r", work_id, e)
            return None
    
    # Alias for PaperRepositoryInterface compatibility
    def fetch_metadata(self, paper_id: str) -> Optional[PaperMetadata]:
        """Alias for get_work to implement PaperRepositoryInterface."""
        return self.get_work(paper_id)
    
    def fetch_fulltext(self, paper_id: str) -> Optional[str]:
        """
        Fetch full text for a paper.

        Note: OpenAlex API does not provide full text content directly.
        It may provide URLs to open access PDFs, but text extraction is not included.

        Args:
            paper_id: OpenAlex work ID

        Returns:
            None (full text not directly available)
        """
        logger.debug("[OPENALEX] Full text not available from OpenAlex API")
        return None
    
    def get_author(self, author_id: str) -> Optional[Dict]:
        """
        Fetch author information.

        Args:
            author_id: OpenAlex author ID (URL like https://openalex.org/A123 or just A123)

        Returns:
            Author metadata dictionary, or None if not found
        """
        if not author_id.startswith("http"):
            author_id = f"https://openalex.org/{author_id}"
        
        try:
            logger.debug("[OPENALEX] Fetching author: %s", author_id)
            
            self._rate_limit()
            response = requests.get(author_id, headers=self.headers, timeout=30)
            
            if response.status_code == 404:
                logger.warning("[OPENALEX] Author not found: %s", author_id)
                return None
            
            response.raise_for_status()
            return response.json()
        
        except requests.RequestException as e:
            logger.error("[OPENALEX] Author fetch failed for %s: %r", author_id, e)
            return None
    
    def get_institution(self, institution_id: str) -> Optional[Dict]:
        """
        Fetch institution information.

        Args:
            institution_id: OpenAlex institution ID (URL like https://openalex.org/I123 or just I123)

        Returns:
            Institution metadata dictionary, or None if not found
        """
        if not institution_id.startswith("http"):
            institution_id = f"https://openalex.org/{institution_id}"
        
        try:
            logger.debug("[OPENALEX] Fetching institution: %s", institution_id)
            
            self._rate_limit()
            response = requests.get(institution_id, headers=self.headers, timeout=30)
            
            if response.status_code == 404:
                logger.warning("[OPENALEX] Institution not found: %s", institution_id)
                return None
            
            response.raise_for_status()
            return response.json()
        
        except requests.RequestException as e:
            logger.error("[OPENALEX] Institution fetch failed for %s: %r", institution_id, e)
            return None
    
    def _parse_openalex_work(self, data: Dict) -> PaperMetadata:
        """Parse OpenAlex work data to PaperMetadata."""
        # Extract IDs
        work_id = data.get("id", "")
        doi = data.get("doi")
        if doi and doi.startswith("https://doi.org/"):
            doi = doi.replace("https://doi.org/", "")
        
        # Extract external IDs for PMID/PMC
        ids = data.get("ids", {}) or {}
        pmid = ids.get("pmid")
        if pmid and pmid.startswith("https://pubmed.ncbi.nlm.nih.gov/"):
            pmid = pmid.split("/")[-2]  # Extract number from URL
        
        pmc_id = ids.get("pmcid")
        if pmc_id and pmc_id.startswith("https://www.ncbi.nlm.nih.gov/pmc/articles/"):
            pmc_id = pmc_id.split("/")[-2]  # Extract PMC ID
        
        # Extract authors
        authorships = data.get("authorships", []) or []
        authors = []
        for authorship in authorships:
            author = authorship.get("author", {})
            display_name = author.get("display_name")
            if display_name:
                authors.append(display_name)
        
        # Extract publication info
        biblio = data.get("biblio", {}) or {}
        year = data.get("publication_year")
        
        # Extract journal/venue
        primary_location = data.get("primary_location", {}) or {}
        source = primary_location.get("source", {}) or {}
        journal = source.get("display_name")
        
        # Extract abstract (inverted index format)
        abstract_inverted = data.get("abstract_inverted_index")
        abstract = self._reconstruct_abstract(abstract_inverted) if abstract_inverted else None
        
        return PaperMetadata(
            paper_id=work_id,
            title=data.get("title", ""),
            abstract=abstract,
            authors=authors,
            journal=journal,
            year=year,
            doi=doi,
            pmid=pmid,
            pmc_id=pmc_id,
            mesh_terms=[],  # OpenAlex doesn't provide MeSH terms directly
            source="openalex",
            url=work_id,  # OpenAlex ID is the URL
        )
    
    def _reconstruct_abstract(self, inverted_index: Dict[str, List[int]]) -> str:
        """
        Reconstruct abstract text from OpenAlex inverted index format.

        Args:
            inverted_index: Dictionary mapping words to their positions

        Returns:
            Reconstructed abstract text
        """
        if not inverted_index:
            return ""
        
        # Build position -> word mapping
        word_positions = []
        for word, positions in inverted_index.items():
            for pos in positions:
                word_positions.append((pos, word))
        
        # Sort by position and join
        word_positions.sort(key=lambda x: x[0])
        abstract = " ".join([word for _, word in word_positions])
        
        return abstract

