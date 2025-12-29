"""
Semantic Scholar API client for searching and fetching scientific papers.

Uses Semantic Scholar Graph API to search papers, fetch metadata, citations,
and get recommendations.
"""

from __future__ import annotations

import os
import time
from typing import Dict, List, Optional

import requests

from amprenta_rag.ingestion.papers.base import PaperMetadata, PaperRepositoryInterface
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Semantic Scholar API base URL
S2_API_BASE = "https://api.semanticscholar.org/graph/v1"

class SemanticScholarRepository(PaperRepositoryInterface):
    """
    Semantic Scholar repository client using Graph API.

    Implements search, metadata fetch, citations, and recommendations
    for scientific papers from Semantic Scholar's comprehensive database.

    API Reference: https://api.semanticscholar.org/

    Example:
        >>> repo = SemanticScholarRepository()
        >>> paper_ids = repo.search_papers("machine learning", max_results=10)
        >>> metadata = repo.fetch_metadata(paper_ids[0])
        >>> citations = repo.get_citations(paper_ids[0])
    """

    def __init__(self):
        """Initialize Semantic Scholar repository client."""
        self.api_key = os.getenv("S2_API_KEY")
        self.headers = {}
        self._last_request_time = 0.0
        
        # Set rate limit interval
        if self.api_key:
            self.headers["x-api-key"] = self.api_key
            self.min_interval = 0.1  # 10 req/sec with API key
        else:
            self.min_interval = 1.0  # 1 req/sec without API key
        
        logger.info(
            "[S2] Initialized Semantic Scholar client (API key: %s, rate: %.1f req/sec)",
            "present" if self.api_key else "absent",
            1.0 / self.min_interval,
        )
    
    def _rate_limit(self) -> None:
        """Apply simple rate limiting with sleep."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last_request_time = time.time()

    def search_papers(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search Semantic Scholar for papers matching a query.

        Args:
            query: Search query string
            max_results: Maximum number of paper IDs to return (max 100)

        Returns:
            List of Semantic Scholar paper IDs

        Raises:
            requests.RequestException: If API call fails
        """
        # Limit to API maximum
        limit = min(max_results, 100)
        
        params = {
            "query": query,
            "limit": limit,
            "fields": "paperId,title",
        }
        
        try:
            logger.debug("[S2] Searching: %s (max_results=%d)", query, limit)
            url = f"{S2_API_BASE}/paper/search"
            
            self._rate_limit()
            response = requests.get(url, params=params, headers=self.headers, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            papers = data.get("data", [])
            
            paper_ids = [p.get("paperId") for p in papers if p.get("paperId")]
            logger.info("[S2] Found %d papers for query: %s", len(paper_ids), query)
            return paper_ids
        
        except requests.RequestException as e:
            logger.error("[S2] Search failed for query '%s': %r", query, e)
            raise
    
    def fetch_metadata(self, paper_id: str) -> Optional[PaperMetadata]:
        """
        Fetch metadata for a specific Semantic Scholar paper.

        Args:
            paper_id: Semantic Scholar paper ID (can also accept DOI or PMID)

        Returns:
            PaperMetadata object, or None if paper not found
        """
        fields = [
            "paperId",
            "title",
            "abstract",
            "authors",
            "year",
            "venue",
            "citationCount",
            "influentialCitationCount",
            "externalIds",
            "url",
            "tldr",
        ]
        
        params = {"fields": ",".join(fields)}
        
        try:
            logger.debug("[S2] Fetching metadata for paper: %s", paper_id)
            url = f"{S2_API_BASE}/paper/{paper_id}"
            
            self._rate_limit()
            response = requests.get(url, params=params, headers=self.headers, timeout=30)
            
            if response.status_code == 404:
                logger.warning("[S2] Paper not found: %s", paper_id)
                return None
            
            response.raise_for_status()
            data = response.json()
            
            # Extract metadata
            metadata = self._parse_s2_paper(data)
            logger.info("[S2] Fetched metadata for %s: %s", paper_id, metadata.title[:50] if metadata.title else "")
            return metadata
        
        except requests.RequestException as e:
            logger.error("[S2] Metadata fetch failed for paper %s: %r", paper_id, e)
            return None
    
    def fetch_fulltext(self, paper_id: str) -> Optional[str]:
        """
        Fetch full text for a paper.

        Note: Semantic Scholar API does not provide full text content.
        This method returns None. Full text should be fetched from PMC or publisher.

        Args:
            paper_id: Semantic Scholar paper ID

        Returns:
            None (full text not available from Semantic Scholar)
        """
        logger.debug("[S2] Full text not available from Semantic Scholar API")
        return None
    
    def get_citations(self, paper_id: str, limit: int = 100) -> List[Dict]:
        """
        Get papers that cite this paper.

        Args:
            paper_id: Semantic Scholar paper ID
            limit: Maximum number of citations to return

        Returns:
            List of citing paper metadata dictionaries
        """
        params = {
            "fields": "paperId,title,year,authors,citationCount",
            "limit": min(limit, 1000),  # API max
        }
        
        try:
            logger.debug("[S2] Fetching citations for paper: %s", paper_id)
            url = f"{S2_API_BASE}/paper/{paper_id}/citations"
            
            self._rate_limit()
            response = requests.get(url, params=params, headers=self.headers, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            citations = data.get("data", [])
            
            # Extract citing papers
            citing_papers = []
            for item in citations:
                citing_paper = item.get("citingPaper", {})
                if citing_paper:
                    citing_papers.append(citing_paper)
            
            logger.info("[S2] Found %d citations for paper: %s", len(citing_papers), paper_id)
            return citing_papers
        
        except requests.RequestException as e:
            logger.error("[S2] Citations fetch failed for paper %s: %r", paper_id, e)
            return []
    
    def get_references(self, paper_id: str, limit: int = 100) -> List[Dict]:
        """
        Get papers referenced by this paper.

        Args:
            paper_id: Semantic Scholar paper ID
            limit: Maximum number of references to return

        Returns:
            List of referenced paper metadata dictionaries
        """
        params = {
            "fields": "paperId,title,year,authors,citationCount",
            "limit": min(limit, 1000),  # API max
        }
        
        try:
            logger.debug("[S2] Fetching references for paper: %s", paper_id)
            url = f"{S2_API_BASE}/paper/{paper_id}/references"
            
            self._rate_limit()
            response = requests.get(url, params=params, headers=self.headers, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            references = data.get("data", [])
            
            # Extract cited papers
            cited_papers = []
            for item in references:
                cited_paper = item.get("citedPaper", {})
                if cited_paper:
                    cited_papers.append(cited_paper)
            
            logger.info("[S2] Found %d references for paper: %s", len(cited_papers), paper_id)
            return cited_papers
        
        except requests.RequestException as e:
            logger.error("[S2] References fetch failed for paper %s: %r", paper_id, e)
            return []
    
    def get_recommendations(self, paper_id: str, limit: int = 10) -> List[str]:
        """
        Get recommended papers based on this paper.

        Args:
            paper_id: Semantic Scholar paper ID
            limit: Maximum number of recommendations to return

        Returns:
            List of recommended paper IDs
        """
        params = {
            "fields": "paperId,title",
            "limit": min(limit, 100),
        }
        
        try:
            logger.debug("[S2] Fetching recommendations for paper: %s", paper_id)
            url = f"{S2_API_BASE}/paper/{paper_id}/recommendations"
            
            self._rate_limit()
            response = requests.get(url, params=params, headers=self.headers, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            recommendations = data.get("recommendedPapers", [])
            
            paper_ids = [p.get("paperId") for p in recommendations if p.get("paperId")]
            logger.info("[S2] Found %d recommendations for paper: %s", len(paper_ids), paper_id)
            return paper_ids
        
        except requests.RequestException as e:
            logger.error("[S2] Recommendations fetch failed for paper %s: %r", paper_id, e)
            return []
    
    def _parse_s2_paper(self, data: Dict) -> PaperMetadata:
        """Parse Semantic Scholar API response to PaperMetadata."""
        # Extract external IDs
        external_ids = data.get("externalIds", {}) or {}
        doi = external_ids.get("DOI")
        pmid = external_ids.get("PubMed")
        pmc_id = external_ids.get("PubMedCentral")
        
        # Extract authors
        authors_data = data.get("authors", []) or []
        authors = [a.get("name", "") for a in authors_data if a.get("name")]
        
        # Use TLDR as enhanced abstract if available
        abstract = data.get("abstract")
        tldr = data.get("tldr", {})
        if not abstract and tldr:
            abstract = tldr.get("text")
        
        paper_id = data.get("paperId", "")
        
        return PaperMetadata(
            paper_id=paper_id,
            title=data.get("title", ""),
            abstract=abstract,
            authors=authors,
            journal=data.get("venue"),
            year=data.get("year"),
            doi=doi,
            pmid=pmid,
            pmc_id=f"PMC{pmc_id}" if pmc_id else None,
            mesh_terms=[],  # S2 doesn't provide MeSH terms
            source="semantic_scholar",
            url=data.get("url") or f"https://www.semanticscholar.org/paper/{paper_id}",
        )

