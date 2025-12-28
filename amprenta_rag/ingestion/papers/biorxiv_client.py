"""
bioRxiv/medRxiv client for searching and fetching preprints.

Uses bioRxiv REST API to search and fetch preprint metadata.
"""

from __future__ import annotations

from datetime import datetime
from typing import List, Optional

import requests

from amprenta_rag.ingestion.papers.base import PaperMetadata, PaperRepositoryInterface
from amprenta_rag.ingestion.papers.rate_limiter import RateLimiter
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# bioRxiv API base URL
BIORXIV_API_BASE = "https://api.biorxiv.org"

# Rate limiter for bioRxiv (more lenient than NCBI)
biorxiv_rate_limiter = RateLimiter(requests_per_second=5.0)


class BioRxivRepository(PaperRepositoryInterface):
    """
    bioRxiv/medRxiv repository client.

    Implements search and metadata fetch for preprints from bioRxiv and medRxiv.
    bioRxiv is a preprint repository for biology, medRxiv for medicine.

    API Reference: https://api.biorxiv.org/

    Example:
        >>> repo = BioRxivRepository()
        >>> dois = repo.search_papers("CRISPR", max_results=10)
        >>> metadata = repo.fetch_metadata(dois[0])
        >>> print(metadata.title)
    """

    def __init__(self, source: str = "biorxiv"):
        """
        Initialize bioRxiv repository client.

        Args:
            source: "biorxiv" or "medrxiv"
        """
        if source not in ("biorxiv", "medrxiv"):
            raise ValueError(f"Invalid source: {source}. Must be 'biorxiv' or 'medrxiv'")
        self.source = source
        logger.info("[BIORXIV] Initialized %s client", source)

    @biorxiv_rate_limiter
    def search_papers(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search bioRxiv/medRxiv for papers.

        Note: The bioRxiv API doesn't have a keyword search endpoint.
        This implementation returns recent papers and filters by title/abstract.

        Args:
            query: Search query string
            max_results: Maximum number of DOIs to return

        Returns:
            List of DOIs

        Note:
            For production use, consider implementing a more sophisticated
            search using the collection API or external indexing.
        """
        try:
            # Get recent papers from last 30 days
            start_date = (datetime.now().replace(day=1)).strftime("%Y-%m-%d")
            end_date = datetime.now().strftime("%Y-%m-%d")

            url = f"{BIORXIV_API_BASE}/details/{self.source}/{start_date}/{end_date}"
            logger.debug("[BIORXIV] Fetching recent papers: %s to %s", start_date, end_date)

            response = requests.get(url, timeout=30)
            response.raise_for_status()

            data = response.json()
            collection = data.get("collection", [])

            # Filter by query in title or abstract
            query_lower = query.lower()
            matching_dois = []

            for paper in collection:
                title = paper.get("title", "").lower()
                abstract = paper.get("abstract", "").lower()

                if query_lower in title or query_lower in abstract:
                    doi = paper.get("doi")
                    if doi:
                        matching_dois.append(doi)

                if len(matching_dois) >= max_results:
                    break

            logger.info(
                "[BIORXIV] Found %d papers matching '%s'",
                len(matching_dois),
                query,
            )
            return matching_dois

        except requests.RequestException as e:
            logger.error("[BIORXIV] Search failed for query '%s': %r", query, e)
            return []

    @biorxiv_rate_limiter
    def fetch_metadata(self, paper_id: str) -> Optional[PaperMetadata]:
        """
        Fetch metadata for a specific paper by DOI.

        Args:
            paper_id: DOI (e.g., "10.1101/2024.01.01.123456")

        Returns:
            PaperMetadata object, or None if paper not found
        """
        try:
            # Clean DOI (remove prefix if present)
            doi = paper_id.replace("https://doi.org/", "").strip()

            url = f"{BIORXIV_API_BASE}/details/{self.source}/{doi}"
            logger.debug("[BIORXIV] Fetching metadata for DOI: %s", doi)

            response = requests.get(url, timeout=30)
            response.raise_for_status()

            data = response.json()
            collection = data.get("collection", [])

            if not collection:
                logger.warning("[BIORXIV] No paper found for DOI: %s", doi)
                return None

            # Get the latest version
            paper = collection[-1]  # Last item is latest version

            # Extract metadata
            title = paper.get("title", "No title")
            abstract = paper.get("abstract")
            authors = paper.get("authors", "").split(";")
            authors = [a.strip() for a in authors if a.strip()]

            # bioRxiv uses posted_date
            posted_date = paper.get("date")
            year = None
            if posted_date:
                try:
                    year = int(posted_date.split("-")[0])
                except (ValueError, IndexError):
                    pass

            # bioRxiv preprints don't have journal, use category
            category = paper.get("category", "")

            # URL
            url = f"https://www.biorxiv.org/content/{doi}"

            logger.info("[BIORXIV] Fetched metadata for DOI %s: %s", doi, title[:50])

            return PaperMetadata(
                paper_id=doi,
                title=title,
                abstract=abstract,
                authors=authors,
                journal=f"{self.source.title()} ({category})" if category else self.source.title(),
                year=year,
                doi=doi,
                source=self.source,
                url=url,
            )

        except requests.RequestException as e:
            logger.error("[BIORXIV] Metadata fetch failed for DOI %s: %r", paper_id, e)
            return None

    def fetch_fulltext(self, paper_id: str) -> Optional[str]:
        """
        Fetch full text for a paper.

        Args:
            paper_id: DOI

        Returns:
            Full text content, or None if not available

        Note:
            bioRxiv full text is available as PDF. Text extraction
            would require PDF parsing, which is deferred.
        """
        logger.debug("[BIORXIV] Full text fetch not yet implemented for DOI: %s", paper_id)
        return None

