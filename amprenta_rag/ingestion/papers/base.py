"""
Base repository interface for scientific paper sources.

Defines the common interface that all paper repository implementations must follow.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class PaperMetadata:
    """
    Metadata for a scientific paper from any source.

    Attributes:
        paper_id: Repository-specific identifier (e.g., PMID, DOI)
        title: Paper title
        abstract: Paper abstract
        authors: List of author names
        journal: Journal name
        year: Publication year
        doi: Digital Object Identifier
        pmid: PubMed ID (if available)
        pmc_id: PubMed Central ID (if available)
        mesh_terms: MeSH terms for categorization
        source: Repository source (e.g., "pubmed", "biorxiv")
        url: URL to paper
    """

    paper_id: str
    title: str
    abstract: Optional[str] = None
    authors: List[str] = field(default_factory=list)
    journal: Optional[str] = None
    year: Optional[int] = None
    doi: Optional[str] = None
    pmid: Optional[str] = None
    pmc_id: Optional[str] = None
    mesh_terms: List[str] = field(default_factory=list)
    source: str = "unknown"
    url: Optional[str] = None


class PaperRepositoryInterface(ABC):
    """
    Base interface for all scientific paper repositories.

    All paper repository implementations (PubMed, PMC, bioRxiv, etc.) must
    implement these methods to provide a unified interface for discovery
    and paper ingestion.
    """

    @abstractmethod
    def search_papers(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search for papers matching a query.

        Args:
            query: Search query string
            max_results: Maximum number of results to return

        Returns:
            List of paper IDs matching the search criteria
        """

    @abstractmethod
    def fetch_metadata(self, paper_id: str) -> Optional[PaperMetadata]:
        """
        Fetch metadata for a specific paper.

        Args:
            paper_id: Repository-specific paper identifier

        Returns:
            PaperMetadata object, or None if paper not found
        """

    @abstractmethod
    def fetch_fulltext(self, paper_id: str) -> Optional[str]:
        """
        Fetch full text content for a specific paper.

        Args:
            paper_id: Repository-specific paper identifier

        Returns:
            Full text content as string, or None if not available
        """

