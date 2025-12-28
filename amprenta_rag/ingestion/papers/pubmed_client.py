"""
PubMed/PMC client for searching and fetching scientific papers.

Uses NCBI E-utilities API to search PubMed and fetch paper metadata.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from typing import List, Optional

import requests

from amprenta_rag.ingestion.papers.base import PaperMetadata, PaperRepositoryInterface
from amprenta_rag.ingestion.papers.rate_limiter import ncbi_rate_limiter
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# NCBI E-utilities base URLs
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


class PubMedRepository(PaperRepositoryInterface):
    """
    PubMed/PMC repository client using NCBI E-utilities.

    Implements search and metadata fetch for scientific papers from PubMed.
    Rate-limited to respect NCBI usage guidelines (3 req/sec without API key,
    10 req/sec with NCBI_API_KEY).

    Example:
        >>> repo = PubMedRepository()
        >>> pmids = repo.search_papers("sphingolipid cancer", max_results=10)
        >>> metadata = repo.fetch_metadata(pmids[0])
        >>> print(metadata.title)
    """

    def __init__(self):
        """Initialize PubMed repository client."""
        self.api_key = os.getenv("NCBI_API_KEY")
        logger.info(
            "[PUBMED] Initialized PubMed client (API key: %s)",
            "present" if self.api_key else "absent",
        )

    @ncbi_rate_limiter
    def search_papers(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search PubMed for papers matching a query.

        Args:
            query: Search query string (uses PubMed search syntax)
            max_results: Maximum number of PMIDs to return

        Returns:
            List of PMIDs as strings

        Raises:
            requests.RequestException: If API call fails
        """
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "xml",
        }
        
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            logger.debug("[PUBMED] Searching: %s (max_results=%d)", query, max_results)
            response = requests.get(ESEARCH_URL, params=params, timeout=30)
            response.raise_for_status()

            # Parse XML response
            root = ET.fromstring(response.content)
            id_list = root.find("IdList")
            
            if id_list is None:
                logger.warning("[PUBMED] No IdList in search response")
                return []

            pmids = [id_elem.text for id_elem in id_list.findall("Id") if id_elem.text]
            logger.info("[PUBMED] Found %d PMIDs for query: %s", len(pmids), query)
            return pmids

        except requests.RequestException as e:
            logger.error("[PUBMED] Search failed for query '%s': %r", query, e)
            raise
        except ET.ParseError as e:
            logger.error("[PUBMED] Failed to parse search XML: %r", e)
            return []

    @ncbi_rate_limiter
    def fetch_metadata(self, paper_id: str) -> Optional[PaperMetadata]:
        """
        Fetch metadata for a specific PubMed paper.

        Args:
            paper_id: PubMed ID (PMID)

        Returns:
            PaperMetadata object, or None if paper not found or error occurs
        """
        params = {
            "db": "pubmed",
            "id": paper_id,
            "retmode": "xml",
        }
        
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            logger.debug("[PUBMED] Fetching metadata for PMID: %s", paper_id)
            response = requests.get(EFETCH_URL, params=params, timeout=30)
            response.raise_for_status()

            # Parse XML response
            root = ET.fromstring(response.content)
            article = root.find(".//PubmedArticle")
            
            if article is None:
                logger.warning("[PUBMED] No article found for PMID: %s", paper_id)
                return None

            # Extract metadata
            metadata = self._parse_pubmed_article(article, paper_id)
            logger.info("[PUBMED] Fetched metadata for PMID %s: %s", paper_id, metadata.title[:50])
            return metadata

        except requests.RequestException as e:
            logger.error("[PUBMED] Metadata fetch failed for PMID %s: %r", paper_id, e)
            return None
        except ET.ParseError as e:
            logger.error("[PUBMED] Failed to parse metadata XML: %r", e)
            return None

    def fetch_fulltext(self, paper_id: str) -> Optional[str]:
        """
        Fetch full text for a paper (if available in PMC Open Access).

        Args:
            paper_id: PubMed ID (PMID)

        Returns:
            Full text content, or None if not available

        Note:
            Full text fetching via PMC will be implemented in a later batch.
            This is a placeholder that returns None.
        """
        logger.debug("[PUBMED] Full text fetch not yet implemented for PMID: %s", paper_id)
        return None

    def _parse_pubmed_article(self, article: ET.Element, pmid: str) -> PaperMetadata:
        """
        Parse PubMed XML article element into PaperMetadata.

        Args:
            article: PubmedArticle XML element
            pmid: PubMed ID

        Returns:
            PaperMetadata object
        """
        # Title
        title_elem = article.find(".//ArticleTitle")
        title = title_elem.text if title_elem is not None and title_elem.text else "No title"

        # Abstract
        abstract_elem = article.find(".//AbstractText")
        abstract = abstract_elem.text if abstract_elem is not None else None

        # Authors
        authors = []
        author_list = article.find(".//AuthorList")
        if author_list is not None:
            for author in author_list.findall("Author"):
                last_name = author.find("LastName")
                fore_name = author.find("ForeName")
                if last_name is not None and last_name.text:
                    name = last_name.text
                    if fore_name is not None and fore_name.text:
                        name = f"{fore_name.text} {name}"
                    authors.append(name)

        # Journal
        journal_elem = article.find(".//Journal/Title")
        journal = journal_elem.text if journal_elem is not None else None

        # Publication year
        year = None
        pub_date = article.find(".//PubDate/Year")
        if pub_date is not None and pub_date.text:
            try:
                year = int(pub_date.text)
            except ValueError:
                pass

        # DOI
        doi = None
        article_id_list = article.find(".//ArticleIdList")
        if article_id_list is not None:
            for article_id in article_id_list.findall("ArticleId"):
                if article_id.get("IdType") == "doi" and article_id.text:
                    doi = article_id.text
                    break

        # PMC ID
        pmc_id = None
        if article_id_list is not None:
            for article_id in article_id_list.findall("ArticleId"):
                if article_id.get("IdType") == "pmc" and article_id.text:
                    pmc_id = article_id.text
                    break

        # MeSH terms
        mesh_terms = []
        mesh_heading_list = article.find(".//MeshHeadingList")
        if mesh_heading_list is not None:
            for mesh_heading in mesh_heading_list.findall("MeshHeading"):
                descriptor = mesh_heading.find("DescriptorName")
                if descriptor is not None and descriptor.text:
                    mesh_terms.append(descriptor.text)

        # URL
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

        return PaperMetadata(
            paper_id=pmid,
            title=title,
            abstract=abstract,
            authors=authors,
            journal=journal,
            year=year,
            doi=doi,
            pmid=pmid,
            pmc_id=pmc_id,
            mesh_terms=mesh_terms,
            source="pubmed",
            url=url,
        )

