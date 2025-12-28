"""
Scientific paper ingestion from PubMed, PMC, bioRxiv, and other sources.
"""

from amprenta_rag.ingestion.papers.base import PaperMetadata, PaperRepositoryInterface
from amprenta_rag.ingestion.papers.jats_parser import PaperContent, PaperSection, parse_jats_xml
from amprenta_rag.ingestion.papers.pubmed_client import PubMedRepository
from amprenta_rag.ingestion.papers.rate_limiter import RateLimiter, ncbi_rate_limiter

__all__ = [
    "PaperContent",
    "PaperMetadata",
    "PaperRepositoryInterface",
    "PaperSection",
    "PubMedRepository",
    "RateLimiter",
    "ncbi_rate_limiter",
    "parse_jats_xml",
]

