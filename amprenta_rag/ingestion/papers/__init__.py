"""
Scientific paper ingestion from PubMed, PMC, bioRxiv, and other sources.
"""

from amprenta_rag.ingestion.papers.base import PaperMetadata, PaperRepositoryInterface
from amprenta_rag.ingestion.papers.rate_limiter import RateLimiter, ncbi_rate_limiter

__all__ = ["PaperMetadata", "PaperRepositoryInterface", "RateLimiter", "ncbi_rate_limiter"]

