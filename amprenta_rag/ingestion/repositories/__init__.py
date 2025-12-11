"""
Repository implementations for public omics data repositories.

Provides implementations for:
- GEO (Gene Expression Omnibus) - Transcriptomics
- PRIDE (Proteomics Identifications Database) - Proteomics
- MetaboLights - Metabolomics (fallback)
- MW (Metabolomics Workbench) - Metabolomics (primary)
- ENA (European Nucleotide Archive) - Genomics (raw FASTQ files)
- ArrayExpress (via EBI BioStudies API) - Transcriptomics
"""

from __future__ import annotations

import os
import time

# Standard User-Agent header for all repository requests (Master Protocol compliance)
REPOSITORY_USER_AGENT = "ResearchBot/1.0 (Bioinformatics Data Pipeline)"
REPOSITORY_RATE_LIMIT_SECONDS = float(os.getenv("REPOSITORY_RATE_LIMIT_SECONDS", "0.34"))


def repo_rate_limit() -> None:
    """Global repository rate limiter (default ~3 req/sec)."""
    time.sleep(REPOSITORY_RATE_LIMIT_SECONDS)

from amprenta_rag.ingestion.repositories.base import RepositoryInterface
from amprenta_rag.ingestion.repositories.discovery import (
    discover_studies,
    fetch_study_metadata,
    get_repository,
)
from amprenta_rag.ingestion.repositories.arrayexpress import ArrayExpressRepository
from amprenta_rag.ingestion.repositories.ena import ENARepository
from amprenta_rag.ingestion.repositories.geo import GEORepository
from amprenta_rag.ingestion.repositories.metabolights import MetaboLightsRepository
from amprenta_rag.ingestion.repositories.mw import MWRepository
from amprenta_rag.ingestion.repositories.pride import PRIDERepository

# Alias for common misspelling
MetabolightsRepository = MetaboLightsRepository

__all__ = [
    "RepositoryInterface",
    "GEORepository",
    "PRIDERepository",
    "MetaboLightsRepository",
    "MetabolightsRepository",  # Alias
    "MWRepository",
    "ENARepository",
    "ArrayExpressRepository",
    "get_repository",
    "discover_studies",
    "fetch_study_metadata",
    "REPOSITORY_USER_AGENT",
]
