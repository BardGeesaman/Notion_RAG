"""
Unified discovery across all omics repositories.

Provides functions to search for studies across multiple repositories
(GEO, PRIDE, MetaboLights, MW) with a unified interface.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional

from amprenta_rag.ingestion.repositories.base import RepositoryInterface
from amprenta_rag.ingestion.repositories.ena import ENARepository
from amprenta_rag.ingestion.repositories.geo import GEORepository
from amprenta_rag.ingestion.repositories.metabolights import MetaboLightsRepository
from amprenta_rag.ingestion.repositories.pride import PRIDERepository
from amprenta_rag.ingestion.repositories.mw import MWRepository
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.repository import StudyMetadata

logger = get_logger(__name__)

# Repository registry - will be populated as implementations are added
_REPOSITORY_REGISTRY: Dict[str, type[RepositoryInterface] | callable] = {
    "MW": MWRepository,
    "MW_LIPIDOMICS": lambda: MWRepository(omics_type="lipidomics"),
    "MW_METABOLOMICS": lambda: MWRepository(omics_type="metabolomics"),
    "GEO": GEORepository,
    "PRIDE": PRIDERepository,
    "MetaboLights": MetaboLightsRepository,
    "ENA": ENARepository,
}


def register_repository(name: str, repository_class: type[RepositoryInterface]) -> None:
    """
    Register a repository implementation.

    Args:
        name: Repository name (e.g., "GEO", "PRIDE")
        repository_class: Repository class that implements RepositoryInterface
    """
    _REPOSITORY_REGISTRY[name] = repository_class
    logger.info("[REPO][DISCOVERY] Registered repository: %s", name)


def get_repository(
    name: str,
    api_key: Optional[str] = None,
    email: Optional[str] = None,
) -> Optional[RepositoryInterface]:
    """
    Get a repository instance by name.

    Args:
        name: Repository name (e.g., "GEO", "PRIDE", "MW")
        api_key: Optional API key for repositories that support it (e.g., GEO/NCBI)
        email: Optional email for repositories that require it (e.g., GEO/NCBI)

    Returns:
        RepositoryInterface instance, or None if not found
    """
    repo_class = _REPOSITORY_REGISTRY.get(name)
    if repo_class is None:
        logger.warning("[REPO][DISCOVERY] Repository '%s' not found", name)
        return None

    # Handle callable factories (like MW_LIPIDOMICS)
    if callable(repo_class):
        return repo_class()

    # Handle regular classes
    if isinstance(repo_class, type):
        # Check if repository supports API key/email initialization
        if name == "GEO":
            from amprenta_rag.ingestion.repositories.geo import GEORepository
            return GEORepository(api_key=api_key, email=email)
        # Add other repositories with API key support here
        return repo_class()

    # Already an instance
    return repo_class


def list_available_repositories() -> List[str]:
    """
    List all available repository names.

    Returns:
        List of repository names
    """
    return list(_REPOSITORY_REGISTRY.keys())


def discover_studies(
    keywords: List[str],
    omics_type: Optional[str] = None,
    repository: Optional[str] = None,
    filters: Optional[Dict[str, any]] = None,
    max_results: int = 100,
) -> Dict[str, List[str]]:
    """
    Discover studies across repositories matching keywords and filters.

    Args:
        keywords: List of search keywords
        omics_type: Optional filter by omics type (e.g., "transcriptomics")
        repository: Optional specific repository name (e.g., "GEO")
        filters: Optional additional filters (disease, organism, sample_type)
        max_results: Maximum results per repository

    Returns:
        Dictionary mapping repository names to lists of study IDs
    """
    logger.info(
        "[REPO][DISCOVERY] Discovering studies with keywords: %s (omics_type=%s, repository=%s)",
        keywords,
        omics_type,
        repository,
    )

    results: Dict[str, List[str]] = {}

    # Determine which repositories to search
    repositories_to_search: List[str] = []

    if repository:
        # Search specific repository
        if repository in _REPOSITORY_REGISTRY:
            repositories_to_search = [repository]
        else:
            logger.warning(
                "[REPO][DISCOVERY] Repository '%s' not found, skipping",
                repository,
            )
            return results
    else:
        # Search all repositories, optionally filtered by omics type
        for repo_name in _REPOSITORY_REGISTRY.keys():
            repo_instance = get_repository(repo_name)
            if repo_instance:
                if omics_type is None or repo_instance.get_omics_type() == omics_type:
                    repositories_to_search.append(repo_name)

    # Load API keys and email from config if available
    try:
        from amprenta_rag.config import GEO_API_KEY, NCBI_EMAIL
        geo_api_key = GEO_API_KEY or None
        ncbi_email = NCBI_EMAIL or None
    except Exception:
        geo_api_key = os.getenv("GEO_API_KEY", "") or None
        ncbi_email = os.getenv("NCBI_EMAIL", "") or None

    # Search each repository
    for repo_name in repositories_to_search:
        # Pass API key and email for GEO repository
        if repo_name == "GEO":
            repo_instance = get_repository(repo_name, api_key=geo_api_key, email=ncbi_email)
        else:
            repo_instance = get_repository(repo_name)
        if not repo_instance:
            continue

        try:
            study_ids = repo_instance.search_studies(
                keywords=keywords,
                filters=filters,
                max_results=max_results,
            )
            results[repo_name] = study_ids

            logger.info(
                "[REPO][DISCOVERY] Found %d studies in %s",
                len(study_ids),
                repo_name,
            )
        except Exception as e:
            logger.error(
                "[REPO][DISCOVERY] Error searching repository %s: %r",
                repo_name,
                e,
            )
            results[repo_name] = []

    total_studies = sum(len(ids) for ids in results.values())
    logger.info(
        "[REPO][DISCOVERY] Total studies found across all repositories: %d",
        total_studies,
    )

    return results


def fetch_study_metadata(
    study_id: str,
    repository: str,
) -> Optional[StudyMetadata]:
    """
    Fetch metadata for a study from a specific repository.

    Args:
        study_id: Repository-specific study identifier
        repository: Repository name

    Returns:
        StudyMetadata object, or None if not found
    """
    # Load API keys and email from config if available
    geo_api_key = None
    ncbi_email = None
    if repository == "GEO":
        try:
            from amprenta_rag.config import GEO_API_KEY, NCBI_EMAIL
            geo_api_key = GEO_API_KEY or None
            ncbi_email = NCBI_EMAIL or None
        except Exception:
            geo_api_key = os.getenv("GEO_API_KEY", "") or None
            ncbi_email = os.getenv("NCBI_EMAIL", "") or None

    repo_instance = get_repository(repository, api_key=geo_api_key, email=ncbi_email)
    if not repo_instance:
        logger.warning(
            "[REPO][DISCOVERY] Repository '%s' not found",
            repository,
        )
        return None

    try:
        return repo_instance.fetch_study_metadata(study_id)
    except Exception as e:
        logger.error(
            "[REPO][DISCOVERY] Error fetching metadata for %s from %s: %r",
            study_id,
            repository,
            e,
        )
        return None

