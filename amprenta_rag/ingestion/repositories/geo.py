"""
Gene Expression Omnibus (GEO) repository implementation.

Provides GEO-specific implementation of the RepositoryInterface for
harvesting transcriptomics data from GEO.

STRICT PROTOCOL: Uses Biopython's Bio.Entrez module following NCBI guidelines.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional
from urllib.error import HTTPError

from Bio import Entrez

from amprenta_rag.ingestion.repositories.base import RepositoryInterface
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.repository import DataFile, StudyMetadata

logger = get_logger(__name__)

# GEO web URL for direct access
GEO_WEB_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


class GEORepository(RepositoryInterface):
    """
    Gene Expression Omnibus (GEO) repository implementation.

    Supports discovery and harvesting of transcriptomics studies from GEO.
    Uses Biopython's Bio.Entrez module following NCBI protocol guidelines.

    STRICT PROTOCOL COMPLIANCE:
    - Uses Bio.Entrez module (not requests)
    - Sets Entrez.email and Entrez.api_key globally
    - Follows search-then-fetch pattern
    - Respects NCBI rate limits
    - Handles 429 errors with retry
    """

    def __init__(self, api_key: Optional[str] = None, email: Optional[str] = None):
        """
        Initialize GEO repository.

        Args:
            api_key: Optional NCBI API key for higher rate limits
            email: Optional email address for NCBI (required by NCBI policy)
        """
        self._api_key = api_key
        self._email = email

        # Set global Entrez parameters (MANDATORY)
        if email:
            Entrez.email = email  # type: ignore[assignment]
            logger.debug("[REPO][GEO] Set Entrez.email = %s", email)
        else:
            logger.warning(
                "[REPO][GEO] No email provided. NCBI requires email for API access. "
                "Set NCBI_EMAIL in config or pass email parameter."
            )

        if api_key:
            Entrez.api_key = api_key  # type: ignore[assignment]
            self._rate_limit_delay = 0.1  # 10 requests/second with API key
            logger.debug("[REPO][GEO] Set Entrez.api_key and rate limit to 10 req/sec")
        else:
            self._rate_limit_delay = 0.34  # ~3 requests/second without API key
            logger.debug("[REPO][GEO] No API key, rate limit set to 3 req/sec")

    def get_repository_name(self) -> str:
        """Get repository name."""
        return "GEO"

    def get_omics_type(self) -> str:
        """Get primary omics type."""
        return "transcriptomics"

    def _rate_limit(self) -> None:
        """Respect NCBI rate limits."""
        time.sleep(self._rate_limit_delay)

    def _handle_api_error(self, error: Exception, retry_count: int = 0) -> bool:
        """
        Handle API errors with retry logic for rate limiting.

        Args:
            error: The exception that occurred
            retry_count: Current retry attempt

        Returns:
            True if should retry, False otherwise
        """
        if isinstance(error, HTTPError):
            if error.code == 429:  # Too Many Requests
                if retry_count < 3:
                    wait_time = 5 * (retry_count + 1)  # Exponential backoff
                    logger.warning(
                        "[REPO][GEO] Rate limit exceeded (429). Waiting %d seconds...",
                        wait_time,
                    )
                    time.sleep(wait_time)
                    return True
                else:
                    logger.error("[REPO][GEO] Rate limit exceeded. Max retries reached.")
                    return False
            else:
                logger.error(
                    "[REPO][GEO] HTTP error %d: %s",
                    error.code,
                    error.reason,
                )
                return False
        else:
            logger.error("[REPO][GEO] Error: %r", error)
            return False

    def _search_geo(
        self,
        query: str,
        max_results: int = 100,
        retry_count: int = 0,
    ) -> List[str]:
        """
        Search GEO using Entrez E-utilities (STRICT PROTOCOL).

        Uses the search-then-fetch pattern with Bio.Entrez module.

        Args:
            query: Entrez search query
            max_results: Maximum results to return
            retry_count: Current retry attempt (for error handling)

        Returns:
            List of GEO Series IDs (GSE numbers)
        """
        self._rate_limit()

        try:
            # Step 1: ESearch - Find study IDs
            logger.debug("[REPO][GEO] ESearch query: %s", query)
            handle = Entrez.esearch(
                db="gds",
                term=query,
                retmax=min(max_results, 10000),  # NCBI limit
                retmode="xml",
            )
            record = Entrez.read(handle)
            handle.close()

            id_list = record.get("IdList", [])
            if not id_list:
                logger.debug("[REPO][GEO] No IDs found for query: %s", query)
                return []

            logger.debug("[REPO][GEO] Found %d ID(s) from ESearch", len(id_list))

            # Step 2: ESummary - Get summaries to extract GSE accession numbers
            self._rate_limit()

            # Limit to max_results
            ids_to_fetch = id_list[:max_results]
            id_string = ",".join(ids_to_fetch)

            logger.debug("[REPO][GEO] ESummary for %d ID(s)", len(ids_to_fetch))
            handle = Entrez.esummary(db="gds", id=id_string, retmode="xml")
            summaries = Entrez.read(handle)
            handle.close()

            # Extract GSE IDs from summaries
            gse_ids: List[str] = []

            # Summaries is a list of dictionaries
            if isinstance(summaries, list):
                for summary in summaries:
                    accession = summary.get("Accession", "")
                    if accession and accession.startswith("GSE"):
                        gse_ids.append(accession)
            elif isinstance(summaries, dict):
                # Sometimes it's a dict keyed by ID
                for key, summary in summaries.items():
                    if isinstance(summary, dict):
                        accession = summary.get("Accession", "")
                        if accession and accession.startswith("GSE"):
                            gse_ids.append(accession)

            logger.info(
                "[REPO][GEO] Search query '%s' returned %d GSE IDs",
                query,
                len(gse_ids),
            )
            return gse_ids[:max_results]

        except HTTPError as e:
            if self._handle_api_error(e, retry_count):
                # Retry with increased retry count
                return self._search_geo(query, max_results, retry_count + 1)
            return []
        except Exception as e:
            logger.error("[REPO][GEO] Error searching GEO: %r", e)
            return []

    def search_studies(
        self,
        keywords: List[str],
        filters: Optional[Dict[str, Any]] = None,
        max_results: int = 100,
    ) -> List[str]:
        """
        Search for GEO studies matching keywords and filters.

        Args:
            keywords: List of search keywords
            filters: Optional filters (disease, organism, sample_type)
            max_results: Maximum number of results

        Returns:
            List of GEO Series IDs (e.g., ["GSE12345", "GSE67890"])
        """
        logger.info(
            "[REPO][GEO] Searching for studies with keywords: %s (max_results=%d)",
            keywords,
            max_results,
        )

        # Build Entrez query
        query_parts: List[str] = []

        # Add keywords
        for kw in keywords:
            query_parts.append(kw)

        # Add filters
        if filters:
            if "disease" in filters:
                disease = filters["disease"]
                query_parts.append(disease)

            if "organism" in filters:
                organism = filters["organism"]
                query_parts.append(f'"{organism}"[Organism]')

        # Combine with AND
        query = " AND ".join(query_parts)

        logger.debug("[REPO][GEO] Entrez query: %s", query)

        # Search GEO using strict protocol
        try:
            gse_ids = self._search_geo(query, max_results=max_results)
            return gse_ids
        except Exception as e:
            logger.error("[REPO][GEO] Error in search_studies: %r", e)
            return []

    def fetch_study_metadata(self, study_id: str) -> Optional[StudyMetadata]:
        """
        Fetch metadata for a specific GEO study (STRICT PROTOCOL).

        Uses search-then-fetch pattern with Bio.Entrez module.

        Args:
            study_id: GEO Series ID (e.g., "GSE12345")

        Returns:
            StudyMetadata object, or None if study not found
        """
        logger.info("[REPO][GEO] Fetching metadata for study %s", study_id)

        # Ensure study_id starts with GSE
        if not study_id.startswith("GSE"):
            study_id = f"GSE{study_id}"

        self._rate_limit()

        try:
            # Step 1: ESearch - Find the study by accession
            search_term = f"{study_id}[Accession]"
            logger.debug("[REPO][GEO] ESearch for: %s", search_term)

            handle = Entrez.esearch(db="gds", term=search_term, retmax=1, retmode="xml")
            record = Entrez.read(handle)
            handle.close()

            id_list = record.get("IdList", [])
            if not id_list:
                logger.warning("[REPO][GEO] Study %s not found", study_id)
                return None

            numeric_id = id_list[0]
            logger.debug("[REPO][GEO] Found numeric ID: %s", numeric_id)

            # Step 2: ESummary - Get summary metadata
            self._rate_limit()

            handle = Entrez.esummary(db="gds", id=numeric_id, retmode="xml")
            summaries = Entrez.read(handle)
            handle.close()

            # Extract summary - can be list or dict
            summary = None
            if isinstance(summaries, list) and len(summaries) > 0:
                summary = summaries[0]
            elif isinstance(summaries, dict):
                summary = summaries.get(numeric_id) or summaries.get(str(numeric_id))

            if not summary:
                logger.warning("[REPO][GEO] No summary found for %s", study_id)
                return None

            # Extract fields from summary
            # Note: Field names may vary, so we access safely
            title = summary.get("Title", summary.get("title", "")) or f"GEO Study {study_id}"
            summary_text = summary.get("Summary", summary.get("summary", "")) or ""
            organism = summary.get("Taxon", summary.get("taxon", "")) or ""
            platform = summary.get("Platform", summary.get("platform", "")) or ""
            samples_str = summary.get("Samples", summary.get("samples", "")) or ""

            # Parse samples count
            num_samples = None
            if samples_str and isinstance(samples_str, (str, int)):
                try:
                    num_samples = int(samples_str)
                except (ValueError, TypeError):
                    pass

            # Build StudyMetadata
            metadata = StudyMetadata(
                study_id=study_id,
                repository="GEO",
                title=title,
                summary=summary_text,
                omics_type="transcriptomics",
                organism=[organism] if organism else [],
                platform=platform if platform else None,
                num_samples=num_samples,
                raw_metadata=dict(summary),  # Store full summary as dict
            )

            # Try to extract disease from title/summary (heuristic)
            title_lower = title.lower()
            summary_lower = summary_text.lower()
            disease_keywords = [
                "als", "alzheimer", "parkinson", "cancer", "diabetes",
                "fragile x", "fxs", "autism", "depression",
            ]
            for disease_kw in disease_keywords:
                if disease_kw in title_lower or disease_kw in summary_lower:
                    # Map to canonical disease name
                    disease_map = {
                        "als": "ALS",
                        "alzheimer": "Alzheimer's disease",
                        "parkinson": "Parkinson's disease",
                        "fragile x": "Fragile X syndrome",
                        "fxs": "Fragile X syndrome",
                    }
                    canonical_disease = disease_map.get(disease_kw, disease_kw.title())
                    if canonical_disease not in metadata.disease:
                        metadata.disease.append(canonical_disease)

            logger.info(
                "[REPO][GEO] Successfully fetched metadata for study %s: %s",
                study_id,
                title,
            )
            return metadata

        except HTTPError as e:
            logger.error(
                "[REPO][GEO] HTTP error fetching metadata for study %s: %d %s",
                study_id,
                e.code,
                e.reason,
            )
            if e.code == 429:
                logger.warning("[REPO][GEO] Rate limited. Wait 5 seconds and retry.")
                time.sleep(5)
                # Retry once
                return self.fetch_study_metadata(study_id)
            return None
        except Exception as e:
            logger.error(
                "[REPO][GEO] Error fetching metadata for study %s: %r",
                study_id,
                e,
            )
            return None

    def fetch_study_data_files(
        self,
        study_id: str,
        file_types: Optional[List[str]] = None,
    ) -> List[DataFile]:
        """
        Fetch list of available data files for a GEO study.

        Args:
            study_id: GEO Series ID
            file_types: Optional filter by file types

        Returns:
            List of DataFile objects
        """
        logger.info("[REPO][GEO] Fetching data files for study %s", study_id)

        # Ensure study_id starts with GSE
        if not study_id.startswith("GSE"):
            study_id = f"GSE{study_id}"

        data_files: List[DataFile] = []

        # GEO provides data through FTP and web interfaces
        # Create a generic data file entry pointing to GEO's download page
        geo_download_url = f"{GEO_WEB_URL}?acc={study_id}&targ=self&form=text&view=full"

        data_file = DataFile(
            file_id=f"{study_id}_full",
            filename=f"{study_id}_full.txt",
            file_type="TXT",
            download_url=geo_download_url,
            description="GEO full metadata and data",
        )
        data_files.append(data_file)

        logger.info(
            "[REPO][GEO] Found %d data files for study %s",
            len(data_files),
            study_id,
        )
        return data_files


# Dashboard compatibility wrapper
def search_geo_studies(
    keyword: str = "",
    disease: str = "",
    species: str = "",
    omics_type: str = "",
    platform: str = "",
    max_results: int = 25,
) -> List[Dict[str, Any]]:
    """
    Lightweight GEO search returning a list of studies for dashboard selection.
    """
    repo = GEORepository()
    query_parts = []
    if keyword:
        query_parts.append(keyword)
    if disease:
        query_parts.append(disease)
    if species:
        query_parts.append(species)
    if platform:
        query_parts.append(platform)
    query = " AND ".join(query_parts) if query_parts else keyword or "GSE"

    try:
        gse_ids = repo._search_geo(query=query, max_results=max_results)
    except Exception as e:
        logger.error("[REPO][GEO] search_geo_studies failed: %r", e)
        return []

    results = []
    for gid in gse_ids:
        results.append(
            {
                "id": gid,
                "accession": gid,
                "title": f"GEO Series {gid}",
                "description": "",
                "disease": disease or "",
                "organism": species or "",
                "omics_type": omics_type or "transcriptomics",
            }
        )
        if len(results) >= max_results:
            break
    return results
