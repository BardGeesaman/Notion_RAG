"""
PRIDE Archive repository implementation.

Provides PRIDE-specific implementation of the RepositoryInterface for
harvesting proteomics data from PRIDE Archive.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.ingestion.repositories import (
    REPOSITORY_USER_AGENT,
    repo_rate_limit,
)
from amprenta_rag.ingestion.repositories.base import RepositoryInterface
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.repository import DataFile, StudyMetadata

logger = get_logger(__name__)

# PRIDE API endpoints (v2 - strict protocol)
# Base URL for PRIDE Archive REST API v2
PRIDE_BASE_URL = "https://www.ebi.ac.uk/pride/ws/archive/v2"
PRIDE_SEARCH_PROJECTS_URL = f"{PRIDE_BASE_URL}/search/projects"
PRIDE_PROJECT_DETAILS_URL = f"{PRIDE_BASE_URL}/projects/{{accession}}"
PRIDE_PROJECT_FILES_URL = f"{PRIDE_BASE_URL}/projects/{{accession}}/files"


class PRIDERepository(RepositoryInterface):
    """
    PRIDE Archive repository implementation.

    Supports discovery and harvesting of proteomics studies from PRIDE Archive.
    Uses PRIDE REST API v2.
    """

    def __init__(self):
        """Initialize PRIDE repository."""
        # Rate limiting: 1 second between requests (strict protocol)
        self._rate_limit_delay = 1.0

    def get_repository_name(self) -> str:
        """Get repository name."""
        return "PRIDE"

    def get_omics_type(self) -> str:
        """Get primary omics type."""
        return "proteomics"

    def _rate_limit(self) -> None:
        """Respect PRIDE rate limits - 1 second between requests (strict protocol)."""
        time.sleep(self._rate_limit_delay)

    def _make_request_with_retry(
        self,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        timeout: int = 30,
    ) -> Optional[requests.Response]:
        """
        Make HTTP request with retry logic for 500/503 errors.

        Args:
            url: Request URL
            params: Query parameters
            timeout: Request timeout in seconds

        Returns:
            Response object if successful, None otherwise
        """
        headers = {
            "Accept": "application/json",
            "User-Agent": REPOSITORY_USER_AGENT,
        }

        try:
            repo_rate_limit()
            response = requests.get(url, params=params, headers=headers, timeout=timeout)

            # Retry on 500/503 after 5 seconds
            if response.status_code in [500, 503]:
                logger.warning(
                    "[REPO][PRIDE] Received %d status code, retrying after 5 seconds...",
                    response.status_code,
                )
                time.sleep(5)
                repo_rate_limit()
                response = requests.get(url, params=params, headers=headers, timeout=timeout)

            # Check status code before parsing
            if response.status_code != 200:
                logger.warning(
                    "[REPO][PRIDE] Request failed with status %d: %s",
                    response.status_code,
                    url,
                )
                return None

            return response

        except requests.exceptions.RequestException as e:
            logger.error("[REPO][PRIDE] Request error for %s: %r", url, e)
            return None

    def search_studies(
        self,
        keywords: List[str],
        filters: Optional[Dict[str, any]] = None,
        max_results: int = 100,
    ) -> List[str]:
        """
        Search for PRIDE projects matching keywords and filters.

        Args:
            keywords: List of search keywords
            filters: Optional filters (disease, organism, sample_type)
            max_results: Maximum number of results

        Returns:
            List of PRIDE project IDs (e.g., ["PXD012345", "PXD67890"])
        """
        logger.info(
            "[REPO][PRIDE] Searching for projects with keywords: %s (max_results=%d)",
            keywords,
            max_results,
        )

        self._rate_limit()

        # Build search query - combine keywords and filters
        query_parts = keywords.copy()

        # Add filters to query
        if filters:
            if "disease" in filters:
                query_parts.append(filters["disease"])
            if "organism" in filters:
                query_parts.append(filters["organism"])

        # Combine keywords into single search term
        search_keyword = " ".join(query_parts)

        if not search_keyword:
            logger.warning("[REPO][PRIDE] No search keywords provided")
            return []

        self._rate_limit()  # Rate limit before request

        # Use strict protocol endpoint: GET /search/projects?keyword={term}
        params = {
            "keyword": search_keyword,
            "pageSize": min(max_results, 100),
            "sortDirection": "DESC",
        }

        try:
            response = self._make_request_with_retry(PRIDE_SEARCH_PROJECTS_URL, params=params)

            if response is None:
                logger.error("[REPO][PRIDE] Search request failed")
                return []

            # Safe JSON parsing
            try:
                data = response.json()
            except ValueError as e:
                logger.error("[REPO][PRIDE] Error parsing JSON response: %r", e)
                return []

            # Extract project IDs from response
            # Response structure may vary - handle both _embedded and direct list
            projects = []
            if "_embedded" in data:
                projects = data.get("_embedded", {}).get("projects", [])
            elif isinstance(data, list):
                projects = data
            elif "projects" in data:
                projects = data.get("projects", [])

            project_ids = [p.get("accession", "") for p in projects if p.get("accession")]
            project_ids = [pid for pid in project_ids if pid]  # Filter empty IDs

            logger.info(
                "[REPO][PRIDE] Search keyword '%s' returned %d results",
                search_keyword,
                len(project_ids),
            )

            # Paginate if needed
            if len(project_ids) < max_results:
                total_pages = data.get("page", {}).get("totalPages", 1)
                current_page = data.get("page", {}).get("number", 0)

                for page in range(current_page + 1, min(total_pages, (max_results // params["pageSize"]) + 1)):
                    self._rate_limit()  # Rate limit between pages

                    page_params = params.copy()
                    # Add pagination if API supports it
                    if "page" in data.get("page", {}):
                        page_params["page"] = page

                    page_response = self._make_request_with_retry(PRIDE_SEARCH_PROJECTS_URL, params=page_params)
                    if page_response is None:
                        break

                    try:
                        page_data = page_response.json()
                        page_projects = []
                        if "_embedded" in page_data:
                            page_projects = page_data.get("_embedded", {}).get("projects", [])
                        elif isinstance(page_data, list):
                            page_projects = page_data
                        elif "projects" in page_data:
                            page_projects = page_data.get("projects", [])

                        page_ids = [p.get("accession", "") for p in page_projects if p.get("accession")]
                        project_ids.extend([pid for pid in page_ids if pid])

                        if len(project_ids) >= max_results:
                            break
                    except Exception as e:
                        logger.warning("[REPO][PRIDE] Error parsing page %d: %r", page, e)
                        break

            return project_ids[:max_results]

        except Exception as e:
            logger.error("[REPO][PRIDE] Error searching projects: %r", e)
            return []

    def fetch_study_metadata(self, study_id: str) -> Optional[StudyMetadata]:
        """
        Fetch metadata for a specific PRIDE project.

        Args:
            study_id: PRIDE project ID (e.g., "PXD012345")

        Returns:
            StudyMetadata object, or None if project not found
        """
        logger.info("[REPO][PRIDE] Fetching metadata for project %s", study_id)

        # Ensure study_id starts with PXD
        if not study_id.startswith("PXD"):
            study_id = f"PXD{study_id}"

        self._rate_limit()  # Rate limit before request

        # Use strict protocol endpoint: GET /projects/{accession}
        url = PRIDE_PROJECT_DETAILS_URL.format(accession=study_id)

        response = self._make_request_with_retry(url)

        if response is None:
            logger.error(
                "[REPO][PRIDE] Could not fetch metadata for project %s",
                study_id,
            )
            return None

        # Safe JSON parsing
        try:
            data = response.json()
        except ValueError as e:
            logger.error(
                "[REPO][PRIDE] Error parsing JSON response for project %s: %r",
                study_id,
                e,
            )
            return None

        try:
            # Extract fields
            title = data.get("title", "") or f"PRIDE Project {study_id}"
            description = data.get("projectDescription", "") or ""

            # Extract publication info
            publications = data.get("publications", [])
            doi = None
            pubmed_id = None
            if publications:
                pub = publications[0]
                doi = pub.get("doi")
                pubmed_id = pub.get("pubmedId")

            # Extract organism
            organisms = data.get("organisms", [])
            organism_list = [org.get("name", "") for org in organisms if org.get("name")]

            # Extract sample processing protocol (may contain disease info)
            data.get("sampleProcessingProtocol", "") or ""

            # Extract instrument info
            instruments = data.get("instruments", [])
            platform = ", ".join([inst.get("name", "") for inst in instruments if inst.get("name")])

            # Extract submission date
            submission_date = data.get("submissionDate")
            publication_date = None
            if submission_date:
                try:
                    from datetime import datetime
                    publication_date = datetime.fromisoformat(submission_date.replace("Z", "+00:00"))
                except Exception as e:
                    logger.warning("[REPO][PRIDE] Error parsing submission date for %s: %r", study_id, e)

            # Try to extract disease from title/description (heuristic)
            title_lower = title.lower()
            desc_lower = description.lower()
            disease_keywords = [
                "als", "alzheimer", "parkinson", "cancer", "diabetes",
                "fragile x", "fxs", "autism", "depression",
            ]
            disease_list = []
            for disease_kw in disease_keywords:
                if disease_kw in title_lower or disease_kw in desc_lower:
                    disease_map = {
                        "als": "ALS",
                        "alzheimer": "Alzheimer's disease",
                        "parkinson": "Parkinson's disease",
                        "fragile x": "Fragile X syndrome",
                        "fxs": "Fragile X syndrome",
                    }
                    canonical_disease = disease_map.get(disease_kw, disease_kw.title())
                    if canonical_disease not in disease_list:
                        disease_list.append(canonical_disease)

            # Build StudyMetadata
            metadata = StudyMetadata(
                study_id=study_id,
                repository="PRIDE",
                title=title,
                summary=description,
                omics_type="proteomics",
                doi=doi,
                pubmed_id=pubmed_id,
                disease=disease_list,
                organism=organism_list,
                platform=platform if platform else None,
                publication_date=publication_date,
                raw_metadata=data,
            )

            logger.info(
                "[REPO][PRIDE] Successfully fetched metadata for project %s: %s",
                study_id,
                title,
            )
            return metadata

        except Exception as e:
            logger.error(
                "[REPO][PRIDE] Error fetching metadata for project %s: %r",
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
        Fetch list of available data files for a PRIDE project.

        Args:
            study_id: PRIDE project ID
            file_types: Optional filter by file types

        Returns:
            List of DataFile objects
        """
        logger.info("[REPO][PRIDE] Fetching data files for project %s", study_id)

        # Ensure study_id starts with PXD
        if not study_id.startswith("PXD"):
            study_id = f"PXD{study_id}"

        data_files: List[DataFile] = []

        self._rate_limit()  # Rate limit before request

        # Use strict protocol endpoint: GET /projects/{accession}/files
        url = PRIDE_PROJECT_FILES_URL.format(accession=study_id)

        try:
            response = self._make_request_with_retry(url)

            if response is None:
                logger.warning(
                    "[REPO][PRIDE] Could not fetch files for project %s",
                    study_id,
                )
                return []

            # Safe JSON parsing
            try:
                data = response.json()
            except ValueError as e:
                logger.error(
                    "[REPO][PRIDE] Error parsing JSON response for files: %r",
                    e,
                )
                return []

            # Extract files from response
            # PRIDE API v2 returns a list directly or in _embedded structure
            files = []
            if isinstance(data, list):
                files = data
            elif "_embedded" in data:
                files = data.get("_embedded", {}).get("files", [])
            elif "files" in data:
                files = data.get("files", [])

            for file_info in files:
                # PRIDE API v2 file structure
                file_id = file_info.get("accession", "") or file_info.get("id", "")
                filename = file_info.get("fileName", "") or file_info.get("name", "")
                file_category = file_info.get("fileCategory", {}).get("value", "") if isinstance(file_info.get("fileCategory"), dict) else file_info.get("fileCategory", "")
                file_type = file_category.upper() if file_category else file_info.get("fileType", "").upper()
                file_size = file_info.get("fileSizeBytes") or file_info.get("fileSize")

                # Get download URL from publicFileLocations
                download_url = ""
                file_locations = file_info.get("publicFileLocations", [])
                if file_locations:
                    download_url = file_locations[0].get("value", "") if isinstance(file_locations[0], dict) else file_locations[0]

                if not download_url:
                    download_url = file_info.get("downloadLink", "")

                # Filter by file types if specified
                # Also check filename extension for TSV/CSV files
                if file_types:
                    file_type_match = file_type in [ft.upper() for ft in file_types]
                    # Check filename extension as fallback
                    extension_match = any(filename.lower().endswith(f'.{ft.lower()}') for ft in file_types)
                    if not (file_type_match or extension_match):
                        continue

                # Prefer processed data files (e.g., protein tables)
                # Skip raw files unless specifically requested
                if file_type in ["RAW", "MZML", "MZXML"] and "raw" not in [ft.lower() for ft in (file_types or [])]:
                    continue

                data_file = DataFile(
                    file_id=str(file_id),
                    filename=filename,
                    file_type=file_type,
                    download_url=download_url,
                    size_bytes=file_size,
                    description=file_category if file_category else "",
                )
                data_files.append(data_file)

            # Paginate if needed (check if there are more pages)
            # Note: PRIDE API may return all files in one response or paginate
            # If we got 100 files, there might be more pages
            if len(files) == 100:
                # Try to fetch additional pages
                for page in range(1, 10):  # Limit to 10 pages
                    self._rate_limit()  # Rate limit between pages

                    page_url = f"{url}?page={page}"

                    try:
                        page_response = self._make_request_with_retry(page_url)

                        if page_response is None:
                            break

                        try:
                            page_data = page_response.json()
                        except ValueError as e:
                            logger.warning("[REPO][PRIDE] Error parsing page %d JSON: %r", page, e)
                            break

                        page_files = []
                        if isinstance(page_data, list):
                            page_files = page_data
                        elif "_embedded" in page_data:
                            page_files = page_data.get("_embedded", {}).get("files", [])
                        elif "files" in page_data:
                            page_files = page_data.get("files", [])

                        if not page_files:
                            break  # No more files

                        for file_info in page_files:
                            file_id = file_info.get("accession", "") or file_info.get("id", "")
                            filename = file_info.get("fileName", "") or file_info.get("name", "")
                            file_category = file_info.get("fileCategory", {}).get("value", "") if isinstance(file_info.get("fileCategory"), dict) else file_info.get("fileCategory", "")
                            file_type = file_category.upper() if file_category else file_info.get("fileType", "").upper()
                            file_size = file_info.get("fileSizeBytes") or file_info.get("fileSize")

                            download_url = ""
                            file_locations = file_info.get("publicFileLocations", [])
                            if file_locations:
                                download_url = file_locations[0].get("value", "") if isinstance(file_locations[0], dict) else file_locations[0]

                            if not download_url:
                                download_url = file_info.get("downloadLink", "")

                            # Filter by file types if specified
                            if file_types and file_type not in [ft.upper() for ft in file_types]:
                                continue

                            # Skip raw files unless specifically requested
                            if file_type in ["RAW", "MZML", "MZXML"] and "raw" not in [ft.lower() for ft in (file_types or [])]:
                                continue

                            data_file = DataFile(
                                file_id=str(file_id),
                                filename=filename,
                                file_type=file_type,
                                download_url=download_url,
                                size_bytes=file_size,
                                description=file_category if file_category else "",
                            )
                            data_files.append(data_file)
                    except Exception as e:
                        logger.warning("[REPO][PRIDE] Error fetching page %d: %r", page, e)
                        break

            logger.info(
                "[REPO][PRIDE] Found %d data files for project %s",
                len(data_files),
                study_id,
            )
            return data_files
        except Exception as e:
            logger.warning(
                "[REPO][PRIDE] Error fetching data files for project %s: %r",
                study_id,
                e,
            )
            return []


# Dashboard compatibility wrapper
def search_pride_studies(
    keyword: str = "",
    disease: str = "",
    species: str = "",
    omics_type: str = "",
    instrument: str = "",
    modality: str = "",
    max_results: int = 25,
) -> List[Dict[str, Any]]:
    try:
        repo = PRIDERepository()
        filters = {}
        if disease:
            filters["disease"] = disease
        if species:
            filters["organism"] = species
        projects = repo.search_studies(
            keywords=[kw for kw in [keyword] if kw],
            filters=filters,
            max_results=max_results,
        )
        results = []
        for pid in projects[:max_results]:
            results.append(
                {
                    "id": pid,
                    "accession": pid,
                    "title": f"PRIDE Project {pid}",
                    "description": "",
                    "disease": disease or "",
                    "organism": species or "",
                    "omics_type": "proteomics",
                }
            )
        return results
    except Exception as e:
        logger.error("[REPO][PRIDE] Search error: %r", e)
        return []

