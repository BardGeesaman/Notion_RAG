"""
MetaboLights repository implementation.

Provides MetaboLights-specific implementation of the RepositoryInterface for
harvesting metabolomics data from MetaboLights.
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

# MetaboLights API endpoints
METABOLIGHTS_BASE_URL = "https://www.ebi.ac.uk/metabolights"
METABOLIGHTS_API_BASE_URL = f"{METABOLIGHTS_BASE_URL}/ws"
METABOLIGHTS_STUDIES_URL = f"{METABOLIGHTS_API_BASE_URL}/studies"
METABOLIGHTS_STUDY_DETAILS_URL = f"{METABOLIGHTS_API_BASE_URL}/studies/{{study_id}}"
METABOLIGHTS_STUDY_FILES_URL = f"{METABOLIGHTS_API_BASE_URL}/studies/{{study_id}}/files"


class MetaboLightsRepository(RepositoryInterface):
    """
    MetaboLights repository implementation.

    Supports discovery and harvesting of metabolomics studies from MetaboLights.
    Uses MetaboLights REST API.
    """

    def __init__(self):
        """Initialize MetaboLights repository."""
        # Rate limiting: 1 second between requests (strict protocol)
        self._rate_limit_delay = 1.0

    def get_repository_name(self) -> str:
        """Get repository name."""
        return "MetaboLights"

    def get_omics_type(self) -> str:
        """Get primary omics type."""
        return "metabolomics"

    def _rate_limit(self) -> None:
        """Respect MetaboLights rate limits - 1 second between requests (strict protocol)."""
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
                    "[REPO][METABOLIGHTS] Received %d status code, retrying after 5 seconds...",
                    response.status_code,
                )
                time.sleep(5)
                repo_rate_limit()
                response = requests.get(url, params=params, headers=headers, timeout=timeout)

            # Check status code before parsing
            if response.status_code != 200:
                if response.status_code == 404:
                    logger.warning(
                        "[REPO][METABOLIGHTS] Study not found or is private: %s",
                        url,
                    )
                else:
                    logger.warning(
                        "[REPO][METABOLIGHTS] Request failed with status %d: %s",
                        response.status_code,
                        url,
                    )
                return None

            return response

        except requests.exceptions.RequestException as e:
            logger.error("[REPO][METABOLIGHTS] Request error for %s: %r", url, e)
            return None

    def search_studies(
        self,
        keywords: List[str],
        filters: Optional[Dict[str, Any]] = None,
        max_results: int = 100,
    ) -> List[str]:
        """
        Search for MetaboLights studies matching keywords and filters.

        Args:
            keywords: List of search keywords
            filters: Optional filters (disease, organism, sample_type)
            max_results: Maximum number of results

        Returns:
            List of MetaboLights study IDs (e.g., ["MTBLS123", "MTBLS456"])
        """
        logger.info(
            "[REPO][METABOLIGHTS] Searching for studies with keywords: %s (max_results=%d)",
            keywords,
            max_results,
        )

        self._rate_limit()

        try:
            # MetaboLights API: Get all studies and filter client-side
            # Note: MetaboLights API may require authentication or have different endpoints
            # Try the studies endpoint
            resp = requests.get(METABOLIGHTS_STUDIES_URL, timeout=60)

            # Check if we got a valid response
            if resp.status_code == 404:
                # Try alternative endpoint
                alt_url = f"{METABOLIGHTS_BASE_URL}/webservice/study/list"
                resp = requests.get(alt_url, timeout=60)

            resp.raise_for_status()

            try:
                data = resp.json()
            except Exception as e:
                logger.error(
                    "[REPO][METABOLIGHTS] Error parsing JSON: %r. Response: %s",
                    e,
                    resp.text[:200],
                )
                return []

            # Extract study list - API may return different structures
            studies = []
            if isinstance(data, list):
                studies = data
            elif isinstance(data, dict):
                # Check both 'content' and 'studies' keys
                studies = data.get("content", [])
                if not studies:
                    studies = data.get("studies", [])

            logger.info(
                "[REPO][METABOLIGHTS] Fetched %d total studies from MetaboLights",
                len(studies),
            )

            # Filter by keywords
            matching_studies: List[str] = []
            keywords_lower = [kw.lower() for kw in keywords]

            for study in studies:
                # Handle both dict and string formats
                if isinstance(study, str):
                    # If study is just an ID string, we need to fetch metadata to search
                    # For now, include all string IDs and filter later if needed
                    study_id = study
                    matching_studies.append(study_id)
                    if len(matching_studies) >= max_results:
                        break
                    continue
                elif isinstance(study, dict):
                    study_id = str(study.get("studyIdentifier", study.get("id", "")) or "")
                    if not study_id:
                        continue
                    title = str(study.get("title", "")).lower()
                    description = str(study.get("description", "")).lower()
                else:
                    continue

                # Check if any keyword matches
                matches = False
                for kw in keywords_lower:
                    if kw in title or kw in description:
                        matches = True
                        break

                if not matches:
                    continue

                # Apply additional filters if provided
                if filters:
                    # Note: MetaboLights study list doesn't include all metadata
                    # We'd need to fetch individual study details for full filtering
                    # For now, we'll do basic filtering
                    if "disease" in filters:
                        disease = str(filters["disease"]).lower()
                        if disease not in title and disease not in description:
                            continue

                    if "organism" in filters:
                        organism = str(filters["organism"]).lower()
                        if organism not in title and organism not in description:
                            continue

                matching_studies.append(study_id)

                if len(matching_studies) >= max_results:
                    break

            logger.info(
                "[REPO][METABOLIGHTS] Found %d matching studies",
                len(matching_studies),
            )
            return matching_studies[:max_results]

        except Exception as e:
            logger.error("[REPO][METABOLIGHTS] Error searching studies: %r", e)
            return []

    def fetch_study_metadata(self, study_id: str) -> Optional[StudyMetadata]:
        """
        Fetch metadata for a specific MetaboLights study.

        Args:
            study_id: MetaboLights study ID (e.g., "MTBLS123")

        Returns:
            StudyMetadata object, or None if study not found
        """
        logger.info("[REPO][METABOLIGHTS] Fetching metadata for study %s", study_id)

        # Ensure study_id format is correct (must start with MTBLS)
        study_id = study_id.upper()
        if not study_id.startswith("MTBLS"):
            logger.error(
                "[REPO][METABOLIGHTS] Invalid study ID format. Must start with 'MTBLS': %s",
                study_id,
            )
            return None

        self._rate_limit()  # Rate limit before request

        # Use strict protocol endpoint: GET /studies/{study_id}
        url = METABOLIGHTS_STUDY_DETAILS_URL.format(study_id=study_id)

        response = self._make_request_with_retry(url)

        if response is None:
            logger.error(
                "[REPO][METABOLIGHTS] Could not fetch metadata for study %s",
                study_id,
            )
            return None

        # Safe JSON parsing
        try:
            data = response.json()
        except ValueError as e:
            logger.error(
                "[REPO][METABOLIGHTS] Error parsing JSON response for study %s: %r",
                study_id,
                e,
            )
            return None

        try:
            # MetaboLights API returns nested structure
            # Extract from mtblsStudy and isaInvestigation
            mtbls_study = data.get("mtblsStudy", {})
            isa_investigation = data.get("isaInvestigation", {})

            # Extract title and description from isaInvestigation
            title = isa_investigation.get("title", "") or mtbls_study.get("title", "") or f"MetaboLights Study {study_id}"
            description = isa_investigation.get("description", "") or mtbls_study.get("description", "") or ""

            # Extract publication info from isaInvestigation
            publications = isa_investigation.get("publications", [])
            doi = None
            pubmed_id = None
            if publications:
                pub = publications[0] if isinstance(publications, list) else publications
                if isinstance(pub, dict):
                    doi = pub.get("doi") or pub.get("doiString")
                    pubmed_id = pub.get("pubMedID") or pub.get("pubmedId")

            # Extract organism from isaInvestigation (may be in studies)
            organism_list = []
            studies = isa_investigation.get("studies", [])
            if studies and isinstance(studies, list):
                for study in studies:
                    # Organisms may be in study design or samples
                    if isinstance(study, dict):
                        # Try different possible locations
                        study_char = study.get("characteristicCategories", [])
                        for char in study_char:
                            if isinstance(char, dict) and "organism" in str(char).lower():
                                organism_list.append(char.get("name", ""))

            # Extract sample types
            sample_type_list: List[str] = []
            # Sample types may be in ISA-Tab structure

            # Extract submission date
            submission_date = isa_investigation.get("submissionDate") or mtbls_study.get("submissionDate")
            publication_date = None
            if submission_date:
                try:
                    from datetime import datetime
                    publication_date = datetime.fromisoformat(submission_date.replace("Z", "+00:00"))
                except Exception as e:
                    logger.warning("[REPO][METABOLIGHTS] Failed to parse submission_date with timezone: %r", e)
                    try:
                        publication_date = datetime.fromisoformat(submission_date)
                    except Exception as e_inner:
                        logger.warning("[REPO][METABOLIGHTS] Failed to parse submission_date fallback: %r", e_inner)

            # Extract technology/platform (may be in assays or study design)
            platform = None
            # Platform information may be in ISA-Tab structure

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
                repository="MetaboLights",
                title=title,
                summary=description,
                omics_type="metabolomics",
                doi=doi,
                pubmed_id=pubmed_id,
                disease=disease_list,
                organism=organism_list,
                sample_type=sample_type_list,
                platform=platform,
                publication_date=publication_date,
                raw_metadata=data,
            )

            logger.info(
                "[REPO][METABOLIGHTS] Successfully fetched metadata for study %s: %s",
                study_id,
                title,
            )
            return metadata

        except Exception as e:
            logger.error(
                "[REPO][METABOLIGHTS] Error processing metadata for study %s: %r",
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
        Fetch list of available data files for a MetaboLights study.

        Note: The MetaboLights `/studies/{study_id}/files` endpoint currently
        returns 500 Internal Server Error (confirmed via API diagnostic).
        This is a known MetaboLights API limitation.

        Files can be accessed via:
        - HTTP/FTP URLs from study details (studyHttpUrl/studyFtpUrl)
        - ISA-Tab file structure (i_Investigation.txt, s_*.txt, a_*.txt, m_*.tsv)

        Args:
            study_id: MetaboLights study ID
            file_types: Optional filter by file types

        Returns:
            List of DataFile objects (empty list if files endpoint unavailable)
        """
        logger.info("[REPO][METABOLIGHTS] Fetching data files for study %s", study_id)

        # Ensure study_id format is correct
        study_id = study_id.upper()
        if not study_id.startswith("MTBLS"):
            logger.error(
                "[REPO][METABOLIGHTS] Invalid study ID format. Must start with 'MTBLS': %s",
                study_id,
            )
            return []

        data_files: List[DataFile] = []

        self._rate_limit()  # Rate limit before request

        # Use strict protocol endpoint: GET /studies/{study_id}/files
        url = METABOLIGHTS_STUDY_FILES_URL.format(study_id=study_id)

        try:
            response = self._make_request_with_retry(url)

            if response is None:
                logger.warning(
                    "[REPO][METABOLIGHTS] Files endpoint unavailable for study %s "
                    "(returns 500). Files must be accessed via ISA-Tab structure.",
                    study_id,
                )
                return []

            # Safe JSON parsing
            try:
                data = response.json()
            except ValueError as e:
                logger.error(
                    "[REPO][METABOLIGHTS] Error parsing JSON response for files: %r",
                    e,
                )
                return []

            # Extract files from response
            # Response structure may vary - handle different formats
            files = []
            if isinstance(data, list):
                files = data
            elif isinstance(data, dict):
                # Check common keys
                files = data.get("files", [])
                if not files:
                    files = data.get("study", {}).get("files", [])
                if not files and "content" in data:
                    files = data.get("content", [])

            for file_info in files:
                filename = file_info.get("name", "")
                file_type = file_info.get("type", "").upper()
                file_size = file_info.get("size")
                download_url = file_info.get("downloadLink", "")

                # Filter by file types if specified
                if file_types and file_type not in [ft.upper() for ft in file_types]:
                    continue

                # Prefer processed data files (e.g., metabolite tables)
                # Skip raw files unless specifically requested
                if file_type in ["RAW", "MZML", "MZXML"] and "raw" not in [ft.lower() for ft in (file_types or [])]:
                    continue

                data_file = DataFile(
                    file_id=filename,  # Use filename as ID
                    filename=filename,
                    file_type=file_type,
                    download_url=download_url,
                    size_bytes=file_size,
                    description=file_info.get("description", ""),
                )
                data_files.append(data_file)

            logger.info(
                "[REPO][METABOLIGHTS] Found %d data files for study %s",
                len(data_files),
                study_id,
            )
            return data_files
        except Exception as e:
            logger.warning(
                "[REPO][METABOLIGHTS] Error fetching data files for study %s: %r",
                study_id,
                e,
            )
            return []


def search_metabolights_studies(
    keyword: str = "",
    disease: str = "",
    species: str = "",
    omics_type: str = "",
    analytical_platform: str = "",
    matrix: str = "",
    max_results: int = 25,
) -> List[Dict[str, Any]]:
    try:
        repo = MetaboLightsRepository()
        filters = {}
        if disease:
            filters["disease"] = disease
        organism = species or ""
        if organism:
            filters["organism"] = organism
        if analytical_platform:
            filters["analytical_platform"] = analytical_platform
        if matrix:
            filters["matrix"] = matrix

        study_ids = repo.search_studies(
            keywords=[kw for kw in [keyword] if kw],
            filters=filters,
            max_results=max_results,
        )
        results = []
        for sid in study_ids[:max_results]:
            results.append(
                {
                    "id": sid,
                    "accession": sid,
                    "title": f"MetaboLights Study {sid}",
                    "description": "",
                    "disease": disease or "",
                    "organism": organism or "",
                    "omics_type": omics_type or "metabolomics",
                    "analytical_platform": analytical_platform or "",
                    "matrix": matrix or "",
                }
            )
        return results
    except Exception as e:
        logger.error("[REPO][MetaboLights] Search error: %r", e)
        return []


def extract_metabolights_metabolites_from_isa_tab(
    study_id: str,
    download_dir: "Path | None" = None,
) -> "set[str]":
    """
    Extract metabolite features from MetaboLights ISA-Tab files.

    Uses the study details to get HTTP URL and finds metabolite data files (m_*.tsv).
    Parses ISA-Tab format to extract metabolite identifiers.

    Args:
        study_id: MetaboLights study ID (e.g., "MTBLS1")
        download_dir: Directory for temporary files

    Returns:
        Set of normalized metabolite identifiers
    """
    from pathlib import Path
    import re
    import zlib


    if download_dir is None:
        download_dir = Path("/tmp")

    download_dir.mkdir(parents=True, exist_ok=True)

    metabolite_set: set[str] = set()

    # Ensure study_id format
    study_id = study_id.upper()
    if not study_id.startswith("MTBLS"):
        logger.error("[FEATURE-EXTRACT][METABOLIGHTS] Invalid study ID format: %s", study_id)
        return metabolite_set

    logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Extracting metabolites from %s", study_id)

    try:
        repo = MetaboLightsRepository()

        # Fetch study metadata to get HTTP URL
        metadata = repo.fetch_study_metadata(study_id)
        if not metadata or not metadata.raw_metadata:
            logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] Could not fetch metadata for %s", study_id)
            return metabolite_set

        # Get HTTP URL from study details
        mtbls_study = metadata.raw_metadata.get("mtblsStudy", {}) if isinstance(metadata.raw_metadata, dict) else {}
        http_url = str(mtbls_study.get("studyHttpUrl", "") or "")

        if not http_url:
            logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] No HTTP URL found for %s", study_id)
            return metabolite_set

        # Convert FTP to HTTPS
        if http_url.startswith("ftp://"):
            https_url = http_url.replace("ftp://", "https://")
        elif http_url.startswith("http://ftp."):
            https_url = http_url.replace("http://ftp.", "https://ftp.")
        else:
            https_url = http_url

        logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Study directory: %s", https_url)

        metabolite_file_url: str | None = None

        # Check investigation file to find actual file names
        inv_url = f"{https_url}/i_Investigation.txt"
        try:
            headers = {"User-Agent": REPOSITORY_USER_AGENT}
            inv_response = requests.get(inv_url, headers=headers, timeout=30)
            if inv_response.status_code == 200:
                inv_content = inv_response.text

                # Look for MAF files first (m_*_maf.tsv)
                maf_files = re.findall(r"(m_[^\s\t\n]*maf[^\s\t\n]*\.(?:tsv|txt))", inv_content, re.IGNORECASE)
                if maf_files:
                    metabolite_file_url = f"{https_url}/{maf_files[0]}"
                    logger.info(
                        "[FEATURE-EXTRACT][METABOLIGHTS] Found MAF file from investigation: %s",
                        maf_files[0],
                    )

                # Fall back to other m_ files
                if not metabolite_file_url:
                    m_files = re.findall(r"(m_[^\s\t\n]+\.(?:tsv|txt))", inv_content)
                    if m_files:
                        metabolite_file_url = f"{https_url}/{m_files[0]}"
                        logger.info(
                            "[FEATURE-EXTRACT][METABOLIGHTS] Found metabolite data file from investigation: %s",
                            m_files[0],
                        )

                # Last resort: assay files (likely sample metadata)
                if not metabolite_file_url:
                    a_files = re.findall(r"(a_[^\s\t\n]+metabolite[^\s\t\n]+\.txt)", inv_content, re.IGNORECASE)
                    if a_files:
                        metabolite_file_url = f"{https_url}/{a_files[0]}"
                        logger.warning(
                            "[FEATURE-EXTRACT][METABOLIGHTS] Using assay file (may contain sample metadata, not metabolite data): %s",
                            a_files[0],
                        )
        except Exception as e:
            logger.debug("[FEATURE-EXTRACT][METABOLIGHTS] Could not read investigation file: %r", e)

        # If not found, try common patterns
        if not metabolite_file_url:
            metabolite_file_patterns = [
                f"m_{study_id}_LC-MS_positive_maf.tsv",
                f"m_{study_id}_LC-MS_negative_maf.tsv",
                f"m_{study_id}_GC-MS_maf.tsv",
                f"m_{study_id}_NMR_maf.tsv",
                f"m_{study_id}_metabolite_profiling_mass_spectrometry.tsv",
                f"m_{study_id}_metabolite_profiling_NMR_spectroscopy.tsv",
            ]

            for pattern in metabolite_file_patterns:
                test_url = f"{https_url}/{pattern}"
                test_response = requests.head(test_url, timeout=10)
                if test_response.status_code == 200:
                    metabolite_file_url = test_url
                    logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Found metabolite file: %s", pattern)
                    break

        if not metabolite_file_url:
            logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] No metabolite data file found for %s", study_id)
            return metabolite_set

        logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Downloading metabolite file: %s", metabolite_file_url)

        headers = {"User-Agent": REPOSITORY_USER_AGENT}
        response = requests.get(metabolite_file_url, headers=headers, timeout=120, stream=True)
        response.raise_for_status()

        # Parse ISA-Tab metabolite file
        is_gzipped = metabolite_file_url.endswith(".gz") or metabolite_file_url.endswith(".gzip")
        decompressor = zlib.decompressobj(zlib.MAX_WBITS | 32) if is_gzipped else None

        buffer = b""
        lines_processed = 0
        max_rows = 50000
        header_skipped = False
        metabolite_column_idx: int | None = None

        from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name

        logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Streaming file and extracting metabolites...")

        for chunk in response.iter_content(chunk_size=8192):
            buffer += chunk

            if is_gzipped:
                try:
                    assert decompressor is not None
                    decompressed = decompressor.decompress(buffer)
                    buffer = decompressor.unconsumed_tail
                except Exception as e:
                    logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] Decompression error, continuing: %r", e)
                    continue
            else:
                decompressed = buffer
                buffer = b""

            if not decompressed:
                continue

            text = decompressed.decode("utf-8", errors="ignore")

            for line in text.split("\n"):
                if lines_processed >= max_rows:
                    break

                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")

                if not header_skipped:
                    header_skipped = True
                    # Find best metabolite identification column in metadata section
                    for i, col in enumerate(parts):
                        col_lower = col.lower().strip()

                        if col_lower in ["sample name", "sample_name"] or col_lower.startswith("sample_"):
                            break

                        if any(
                            keyword in col_lower
                            for keyword in [
                                "metabolite identification",
                                "metabolite_identification",
                                "refmet_name",
                                "metabolite name",
                                "metabolite_name",
                                "database identifier",
                                "database_identifier",
                                "chemical name",
                                "chemical_name",
                            ]
                        ):
                            metabolite_column_idx = i
                            break

                    if metabolite_column_idx is None:
                        metabolite_column_idx = 0
                    continue

                if metabolite_column_idx is None or metabolite_column_idx >= len(parts):
                    continue

                metabolite_id = parts[metabolite_column_idx].strip()
                if not metabolite_id:
                    continue

                normalized = normalize_metabolite_name(metabolite_id)
                if normalized and len(normalized) > 1:
                    metabolite_set.add(normalized)

                lines_processed += 1

            if lines_processed >= max_rows:
                break

        logger.info(
            "[FEATURE-EXTRACT][METABOLIGHTS] Extracted %d unique metabolites from %d rows",
            len(metabolite_set),
            lines_processed,
        )

        return metabolite_set

    except Exception as e:
        logger.error("[FEATURE-EXTRACT][METABOLIGHTS] Error extracting metabolites: %r", e)
        return metabolite_set

