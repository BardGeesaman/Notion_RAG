"""
European Nucleotide Archive (ENA) repository implementation.

Provides ENA-specific implementation of the RepositoryInterface for
harvesting genomics (raw sequencing reads) data from ENA.

STRICT PROTOCOL: Uses ENA Browser API to avoid complex sra-toolkit configuration.
Following Master Bioinformatics Protocol for Genomics.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
from amprenta_rag.ingestion.repositories.base import RepositoryInterface
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.repository import DataFile, StudyMetadata

logger = get_logger(__name__)

# ENA Browser API endpoints (Master Protocol)
ENA_BASE_URL = "https://www.ebi.ac.uk/ena/portal/api"
ENA_SEARCH_URL = f"{ENA_BASE_URL}/search"
ENA_FIELDS_URL = f"{ENA_BASE_URL}/fields"


class ENARepository(RepositoryInterface):
    """
    European Nucleotide Archive (ENA) repository implementation.

    Supports discovery and harvesting of genomics studies (raw FASTQ sequencing files)
    from ENA Browser API.

    STRICT PROTOCOL COMPLIANCE:
    - Uses ENA Browser API (not NCBI SRA - avoids complex sra-toolkit configuration)
    - Provides direct fastq_ftp links
    - Does not download terabytes of data unless explicitly confirmed
    - Rate limiting: 1 second between requests
    """

    def __init__(self):
        """Initialize ENA repository."""
        # Rate limiting: 1 second between requests (Master Protocol)
        self._rate_limit_delay = 1.0

    def get_repository_name(self) -> str:
        """Get repository name."""
        return "ENA"

    def get_omics_type(self) -> str:
        """Get primary omics type."""
        return "genomics"

    def _rate_limit(self) -> None:
        """Respect ENA rate limits - 1 second between requests (Master Protocol)."""
        time.sleep(self._rate_limit_delay)

    def _make_request(
        self,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        timeout: int = 30,
    ) -> Optional[requests.Response]:
        """
        Make HTTP request with User-Agent header (Master Protocol).

        Args:
            url: Request URL
            params: Query parameters
            timeout: Request timeout in seconds

        Returns:
            Response object if successful, None otherwise
        """
        headers = {
            "User-Agent": REPOSITORY_USER_AGENT,
        }

        try:
            response = requests.get(url, params=params, headers=headers, timeout=timeout)

            # Check status code
            if response.status_code != 200:
                if response.status_code == 404:
                    logger.warning(
                        "[REPO][ENA] Resource not found (or private): %s",
                        url,
                    )
                else:
                    logger.warning(
                        "[REPO][ENA] Request failed with status %d: %s",
                        response.status_code,
                        url,
                    )
                return None

            return response

        except requests.exceptions.RequestException as e:
            logger.error("[REPO][ENA] Request error for %s: %r", url, e)
            return None

    def search_studies(
        self,
        keywords: List[str],
        filters: Optional[Dict[str, Any]] = None,
        max_results: int = 100,
    ) -> List[str]:
        """
        Search for ENA studies (read runs) matching keywords and filters.

        Uses ENA Browser API search endpoint.

        Args:
            keywords: List of search keywords
            filters: Optional filters (organism, library_strategy, library_source)
            max_results: Maximum number of results

        Returns:
            List of ENA study IDs (run_accession or study_accession)
        """
        logger.info(
            "[REPO][ENA] Searching for read runs with keywords: %s (max_results=%d)",
            keywords,
            max_results,
        )

        self._rate_limit()

        # Build search query
        # ENA API uses specific query syntax:
        # - For organism: scientific_name="Homo sapiens" or tax_eq(taxonomy_id)
        # - Keywords can be used in free-text search across fields

        query_parts = []

        # Handle keywords - if they look like organism names, use scientific_name field
        for keyword in keywords:
            keyword_clean = keyword.strip()
            if keyword_clean:
                # If keyword contains common organism patterns, use scientific_name
                if any(common_name in keyword_clean.lower() for common_name in
                       ['homo sapiens', 'mus musculus', 'rattus norvegicus', 'drosophila']):
                    query_parts.append(f'scientific_name="{keyword_clean}"')
                else:
                    # For general keywords, search in multiple fields
                    # ENA supports free-text search - use keyword as-is or in quotes
                    query_parts.append(f'"{keyword_clean}"')

        # Add filters to query
        if filters:
            if "organism" in filters:
                organism = filters["organism"]
                # Check if it's a taxonomy ID (numeric) or organism name
                if str(organism).isdigit():
                    query_parts.append(f'tax_eq({organism})')
                else:
                    query_parts.append(f'scientific_name="{organism}"')
            if "library_strategy" in filters:
                query_parts.append(f'library_strategy="{filters["library_strategy"]}"')
            if "library_source" in filters:
                query_parts.append(f'library_source="{filters["library_source"]}"')

        # Combine query parts with AND
        if query_parts:
            search_query = " AND ".join(query_parts)
        else:
            logger.warning("[REPO][ENA] No search query provided")
            return []

        # Use ENA Browser API search endpoint (Master Protocol)
        # Search for read_run result type
        params = {
            "result": "read_run",
            "query": search_query,
            "limit": min(max_results, 100),
            "fields": "run_accession,study_accession,experiment_accession,sample_accession,fastq_ftp,fastq_aspera,fastq_md5,read_count,base_count",
            "format": "json",
        }

        try:
            response = self._make_request(ENA_SEARCH_URL, params=params)

            if response is None:
                logger.error("[REPO][ENA] Search request failed")
                return []

            # Parse JSON response
            try:
                data = response.json()
            except ValueError as e:
                logger.error("[REPO][ENA] Error parsing JSON response: %r", e)
                return []

            # Extract run accessions (primary study IDs for ENA)
            study_ids = []
            if isinstance(data, list):
                for item in data:
                    run_accession = item.get("run_accession", "")
                    if run_accession:
                        study_ids.append(run_accession)
            elif isinstance(data, dict):
                # Handle different response structures
                if "data" in data:
                    for item in data["data"]:
                        run_accession = item.get("run_accession", "")
                        if run_accession:
                            study_ids.append(run_accession)

            logger.info(
                "[REPO][ENA] Search query '%s' returned %d results",
                search_query,
                len(study_ids),
            )

            return study_ids[:max_results]

        except Exception as e:
            logger.error("[REPO][ENA] Error searching read runs: %r", e)
            return []

    def fetch_study_metadata(self, study_id: str) -> Optional[StudyMetadata]:
        """
        Fetch metadata for a specific ENA read run.

        Args:
            study_id: ENA run accession (e.g., "ERR123456") or study accession

        Returns:
            StudyMetadata object, or None if study not found
        """
        logger.info("[REPO][ENA] Fetching metadata for read run %s", study_id)

        self._rate_limit()

        # Fetch run metadata
        params = {
            "result": "read_run",
            "query": f"run_accession={study_id}",
            "fields": "run_accession,study_accession,experiment_accession,sample_accession,secondary_study_accession,secondary_sample_accession,study_title,study_alias,experiment_title,experiment_alias,run_alias,tax_id,scientific_name,instrument_platform,instrument_model,library_strategy,library_source,library_selection,library_layout,read_count,base_count,fastq_ftp,fastq_aspera,fastq_md5,submitted_ftp,submitted_aspera,submitted_md5",
            "format": "json",
        }

        response = self._make_request(ENA_SEARCH_URL, params=params)

        if response is None:
            logger.error(
                "[REPO][ENA] Could not fetch metadata for run %s",
                study_id,
            )
            return None

        try:
            data = response.json()

            # Extract first result
            run_data = None
            if isinstance(data, list) and data:
                run_data = data[0]
            elif isinstance(data, dict):
                if "data" in data and isinstance(data["data"], list) and data["data"]:
                    run_data = data["data"][0]
                else:
                    run_data = data

            if not run_data:
                logger.warning("[REPO][ENA] No data found for run %s", study_id)
                return None

            # Extract fields
            title = run_data.get("study_title", "") or run_data.get("experiment_title", "") or f"ENA Run {study_id}"
            study_accession = run_data.get("study_accession", "")
            experiment_accession = run_data.get("experiment_accession", "")
            sample_accession = run_data.get("sample_accession", "")
            scientific_name = run_data.get("scientific_name", "")
            library_strategy = run_data.get("library_strategy", "")
            run_data.get("library_source", "")
            instrument_platform = run_data.get("instrument_platform", "")

            # Extract FASTQ FTP links (Master Protocol: generate links, don't download)
            fastq_ftp = run_data.get("fastq_ftp", "")
            if fastq_ftp:
                [link.strip() for link in fastq_ftp.split(";") if link.strip()]

            # Build StudyMetadata
            metadata = StudyMetadata(
                study_id=study_id,
                repository="ENA",
                title=title,
                summary=f"ENA read run: {study_id}\nStudy: {study_accession}\nExperiment: {experiment_accession}\nSample: {sample_accession}\nOrganism: {scientific_name}\nPlatform: {instrument_platform}",
                omics_type="genomics",
                organism=[scientific_name] if scientific_name else [],
                platform=f"{instrument_platform} ({library_strategy})" if instrument_platform else None,
                raw_metadata=run_data,
            )

            logger.info(
                "[REPO][ENA] Successfully fetched metadata for run %s: %s",
                study_id,
                title,
            )
            return metadata

        except Exception as e:
            logger.error(
                "[REPO][ENA] Error fetching metadata for run %s: %r",
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
        Fetch list of available FASTQ files for an ENA read run.

        Master Protocol: Generates FTP links but does not download terabytes of data.

        Args:
            study_id: ENA run accession
            file_types: Optional filter by file types (e.g., ["FASTQ"])

        Returns:
            List of DataFile objects with FTP download links
        """
        logger.info("[REPO][ENA] Fetching data files for run %s", study_id)

        data_files: List[DataFile] = []

        self._rate_limit()

        # Fetch run metadata to get FASTQ links
        metadata = self.fetch_study_metadata(study_id)
        if not metadata or not metadata.raw_metadata:
            return []

        run_data = metadata.raw_metadata

        # Extract FASTQ FTP links (Master Protocol)
        fastq_ftp = run_data.get("fastq_ftp", "")
        fastq_aspera = run_data.get("fastq_aspera", "")
        fastq_md5 = run_data.get("fastq_md5", "")

        if fastq_ftp:
            # Parse FTP links (semicolon-separated)
            ftp_links = [link.strip() for link in fastq_ftp.split(";") if link.strip()]
            md5_hashes = [h.strip() for h in fastq_md5.split(";") if fastq_md5 and h.strip()] if fastq_md5 else []

            for i, ftp_link in enumerate(ftp_links):
                filename = ftp_link.split("/")[-1] if "/" in ftp_link else ftp_link
                md5_hash = md5_hashes[i] if i < len(md5_hashes) else None

                # Create DataFile with both FTP and HTTP links
                # Pipeline Protocol: Convert FTP to HTTP for requests library
                http_url = f"http://{ftp_link}" if not ftp_link.startswith("http") else ftp_link
                ftp_url = f"ftp://{ftp_link}" if not ftp_link.startswith("ftp://") else ftp_link

                # Primary download URL: HTTP (for requests library compatibility)
                data_file = DataFile(
                    file_id=f"{study_id}_fastq_{i+1}",
                    filename=filename,
                    file_type="FASTQ",
                    download_url=http_url,  # HTTP for requests compatibility
                    size_bytes=None,  # ENA doesn't always provide size
                    description=f"FASTQ file {i+1} for run {study_id}" + (f" (MD5: {md5_hash})" if md5_hash else "") + f" | FTP: {ftp_url}",
                )
                data_files.append(data_file)

        # Also include Aspera links if available (faster transfer)
        if fastq_aspera:
            aspera_links = [link.strip() for link in fastq_aspera.split(";") if link.strip()]
            for i, aspera_link in enumerate(aspera_links):
                filename = aspera_link.split("/")[-1] if "/" in aspera_link else aspera_link
                data_file = DataFile(
                    file_id=f"{study_id}_fastq_aspera_{i+1}",
                    filename=filename,
                    file_type="FASTQ",
                    download_url=f"fasp://{aspera_link}" if not aspera_link.startswith("fasp://") else aspera_link,
                    size_bytes=None,
                    description=f"FASTQ file {i+1} (Aspera) for run {study_id}",
                )
                data_files.append(data_file)

        # Filter by file types if specified
        if file_types:
            data_files = [
                f for f in data_files
                if f.file_type.upper() in [ft.upper() for ft in file_types]
            ]

        logger.info(
            "[REPO][ENA] Found %d data files for run %s",
            len(data_files),
            study_id,
        )

        # Master Protocol: Log warning about large file sizes
        if data_files:
            logger.info(
                "[REPO][ENA] Note: FASTQ files can be very large (gigabytes to terabytes). "
                "FTP links generated but files not downloaded. Use download_data_file() to download specific files."
            )

        return data_files

