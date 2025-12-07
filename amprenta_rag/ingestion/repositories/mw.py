"""
Metabolomics Workbench (MW) repository implementation.

Provides MW-specific implementation of the RepositoryInterface for
harvesting lipidomics and metabolomics data from Metabolomics Workbench.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.ingestion.repositories.base import RepositoryInterface
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.repository import DataFile, StudyMetadata

logger = get_logger(__name__)

MW_BASE_URL = "https://www.metabolomicsworkbench.org/rest"
MW_STUDY_SUMMARY_URL = f"{MW_BASE_URL}/study/study_id/ST/summary"


class MWRepository(RepositoryInterface):
    """
    Metabolomics Workbench repository implementation.
    
    Supports discovery and harvesting of lipidomics and metabolomics
    studies from Metabolomics Workbench.
    """
    
    def __init__(self, omics_type: str = "metabolomics"):
        """
        Initialize MW repository.
        
        Args:
            omics_type: Primary omics type ("metabolomics" or "lipidomics")
        """
        self._omics_type = omics_type
        self._study_summaries_cache: Optional[List[Dict[str, Any]]] = None
    
    def get_repository_name(self) -> str:
        """Get repository name."""
        return "MW"
    
    def get_omics_type(self) -> str:
        """Get primary omics type."""
        return self._omics_type
    
    def _fetch_all_study_summaries(self) -> List[Dict[str, Any]]:
        """
        Fetch all public study summaries from MW (with caching).
        
        Returns:
            List of study summary dictionaries
        """
        if self._study_summaries_cache is not None:
            return self._study_summaries_cache
        
        try:
            from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
            headers = {"User-Agent": REPOSITORY_USER_AGENT}
            resp = requests.get(MW_STUDY_SUMMARY_URL, params={"format": "json"}, headers=headers, timeout=60)
            resp.raise_for_status()
            data = resp.json()
            
            # Normalize data structure
            if isinstance(data, dict):
                studies = []
                for sid, summary in data.items():
                    if isinstance(summary, dict):
                        summary.setdefault("study_id", sid)
                        studies.append(summary)
                data = studies
            
            self._study_summaries_cache = data
            logger.info("[REPO][MW] Cached %d study summaries", len(data))
            return data
            
        except Exception as e:
            logger.error("[REPO][MW] Error fetching all study summaries: %r", e)
            raise
    
    def search_studies(
        self,
        keywords: List[str],
        filters: Optional[Dict[str, any]] = None,
        max_results: int = 100,
    ) -> List[str]:
        """
        Search for MW studies matching keywords and filters.
        
        Args:
            keywords: List of search keywords
            filters: Optional filters (disease, organism, sample_type)
            max_results: Maximum number of results
            
        Returns:
            List of MW study IDs (e.g., ["ST001111", "ST001112"])
        """
        logger.info(
            "[REPO][MW] Searching for studies with keywords: %s (max_results=%d)",
            keywords,
            max_results,
        )
        
        summaries = self._fetch_all_study_summaries()
        matching_ids: List[str] = []
        seen_ids: set[str] = set()
        
        # Search by keywords
        for kw in keywords:
            kw_lower = kw.lower()
            for study in summaries:
                study_id = study.get("study_id", "")
                if not study_id or study_id in seen_ids:
                    continue
                
                title = str(study.get("study_title", "") or "").lower()
                summary_text = str(study.get("study_summary", "") or "").lower()
                
                if kw_lower in title or kw_lower in summary_text:
                    matching_ids.append(study_id)
                    seen_ids.add(study_id)
                    
                    if len(matching_ids) >= max_results:
                        break
            
            if len(matching_ids) >= max_results:
                break
        
        # Apply filters
        if filters:
            filtered_ids: List[str] = []
            for study_id in matching_ids:
                # Find study in summaries
                study = next((s for s in summaries if s.get("study_id") == study_id), None)
                if not study:
                    continue
                
                # Check disease filter
                if "disease" in filters:
                    study_disease = str(study.get("disease", "") or "").lower()
                    filter_disease = str(filters["disease"]).lower()
                    if filter_disease not in study_disease:
                        continue
                
                # Check organism filter
                if "organism" in filters:
                    study_organism = str(study.get("organism", "") or "").lower()
                    filter_organism = str(filters["organism"]).lower()
                    if filter_organism not in study_organism:
                        continue
                
                # Check sample_type filter
                if "sample_type" in filters:
                    study_sample_type = str(study.get("sample_type", "") or "").lower()
                    filter_sample_type = str(filters["sample_type"]).lower()
                    if filter_sample_type not in study_sample_type:
                        continue
                
                filtered_ids.append(study_id)
            
            matching_ids = filtered_ids[:max_results]
        
        logger.info(
            "[REPO][MW] Found %d matching studies",
            len(matching_ids),
        )
        return matching_ids[:max_results]
    
    def fetch_study_metadata(self, study_id: str) -> Optional[StudyMetadata]:
        """
        Fetch metadata for a specific MW study.
        
        Args:
            study_id: MW study ID (e.g., "ST001111")
            
        Returns:
            StudyMetadata object, or None if study not found
        """
        logger.info("[REPO][MW] Fetching metadata for study %s", study_id)
        
        try:
            from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
            headers = {"User-Agent": REPOSITORY_USER_AGENT}
            
            # Fetch study summary
            url = f"{MW_BASE_URL}/study/study_id/{study_id}/summary"
            resp = requests.get(url, headers=headers, timeout=30)
            resp.raise_for_status()
            
            content_type = resp.headers.get("content-type", "")
            if "application/json" in content_type:
                data = resp.json()
                if isinstance(data, list) and data:
                    summary = data[0]
                elif isinstance(data, dict) and "summary" in data:
                    summary = data["summary"]
                else:
                    summary = data if isinstance(data, dict) else {}
            else:
                try:
                    summary = resp.json()
                except:
                    logger.warning("[REPO][MW] Study summary for %s not in JSON format", study_id)
                    return None
            
            # Extract metadata
            title = summary.get("study_title", "") or f"MW Study {study_id}"
            summary_text = summary.get("study_summary", "") or ""
            doi = summary.get("doi", "")
            disease = summary.get("disease", "")
            organism = summary.get("organism", "")
            sample_type = summary.get("sample_type", "")
            
            # Build StudyMetadata
            metadata = StudyMetadata(
                study_id=study_id,
                repository="MW",
                title=title,
                summary=summary_text,
                omics_type=self._omics_type,
                doi=doi if doi else None,
                disease=[disease] if disease else [],
                organism=[organism] if organism else [],
                sample_type=[sample_type] if sample_type else [],
                raw_metadata=summary,
            )
            
            # Try to fetch mwTab to get more details
            try:
                mwtab_url = f"{MW_BASE_URL}/study/study_id/{study_id}/mwtab"
                mwtab_resp = requests.get(mwtab_url, headers=headers, timeout=30)
                if mwtab_resp.status_code == 200:
                    # mwTab contains structured metadata
                    # Could parse for additional fields if needed
                    pass
            except:
                pass  # mwTab fetch is optional
            
            logger.info(
                "[REPO][MW] Successfully fetched metadata for study %s: %s",
                study_id,
                title,
            )
            return metadata
            
        except Exception as e:
            logger.error(
                "[REPO][MW] Error fetching metadata for study %s: %r",
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
        Fetch list of available data files for an MW study.
        
        Args:
            study_id: MW study ID
            file_types: Optional filter by file types
            
        Returns:
            List of DataFile objects
        """
        logger.info("[REPO][MW] Fetching data files for study %s", study_id)
        
        data_files: List[DataFile] = []
        
        try:
            from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
            headers = {"User-Agent": REPOSITORY_USER_AGENT}
            
            # MW provides data through various endpoints
            # Try to get file listing if available
            url = f"{MW_BASE_URL}/study/study_id/{study_id}/data"
            resp = requests.get(url, headers=headers, timeout=30)
            
            if resp.status_code == 200:
                # MW may return file information in various formats
                # For now, we'll create a generic data file entry
                data_file = DataFile(
                    file_id=f"{study_id}_data",
                    filename=f"{study_id}_data.txt",
                    file_type="TXT",
                    download_url=url,
                    description="MW study data",
                )
                data_files.append(data_file)
            
            # Also try mwTab endpoint
            mwtab_url = f"{MW_BASE_URL}/study/study_id/{study_id}/mwtab"
            mwtab_resp = requests.get(mwtab_url, headers=headers, timeout=30)
            if mwtab_resp.status_code == 200:
                mwtab_file = DataFile(
                    file_id=f"{study_id}_mwtab",
                    filename=f"{study_id}.mwtab",
                    file_type="MWTAB",
                    download_url=mwtab_url,
                    description="MW mwTab formatted data",
                )
                data_files.append(mwtab_file)
            
            # Filter by file types if specified
            if file_types:
                data_files = [
                    f for f in data_files
                    if f.file_type.upper() in [ft.upper() for ft in file_types]
                ]
            
            logger.info(
                "[REPO][MW] Found %d data files for study %s",
                len(data_files),
                study_id,
            )
            return data_files
            
        except Exception as e:
            logger.warning(
                "[REPO][MW] Error fetching data files for study %s: %r",
                study_id,
                e,
            )
            return []


# Convenience function for dashboard compatibility
def search_mw_studies(
    keyword: str = "",
    disease: str = "",
    species: str = "",
    omics_type: str = "",
    matrix: str = "",
    platform: str = "",
    analytical_platform: str = "",
    max_results: int = 25,
) -> List[Dict[str, Any]]:
    """
    Search MW studies with filters (dashboard-compatible interface).
    
    Args:
        keyword: Search keyword
        disease: Disease filter
        species: Species/organism filter
        omics_type: Omics type filter
        matrix: Matrix/sample type filter
        platform: Analytical platform filter (alias)
        analytical_platform: Analytical platform filter
        max_results: Maximum results to return
        
    Returns:
        List of study dictionaries with id, title, description, etc.
    """
    repo = MWRepository(omics_type=omics_type if omics_type else "metabolomics")
    
    # Fetch all study summaries
    try:
        summaries = repo._fetch_all_study_summaries()
    except Exception as e:
        logger.error("[REPO][MW] Failed to fetch study summaries: %r", e)
        return []
    
    # Filter studies
    results = []
    keyword_lower = keyword.lower().strip() if keyword else ""
    disease_lower = disease.lower().strip() if disease else ""
    species_lower = species.lower().strip() if species else ""
    matrix_lower = matrix.lower().strip() if matrix else ""
    
    for study in summaries:
        study_id = study.get("study_id", "")
        if not study_id:
            continue
        
        title = str(study.get("study_title", "") or "")
        summary_text = str(study.get("study_summary", "") or "")
        study_disease = str(study.get("disease", "") or "")
        study_organism = str(study.get("organism", "") or "")
        study_sample_type = str(study.get("sample_type", "") or "")
        
        # Apply keyword filter (searches title and summary)
        if keyword_lower:
            if keyword_lower not in title.lower() and keyword_lower not in summary_text.lower():
                continue
        
        # Apply disease filter
        if disease_lower:
            if disease_lower not in study_disease.lower() and disease_lower not in title.lower() and disease_lower not in summary_text.lower():
                continue
        
        # Apply species filter
        if species_lower:
            if species_lower not in study_organism.lower():
                continue
        
        # Apply matrix filter
        if matrix_lower:
            if matrix_lower not in study_sample_type.lower():
                continue
        
        # Study matches all filters
        results.append({
            "id": study_id,
            "accession": study_id,  # Dashboard expects this field
            "title": title or f"MW Study {study_id}",
            "description": summary_text[:200] + "..." if len(summary_text) > 200 else summary_text,
            "disease": study_disease,
            "organism": study_organism,
            "omics_type": omics_type or "metabolomics",
        })
        
        if len(results) >= max_results:
            break
    
    logger.info("[REPO][MW] Search returned %d results", len(results))
    return results

