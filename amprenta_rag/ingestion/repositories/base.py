"""
Base repository interface for public omics data repositories.

Defines the common interface that all repository implementations must follow.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional

from amprenta_rag.models.repository import DataFile, StudyMetadata


class RepositoryInterface(ABC):
    """
    Base interface for all omics data repositories.
    
    All repository implementations (GEO, PRIDE, MetaboLights, MW) must
    implement these methods to provide a unified interface for discovery
    and data harvesting.
    """
    
    @abstractmethod
    def get_repository_name(self) -> str:
        """
        Get the name of this repository.
        
        Returns:
            Repository name (e.g., "GEO", "PRIDE", "MetaboLights", "MW")
        """
    
    @abstractmethod
    def get_omics_type(self) -> str:
        """
        Get the primary omics type this repository serves.
        
        Returns:
            Omics type (e.g., "transcriptomics", "proteomics", "metabolomics", "lipidomics")
        """
    
    @abstractmethod
    def search_studies(
        self,
        keywords: List[str],
        filters: Optional[Dict[str, Any]] = None,
        max_results: int = 100,
    ) -> List[str]:
        """
        Search for studies matching keywords and filters.
        
        Args:
            keywords: List of search keywords
            filters: Optional filters (e.g., {"disease": "ALS", "organism": "Homo sapiens"})
            max_results: Maximum number of results to return
            
        Returns:
            List of study IDs matching the search criteria
        """
    
    @abstractmethod
    def fetch_study_metadata(self, study_id: str) -> Optional[StudyMetadata]:
        """
        Fetch metadata for a specific study.
        
        Args:
            study_id: Repository-specific study identifier
            
        Returns:
            StudyMetadata object, or None if study not found
        """
    
    @abstractmethod
    def fetch_study_data_files(
        self,
        study_id: str,
        file_types: Optional[List[str]] = None,
    ) -> List[DataFile]:
        """
        Fetch list of available data files for a study.
        
        Args:
            study_id: Repository-specific study identifier
            file_types: Optional filter by file types (e.g., ["CSV", "TSV"])
            
        Returns:
            List of DataFile objects
        """
    
    def download_data_file(
        self,
        study_id: str,
        file_id: str,
        output_path: str,
    ) -> bool:
        """
        Download a specific data file from a study.
        
        Default implementation that can be overridden by repositories
        that need custom download logic.
        
        Args:
            study_id: Repository-specific study identifier
            file_id: File identifier
            output_path: Local path to save the file
            
        Returns:
            True if download successful, False otherwise
        """
        # Default implementation - repositories can override
        metadata = self.fetch_study_metadata(study_id)
        if not metadata:
            return False
        
        # Find the file
        data_file = None
        for f in metadata.data_files:
            if f.file_id == file_id:
                data_file = f
                break
        
        if not data_file or not data_file.download_url:
            return False
        
        # Download using requests
        import requests
        from pathlib import Path
        
        try:
            resp = requests.get(data_file.download_url, timeout=300, stream=True)
            resp.raise_for_status()
            
            output_path_obj = Path(output_path)
            output_path_obj.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            return True
        except Exception as e:
            from amprenta_rag.logging_utils import get_logger
            logger = get_logger(__name__)
            logger.error(
                "[REPO][%s] Error downloading file %s: %r",
                self.get_repository_name(),
                file_id,
                e,
            )
            return False

