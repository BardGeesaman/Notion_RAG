"""
Domain models for public repository ingestion.

Defines data structures for study metadata and data files from various
omics repositories (GEO, PRIDE, MetaboLights, Metabolomics Workbench).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional


@dataclass
class DataFile:
    """
    Represents a data file from a repository study.

    Attributes:
        file_id: Repository-specific file identifier
        filename: Original filename
        file_type: File type (e.g., "CSV", "TSV", "TXT", "XLSX")
        download_url: URL to download the file
        size_bytes: File size in bytes (if available)
        description: Optional file description
    """
    file_id: str
    filename: str
    file_type: str
    download_url: Optional[str] = None
    size_bytes: Optional[int] = None
    description: Optional[str] = None


@dataclass
class StudyMetadata:
    """
    Unified metadata model for studies from any repository.

    Attributes:
        study_id: Repository-specific study identifier (e.g., "GSE12345", "PXD012345")
        repository: Repository name ("GEO", "PRIDE", "MetaboLights", "MW")
        title: Study title
        summary: Study summary/description
        doi: DOI if available
        pubmed_id: PubMed ID if available
        disease: List of disease terms
        organism: List of organism names
        sample_type: List of sample types (CSF, plasma, tissue, etc.)
        omics_type: Omics type ("transcriptomics", "proteomics", "metabolomics", "lipidomics")
        data_files: List of available data files
        publication_date: Publication date if available
        authors: List of author names
        keywords: List of keywords/tags
        platform: Platform/technology used (e.g., "RNA-seq", "LC-MS")
        num_samples: Number of samples (if available)
        raw_metadata: Raw metadata dictionary from repository (for reference)
    """
    study_id: str
    repository: str
    title: str
    summary: str
    omics_type: str

    # Optional fields
    doi: Optional[str] = None
    pubmed_id: Optional[str] = None
    disease: List[str] = field(default_factory=list)
    organism: List[str] = field(default_factory=list)
    sample_type: List[str] = field(default_factory=list)
    data_files: List[DataFile] = field(default_factory=list)
    publication_date: Optional[datetime] = None
    authors: List[str] = field(default_factory=list)
    keywords: List[str] = field(default_factory=list)
    platform: Optional[str] = None
    num_samples: Optional[int] = None
    raw_metadata: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "study_id": self.study_id,
            "repository": self.repository,
            "title": self.title,
            "summary": self.summary,
            "omics_type": self.omics_type,
            "doi": self.doi,
            "pubmed_id": self.pubmed_id,
            "disease": self.disease,
            "organism": self.organism,
            "sample_type": self.sample_type,
            "data_files": [
                {
                    "file_id": f.file_id,
                    "filename": f.filename,
                    "file_type": f.file_type,
                    "download_url": f.download_url,
                    "size_bytes": f.size_bytes,
                    "description": f.description,
                }
                for f in self.data_files
            ],
            "publication_date": self.publication_date.isoformat() if self.publication_date else None,
            "authors": self.authors,
            "keywords": self.keywords,
            "platform": self.platform,
            "num_samples": self.num_samples,
        }

