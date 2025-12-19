"""
Unified Domain Models for the Multi-Omics Platform.

This module provides formalized domain models (Pydantic/dataclasses) for all
core entities in the system. These models serve as:
1. Stable abstraction layer
2. Foundation for Postgres schema design
3. Type-safe data structures across the codebase

These models will eventually map to Postgres tables, but for now they provide
a unified interface that can work with both Notion (current) and Postgres (future).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Dict, List, Optional


# ============================================================================
# Enums
# ============================================================================


class OmicsType(str, Enum):
    """Supported omics types."""
    TRANSCRIPTOMICS = "transcriptomics"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    LIPIDOMICS = "lipidomics"


class FeatureType(str, Enum):
    """Feature types across omics."""
    GENE = "gene"
    PROTEIN = "protein"
    METABOLITE = "metabolite"
    LIPID = "lipid"


class SignatureDirection(str, Enum):
    """Direction of change for signature components."""
    UP = "↑"
    DOWN = "↓"
    NEUTRAL = "neutral"
    COMPLEX = "complex"


# ============================================================================
# Core Domain Models
# ============================================================================


@dataclass
class Program:
    """
    A research program/therapeutic area.

    Attributes:
        id: Unique identifier (Notion page ID now, Postgres ID in future)
        name: Program name
        description: Program description
        disease: Disease/indication focus
        created_at: Creation timestamp
        updated_at: Last update timestamp
        notion_page_id: Current Notion page ID (for migration support)
        external_ids: External identifier mappings (DOI, PubMed, etc.)
    """
    id: Optional[str] = None
    name: str = ""
    description: Optional[str] = None
    disease: List[str] = field(default_factory=list)
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    notion_page_id: Optional[str] = None
    external_ids: Dict[str, str] = field(default_factory=dict)


@dataclass
class Experiment:
    """
    An experimental study/assay.

    Attributes:
        id: Unique identifier
        name: Experiment name/title
        type: Experiment type (e.g., "in_vivo", "in_vitro", "patient")
        description: Experiment description
        disease: Disease context
        matrix: Sample matrix (CSF, plasma, tissue, etc.)
        model_systems: Model systems used
        program_ids: Related program IDs
        dataset_ids: Related dataset IDs
        created_at: Creation timestamp
        updated_at: Last update timestamp
        notion_page_id: Current Notion page ID
    """
    id: Optional[str] = None
    name: str = ""
    type: Optional[str] = None
    description: Optional[str] = None
    disease: List[str] = field(default_factory=list)
    matrix: List[str] = field(default_factory=list)
    model_systems: List[str] = field(default_factory=list)
    program_ids: List[str] = field(default_factory=list)
    dataset_ids: List[str] = field(default_factory=list)
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    notion_page_id: Optional[str] = None


@dataclass
class Dataset:
    """
    An experimental dataset (omics data).

    Attributes:
        id: Unique identifier
        name: Dataset name/title
        omics_type: Omics type (transcriptomics, proteomics, etc.)
        description: Dataset description
        file_paths: List of data file paths
        file_urls: List of data file URLs
        organism: Organism(s) studied
        sample_type: Sample types (CSF, plasma, tissue, etc.)
        disease: Disease context
        program_ids: Related program IDs
        experiment_ids: Related experiment IDs
        feature_ids: Related feature IDs (by type)
        signature_ids: Related signature IDs
        signature_match_score: Highest signature match score
        created_at: Creation timestamp
        updated_at: Last update timestamp
        notion_page_id: Current Notion page ID
        external_ids: External identifiers (repository IDs, etc.)
    """
    id: Optional[str] = None
    name: str = ""
    omics_type: OmicsType = OmicsType.TRANSCRIPTOMICS
    description: Optional[str] = None
    file_paths: List[str] = field(default_factory=list)
    file_urls: List[str] = field(default_factory=list)
    organism: List[str] = field(default_factory=list)
    sample_type: List[str] = field(default_factory=list)
    disease: List[str] = field(default_factory=list)
    program_ids: List[str] = field(default_factory=list)
    experiment_ids: List[str] = field(default_factory=list)
    feature_ids: Dict[FeatureType, List[str]] = field(default_factory=dict)
    signature_ids: List[str] = field(default_factory=list)
    signature_match_score: Optional[float] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    notion_page_id: Optional[str] = None
    external_ids: Dict[str, str] = field(default_factory=dict)


@dataclass
class Feature:
    """
    A biological feature (gene, protein, metabolite, lipid).

    Attributes:
        id: Unique identifier
        name: Feature name (e.g., gene symbol, protein name, metabolite name)
        feature_type: Feature type (gene, protein, metabolite, lipid)
        normalized_name: Canonical/normalized name
        aliases: Alternative names/identifiers
        dataset_ids: Datasets containing this feature
        signature_ids: Signatures containing this feature
        created_at: Creation timestamp
        updated_at: Last update timestamp
        notion_page_id: Current Notion page ID
        external_ids: External identifiers (UniProt, KEGG, etc.)
    """
    id: Optional[str] = None
    name: str = ""
    feature_type: FeatureType = FeatureType.GENE
    normalized_name: Optional[str] = None
    aliases: List[str] = field(default_factory=list)
    dataset_ids: List[str] = field(default_factory=list)
    signature_ids: List[str] = field(default_factory=list)
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    notion_page_id: Optional[str] = None
    external_ids: Dict[str, str] = field(default_factory=dict)


@dataclass
class SignatureComponent:
    """
    A single component of a multi-omics signature.

    Attributes:
        feature_id: Feature ID (or name if feature doesn't exist yet)
        feature_type: Feature type
        direction: Direction of change (↑, ↓, neutral, complex)
        weight: Optional weight for this component
    """
    feature_id: Optional[str] = None
    feature_name: str = ""  # For features not yet in DB
    feature_type: FeatureType = FeatureType.GENE
    direction: Optional[SignatureDirection] = None
    weight: Optional[float] = None


@dataclass
class Signature:
    """
    A multi-omics signature definition.

    Attributes:
        id: Unique identifier
        name: Signature name/identifier
        description: Signature description
        components: List of signature components
        modalities: Set of feature types present (auto-computed)
        dataset_ids: Datasets that match this signature
        program_ids: Programs using this signature
        created_at: Creation timestamp
        updated_at: Last update timestamp
        notion_page_id: Current Notion page ID
    """
    id: Optional[str] = None
    name: str = ""
    description: Optional[str] = None
    components: List[SignatureComponent] = field(default_factory=list)
    modalities: List[FeatureType] = field(default_factory=list)
    dataset_ids: List[str] = field(default_factory=list)
    program_ids: List[str] = field(default_factory=list)
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    notion_page_id: Optional[str] = None

    def __post_init__(self):
        """Auto-compute modalities from components if not provided."""
        if not self.modalities:
            self.modalities = list(
                set(
                    comp.feature_type
                    for comp in self.components
                )
            )


# ============================================================================
# Helper Functions
# ============================================================================


def feature_type_from_string(feature_type_str: str) -> FeatureType:
    """Convert string to FeatureType enum."""
    feature_type_str = feature_type_str.lower().strip()
    for ft in FeatureType:
        if ft.value == feature_type_str:
            return ft
    # Default to GENE if unknown
    return FeatureType.GENE


def omics_type_from_string(omics_type_str: str) -> OmicsType:
    """Convert string to OmicsType enum."""
    omics_type_str = omics_type_str.lower().strip()
    for ot in OmicsType:
        if ot.value == omics_type_str:
            return ot
    # Default to TRANSCRIPTOMICS if unknown
    return OmicsType.TRANSCRIPTOMICS

