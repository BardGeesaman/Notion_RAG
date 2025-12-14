"""
Pydantic schemas for API request/response models.

These schemas define the structure of API requests and responses,
mapping between domain models and API representations.
"""

from __future__ import annotations

from datetime import datetime
from typing import List, Optional
from uuid import UUID

from pydantic import BaseModel, Field, field_validator

from amprenta_rag.models.domain import FeatureType, OmicsType, SignatureDirection


# ============================================================================
# Common schemas
# ============================================================================


class BaseSchema(BaseModel):
    """Base schema with common configuration."""
    
    class Config:
        from_attributes = True
        json_encoders = {
            UUID: str,
            datetime: lambda v: v.isoformat(),
        }


class AnnotationCreate(BaseSchema):
    """Schema for creating an annotation/note on an entity."""

    text: str
    annotation_type: Optional[str] = None


# ============================================================================
# Program schemas
# ============================================================================


class ProgramBase(BaseSchema):
    """Base program schema."""
    name: str
    description: Optional[str] = None
    disease: Optional[List[str]] = Field(default_factory=list)


class ProgramCreate(ProgramBase):
    """Schema for creating a program."""


class ProgramUpdate(BaseSchema):
    """Schema for updating a program."""
    name: Optional[str] = None
    description: Optional[str] = None
    disease: Optional[List[str]] = None


class Program(ProgramBase):
    """Program response schema."""
    id: UUID
    created_at: datetime
    updated_at: datetime
    notion_page_id: Optional[str] = None
    external_ids: Optional[dict] = None
    
    # Override disease to handle None from database
    disease: Optional[List[str]] = None


# ============================================================================
# Experiment schemas
# ============================================================================


class ExperimentBase(BaseSchema):
    """Base experiment schema."""
    name: str
    type: Optional[str] = None
    description: Optional[str] = None
    disease: Optional[List[str]] = None
    matrix: Optional[List[str]] = None
    model_systems: Optional[List[str]] = None

    @field_validator("disease", "matrix", "model_systems", mode="before")
    @classmethod
    def ensure_list(cls, v):
        """Convert None to empty list for nullable array fields."""
        return v if v is not None else []


class ExperimentCreate(ExperimentBase):
    """Schema for creating an experiment."""
    program_ids: List[UUID] = Field(default_factory=list)


class ExperimentUpdate(BaseSchema):
    """Schema for updating an experiment."""
    name: Optional[str] = None
    type: Optional[str] = None
    description: Optional[str] = None
    disease: Optional[List[str]] = None
    matrix: Optional[List[str]] = None
    model_systems: Optional[List[str]] = None
    program_ids: Optional[List[UUID]] = None


class Experiment(ExperimentBase):
    """Experiment response schema."""
    id: UUID
    program_ids: List[UUID] = Field(default_factory=list)
    dataset_ids: List[UUID] = Field(default_factory=list)
    created_at: datetime
    updated_at: datetime
    notion_page_id: Optional[str] = None


# ============================================================================
# Dataset schemas
# ============================================================================


class DatasetBase(BaseSchema):
    """Base dataset schema."""
    name: str
    omics_type: OmicsType
    description: Optional[str] = None
    file_paths: List[str] = Field(default_factory=list)
    file_urls: List[str] = Field(default_factory=list)
    organism: List[str] = Field(default_factory=list)
    sample_type: List[str] = Field(default_factory=list)
    disease: List[str] = Field(default_factory=list)


class DatasetCreate(DatasetBase):
    """Schema for creating a dataset."""
    program_ids: List[UUID] = Field(default_factory=list)
    experiment_ids: List[UUID] = Field(default_factory=list)


class DatasetUpdate(BaseSchema):
    """Schema for updating a dataset."""
    name: Optional[str] = None
    omics_type: Optional[OmicsType] = None
    description: Optional[str] = None
    file_paths: Optional[List[str]] = None
    file_urls: Optional[List[str]] = None
    organism: Optional[List[str]] = None
    sample_type: Optional[List[str]] = None
    disease: Optional[List[str]] = None
    program_ids: Optional[List[UUID]] = None
    experiment_ids: Optional[List[UUID]] = None
    signature_match_score: Optional[float] = None


class Dataset(DatasetBase):
    """Dataset response schema."""
    id: UUID
    program_ids: List[UUID] = Field(default_factory=list)
    experiment_ids: List[UUID] = Field(default_factory=list)
    feature_ids: dict[FeatureType, List[UUID]] = Field(default_factory=dict)
    signature_ids: List[UUID] = Field(default_factory=list)
    signature_match_score: Optional[float] = None
    created_at: datetime
    updated_at: datetime
    notion_page_id: Optional[str] = None
    external_ids: Optional[dict] = None


# ============================================================================
# Feature schemas
# ============================================================================


class FeatureBase(BaseSchema):
    """Base feature schema."""
    name: str
    feature_type: FeatureType
    normalized_name: Optional[str] = None
    aliases: List[str] = Field(default_factory=list)


class FeatureCreate(FeatureBase):
    """Schema for creating a feature."""
    external_ids: Optional[dict] = None


class FeatureUpdate(BaseSchema):
    """Schema for updating a feature."""
    name: Optional[str] = None
    normalized_name: Optional[str] = None
    aliases: Optional[List[str]] = None
    external_ids: Optional[dict] = None


class Feature(FeatureBase):
    """Feature response schema."""
    id: UUID
    dataset_ids: List[UUID] = Field(default_factory=list)
    signature_ids: List[UUID] = Field(default_factory=list)
    created_at: datetime
    updated_at: datetime
    notion_page_id: Optional[str] = None
    external_ids: Optional[dict] = None


# ============================================================================
# Chemistry / Screening schemas
# ============================================================================


class CompoundResponse(BaseSchema):
    compound_id: str
    smiles: str
    inchi_key: Optional[str] = None
    canonical_smiles: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    hbd_count: Optional[int] = None
    hba_count: Optional[int] = None
    rotatable_bonds: Optional[int] = None


class ProgramLinkResponse(BaseSchema):
    program_id: str
    status: Optional[str] = None
    notes: Optional[str] = None
    created_at: Optional[str] = None


class CampaignResponse(BaseSchema):
    campaign_id: str
    campaign_name: str
    description: Optional[str] = None
    assay_type: Optional[str] = None
    target: Optional[str] = None
    library_id: Optional[str] = None
    total_wells: Optional[int] = None
    hit_count: Optional[int] = None
    run_date: Optional[str] = None


class HTSHitResponse(BaseSchema):
    result_id: str
    compound_id: str
    well_position: Optional[str] = None
    raw_value: Optional[float] = None
    normalized_value: Optional[float] = None
    z_score: Optional[float] = None
    hit_flag: int
    hit_category: Optional[str] = None


# ============================================================================
# Signature schemas
# ============================================================================


class SignatureComponentBase(BaseSchema):
    """Base signature component schema."""
    feature_name: str
    feature_type: FeatureType
    direction: Optional[SignatureDirection] = None
    weight: Optional[float] = 1.0


class SignatureComponentCreate(SignatureComponentBase):
    """Schema for creating a signature component."""
    feature_id: Optional[UUID] = None


class SignatureComponent(SignatureComponentBase):
    """Signature component response schema."""
    id: UUID
    signature_id: UUID
    feature_id: Optional[UUID] = None


class SignatureBase(BaseSchema):
    """Base signature schema."""
    name: str
    description: Optional[str] = None


class SignatureCreate(SignatureBase):
    """Schema for creating a signature."""
    components: List[SignatureComponentCreate] = Field(default_factory=list)
    program_ids: List[UUID] = Field(default_factory=list)


class SignatureUpdate(BaseSchema):
    """Schema for updating a signature."""
    name: Optional[str] = None
    description: Optional[str] = None
    components: Optional[List[SignatureComponentCreate]] = None
    program_ids: Optional[List[UUID]] = None


class Signature(SignatureBase):
    """Signature response schema."""
    id: UUID
    components: List[SignatureComponent] = Field(default_factory=list)
    modalities: List[FeatureType] = Field(default_factory=list)
    dataset_ids: List[UUID] = Field(default_factory=list)
    program_ids: List[UUID] = Field(default_factory=list)
    created_at: datetime
    updated_at: datetime
    notion_page_id: Optional[str] = None

