"""
Pydantic schemas for API request/response models.

These schemas define the structure of API requests and responses,
mapping between domain models and API representations.
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import Any, Dict, List, Literal, Optional, Tuple
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field, field_validator

from amprenta_rag.models.domain import FeatureType, OmicsType, SignatureDirection

# Generative chemistry schemas (added inline to avoid circular imports)


# ============================================================================
# Common schemas
# ============================================================================


class BaseSchema(BaseModel):
    """Base schema with common configuration."""

    model_config = ConfigDict(
        from_attributes=True,
        json_encoders={
            UUID: str,
            datetime: lambda v: v.isoformat(),
        }
    )


class StrictBaseSchema(BaseModel):
    """Strict schema for security-sensitive endpoints.
    
    Use this for endpoints that handle:
    - Authentication/passwords
    - URLs/external resources
    - Chemical structures (SMILES)
    - File paths
    
    Strict mode prevents type coercion (e.g., "1" → 1) which can
    bypass validation in edge cases.
    """
    
    model_config = ConfigDict(
        from_attributes=True,
        strict=True,
        json_encoders={
            UUID: str,
            datetime: lambda v: v.isoformat(),
        }
    )


def validate_safe_string(v: str, max_length: int = 10000) -> str:
    """Validate string is safe (no control chars, reasonable length)."""
    if not v:
        return v
    # Remove null bytes and control characters (except newline, tab)
    v = re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]', '', v)
    if len(v) > max_length:
        raise ValueError(f"String exceeds maximum length of {max_length}")
    return v


def validate_smiles(v: str) -> str:
    """Validate SMILES string is safe and well-formed."""
    if not v or not v.strip():
        raise ValueError("SMILES cannot be empty")
    v = v.strip()
    # Basic SMILES character validation (conservative)
    if not re.match(r'^[A-Za-z0-9@+\-\[\]()=#%.$\\/:*]+$', v):
        raise ValueError("Invalid characters in SMILES")
    if len(v) > 5000:
        raise ValueError("SMILES too long (max 5000 chars)")
    return v


def validate_url(v: str) -> str:
    """Validate URL is safe (no javascript:, data:, etc.)."""
    if not v:
        return v
    v = v.strip()
    # Block dangerous URL schemes
    dangerous_schemes = ['javascript:', 'data:', 'vbscript:', 'file:']
    v_lower = v.lower()
    for scheme in dangerous_schemes:
        if v_lower.startswith(scheme):
            raise ValueError(f"Dangerous URL scheme: {scheme}")
    return v


# Retrosynthesis schemas
class SynthesisStepSchema(BaseModel):
    reactants: List[str]
    product: str
    reaction_type: str
    conditions: str
    confidence: float = Field(ge=0.0, le=1.0)
    
    model_config = ConfigDict(from_attributes=True)


class SynthesisRouteSchema(BaseModel):
    id: str
    steps: List[SynthesisStepSchema]
    total_steps: int
    confidence: float = Field(ge=0.0, le=1.0)
    
    model_config = ConfigDict(from_attributes=True)


class SynthesisTreeSchema(BaseModel):
    target: str
    routes: List[SynthesisRouteSchema]
    analysis_time_ms: int
    num_alternatives: int
    
    model_config = ConfigDict(from_attributes=True)


class RetrosynthesisRequest(BaseModel):
    smiles: str
    max_depth: int = Field(default=5, ge=1, le=10)


class RetrosynthesisResponse(BaseModel):
    analysis_id: str
    tree: SynthesisTreeSchema


class RouteScoreSchema(BaseModel):
    total_score: float = Field(ge=0.0, le=100.0)
    step_count: int
    complexity_score: float
    availability_score: float
    estimated_cost: float
    
    model_config = ConfigDict(from_attributes=True)


class BuildingBlockResultSchema(BaseModel):
    smiles: str
    available: bool
    vendors: List[dict]
    
    model_config = ConfigDict(from_attributes=True)


# Target Management Schemas
class TargetCreate(BaseSchema):
    """Schema for creating a new target."""
    name: str = Field(..., min_length=1, max_length=200)
    gene_symbol: Optional[str] = Field(None, max_length=50)
    uniprot_id: Optional[str] = Field(None, max_length=20)
    target_class: Optional[str] = Field(None, max_length=100)
    target_family: Optional[str] = Field(None, max_length=100)
    description: Optional[str] = None


class TargetUpdate(BaseSchema):
    """Schema for updating target fields."""
    name: Optional[str] = Field(None, min_length=1, max_length=200)
    description: Optional[str] = None
    gene_symbol: Optional[str] = Field(None, max_length=50)
    target_class: Optional[str] = Field(None, max_length=100)
    target_family: Optional[str] = Field(None, max_length=100)
    validation_status: Optional[str] = Field(None, max_length=50)
    lifecycle_status: Optional[str] = Field(None, max_length=50)
    druggability_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    pocket_count: Optional[int] = Field(None, ge=0)
    is_synthetic: Optional[bool] = None


class TargetResponse(BaseSchema):
    """Schema for target response data."""
    id: UUID
    name: str
    description: Optional[str] = None
    gene_symbol: Optional[str] = None
    uniprot_id: Optional[str] = None
    target_class: Optional[str] = None
    target_family: Optional[str] = None
    druggability_score: Optional[float] = None
    druggability_source: Optional[str] = None
    pocket_count: Optional[int] = None
    validation_status: Optional[str] = None
    validation_evidence: Optional[dict] = None
    lifecycle_status: str = "active"
    is_synthetic: bool = False
    chembl_id: Optional[str] = None
    ensembl_id: Optional[str] = None
    created_at: datetime
    updated_at: datetime
    
    model_config = ConfigDict(from_attributes=True)


class CompoundActivityResponse(BaseSchema):
    """Schema for compound activity data."""
    compound_id: UUID
    compound_name: Optional[str] = None
    smiles: str
    activity_type: Optional[str] = None
    activity_value: Optional[float] = None
    activity_units: Optional[str] = None
    assay_id: Optional[UUID] = None
    
    model_config = ConfigDict(from_attributes=True)


class DruggabilityResponse(BaseSchema):
    """Schema for druggability calculation results."""
    target_id: UUID
    score: float = Field(ge=0.0, le=1.0)
    factors: dict
    
    model_config = ConfigDict(from_attributes=True)


class TargetAssayResponse(BaseSchema):
    """Schema for target-linked assay data."""
    id: UUID
    name: str
    description: Optional[str] = None
    experiment_type: Optional[str] = None
    created_at: datetime
    
    model_config = ConfigDict(from_attributes=True)


class AnnotationCreate(BaseSchema):
    """Schema for creating an annotation/note on an entity."""

    text: str
    annotation_type: Optional[str] = None


# Import PriorConfig from analysis layer to avoid circular imports


# ============================================================================
# LLM-based Planning schemas
# ============================================================================


class PlanRequest(BaseModel):
    """Request for LLM-based experiment planning."""
    
    goal: str = Field(..., description="Research goal or hypothesis")
    context: Optional[ScoringContextRequest] = Field(None, description="Research context")
    constraints: Optional[List[str]] = Field(None, description="Budget, time, or resource constraints")
    

class CritiqueRequest(BaseModel):
    """Request for LLM-based plan critique."""
    
    plan: str = Field(..., description="Experimental plan to critique")
    criteria: Optional[List[str]] = Field(None, description="Specific criteria to evaluate")
    context: Optional[ScoringContextRequest] = Field(None, description="Research context")


class RefineRequest(BaseModel):
    """Request for LLM-based plan refinement."""
    
    original_plan: str = Field(..., description="Original experimental plan")
    critique: str = Field(..., description="Critique feedback")
    additional_requirements: Optional[str] = Field(None, description="Additional requirements")


class ExecuteRequest(BaseModel):
    """Request for LLM-based execution guidance."""
    
    plan: str = Field(..., description="Experimental plan to execute")
    current_step: Optional[int] = Field(None, description="Current execution step")
    issues: Optional[List[str]] = Field(None, description="Issues encountered during execution")


class PlanStep(BaseModel):
    """Individual step in an experimental plan."""
    
    step_number: int
    description: str
    duration: Optional[str] = None
    resources: Optional[List[str]] = None
    dependencies: Optional[List[int]] = None


class PlanResult(BaseModel):
    """LLM-based planning result."""
    
    plan_id: str
    title: str
    objective: str
    steps: List[PlanStep]
    estimated_duration: Optional[str] = None
    estimated_cost: Optional[str] = None
    risks: Optional[List[str]] = None
    processing_time_seconds: float
    cached: bool = False


class CritiqueResult(BaseModel):
    """LLM-based critique result."""
    
    overall_score: float = Field(..., ge=0, le=10, description="Overall plan quality score")
    strengths: List[str]
    weaknesses: List[str]
    recommendations: List[str]
    feasibility_score: Optional[float] = Field(None, ge=0, le=10)
    processing_time_seconds: float
    cached: bool = False


class RefinementResult(BaseModel):
    """LLM-based refinement result."""
    
    refined_plan: PlanResult
    changes_made: List[str]
    rationale: str
    processing_time_seconds: float
    cached: bool = False


class ExecutionGuidance(BaseModel):
    """LLM-based execution guidance."""
    
    current_step: PlanStep
    next_steps: List[PlanStep]
    recommendations: List[str]
    potential_issues: Optional[List[str]] = None
    troubleshooting: Optional[List[str]] = None
    processing_time_seconds: float
    cached: bool = False


# ============================================================================
# LLM-based Ranking schemas
# ============================================================================


class RankRequest(BaseModel):
    """Request for LLM-based ranking of items."""
    
    items: List[ScoringItemRequest] = Field(..., min_length=1, description="Items to rank")
    criteria: Optional[str] = Field(None, description="Ranking criteria")
    context: Optional[ScoringContextRequest] = Field(None, description="Ranking context")


class RerankRequest(BaseModel):
    """Request for LLM-based re-ranking of items."""
    
    items: List[ScoringItemRequest] = Field(..., min_length=1, description="Items to re-rank")
    previous_ranking: List[str] = Field(..., description="Previous ranking order (item IDs)")
    criteria: Optional[str] = Field(None, description="Re-ranking criteria")
    context: Optional[ScoringContextRequest] = Field(None, description="Re-ranking context")


class RankingItem(BaseModel):
    """Ranked item result."""
    
    item_id: str
    rank: int
    score: Optional[float] = None
    explanation: Optional[str] = None


class RankingResult(BaseModel):
    """LLM-based ranking result."""
    
    ranked_items: List[RankingItem]
    criteria_used: Optional[str] = None
    explanation: Optional[str] = None
    processing_time_seconds: float
    cached: bool = False


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
    def ensure_list(cls, v: Optional[List[str]]) -> List[str]:
        """Convert None to empty list for nullable array fields."""
        return v if v is not None else []


class ExperimentCreate(ExperimentBase):
    """Schema for creating an experiment."""
    program_ids: List[UUID] = Field(default_factory=list)


class ExperimentUpdate(BaseSchema):
    """Schema for updating an experiment."""
    version: int = Field(..., ge=1, description="Current entity version for optimistic locking")
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
    version: int = 1


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
    version: int = Field(..., ge=1, description="Current entity version for optimistic locking")
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
    file_paths: List[str] = Field(default_factory=list)
    file_urls: List[str] = Field(default_factory=list)
    organism: List[str] = Field(default_factory=list)
    version: int = 1
    dvc_version: Optional[str] = None
    dvc_metadata: Optional[dict] = None
    dvc_pushed: bool = False

    @field_validator("file_paths", "file_urls", "organism", "sample_type", "disease", mode="before")
    @classmethod
    def none_to_empty_list(cls, v: Optional[List[str]]) -> List[str]:
        return v if v is not None else []


# ============================================================================
# Catalog / Repository schemas
# ============================================================================


class RepositorySummary(BaseSchema):
    """Summary of datasets available per external repository."""

    name: str
    dataset_count: int
    last_sync_date: Optional[datetime] = None
    health_status: str


class CatalogDataset(BaseSchema):
    """Dataset entry for the external catalog."""

    id: UUID
    accession: Optional[str] = None
    title: Optional[str] = None
    source: Optional[str] = None
    created_at: datetime
    status: Optional[str] = None
    feature_count: int = 0


class CatalogSummaryResponse(BaseSchema):
    """Response wrapper for catalog summary."""

    repositories: List[RepositorySummary]


# ============================================================================
# Dataset Finder schemas
# ============================================================================


class DatasetFinderRequest(BaseModel):
    """Request schema for AI dataset finder."""
    
    query: str = Field(..., min_length=1, max_length=500, description="Natural language search query")
    repositories: Optional[List[str]] = Field(default=None, description="Repositories to search (geo, arrayexpress, metabolomics_workbench)")
    max_results: int = Field(default=50, ge=1, le=200, description="Maximum number of results to return")


class DatasetResultSchema(BaseModel):
    """Schema for a single dataset search result."""
    
    accession: str
    title: str
    description: str
    source: str
    species: Optional[str] = None
    tissue: Optional[str] = None
    disease: Optional[str] = None
    assay_type: Optional[str] = None
    sample_count: Optional[int] = None
    url: Optional[str] = None
    score: float = 0.0


class DatasetFinderResponse(BaseModel):
    """Response schema for AI dataset finder."""
    
    query: str
    extracted_terms: Dict[str, List[str]]
    results: List[DatasetResultSchema]
    total_found: int
    sources_searched: List[str]
    sources_failed: List[str]


class MetadataEnrichmentResponse(BaseModel):
    """Response schema for metadata enrichment."""
    
    dataset_id: UUID
    success: bool
    enriched_fields: List[str]
    extracted_metadata: Dict[str, Any]
    error_message: Optional[str] = None
    processing_time_seconds: Optional[float] = None


class EnrichmentStatusResponse(BaseModel):
    """Response schema for enrichment status."""
    
    enriched: bool
    enrichment_timestamp: Optional[float] = None
    enrichment_version: Optional[str] = None
    available_fields: List[str] = []
    error: Optional[str] = None


# ============================================================================
# Relevance & Novelty Scoring schemas
# ============================================================================


class ScoringItemRequest(BaseModel):
    """Schema for an item to be scored."""
    
    id: str = Field(..., description="Unique identifier for the item")
    title: Optional[str] = None
    description: Optional[str] = None
    species: Optional[str] = None
    assay_type: Optional[str] = None
    sample_count: Optional[int] = None
    # Additional fields can be added as needed


class ScoringContextRequest(BaseModel):
    """Schema for research context used in scoring."""
    
    diseases: Optional[List[str]] = Field(default=[], description="Diseases of interest")
    targets: Optional[List[str]] = Field(default=[], description="Target molecules of interest")
    species: Optional[List[str]] = Field(default=[], description="Preferred species")
    assay_types: Optional[List[str]] = Field(default=[], description="Required assay types")
    min_sample_size: Optional[int] = Field(default=None, description="Minimum sample size")


class RelevanceScoreRequest(BaseModel):
    """Request schema for relevance scoring."""
    
    item: ScoringItemRequest
    context: ScoringContextRequest
    criteria: Optional[Dict[str, float]] = Field(default=None, description="Custom criteria weights")


class NoveltyScoreRequest(BaseModel):
    """Request schema for novelty scoring."""
    
    item: ScoringItemRequest
    existing_items: List[ScoringItemRequest] = Field(..., description="Existing items to compare against")


class BatchScoreRequest(BaseModel):
    """Request schema for batch scoring."""
    
    items: List[ScoringItemRequest] = Field(..., min_length=1, max_length=50, description="Items to score")
    context: ScoringContextRequest
    score_relevance: bool = Field(default=True, description="Whether to compute relevance scores")
    score_novelty: bool = Field(default=True, description="Whether to compute novelty scores")


class RelevanceScoreResponse(BaseModel):
    """Response schema for relevance scoring."""
    
    item_id: str
    overall_score: float = Field(..., ge=0.0, le=1.0, description="Overall relevance score (0-1)")
    disease_match: float = Field(..., ge=0.0, le=1.0, description="Disease match score (0-1)")
    target_overlap: float = Field(..., ge=0.0, le=1.0, description="Target overlap score (0-1)")
    data_quality: float = Field(..., ge=0.0, le=1.0, description="Data quality score (0-1)")
    explanation: str
    processing_time_seconds: float
    cached: bool = Field(default=False, description="Whether result was retrieved from cache")


class NoveltyScoreResponse(BaseModel):
    """Response schema for novelty scoring."""
    
    item_id: str
    novelty_score: float = Field(..., ge=0.0, le=1.0, description="Novelty score (0-1, higher = more novel)")
    max_similarity: float = Field(..., ge=0.0, le=1.0, description="Maximum similarity to existing items")
    most_similar_item_id: Optional[str] = Field(default=None, description="ID of most similar existing item")
    explanation: str
    processing_time_seconds: float
    cached: bool = Field(default=False, description="Whether result was retrieved from cache")


class ScoredItemResponse(BaseModel):
    """Response schema for a scored item."""
    
    item_id: str
    relevance_score: Optional[RelevanceScoreResponse] = None
    novelty_score: Optional[NoveltyScoreResponse] = None


class BatchScoreResponse(BaseModel):
    """Response schema for batch scoring."""
    
    items: List[ScoredItemResponse]
    total_items: int
    processing_time_seconds: float


# ============================================================================
# Assay Predictor schemas
# ============================================================================


class AssayPredictorTrainRequest(BaseModel):
    """Request schema for training assay predictor."""
    
    program_id: UUID = Field(..., description="UUID of the program to train model for")
    assay_type: str = Field(..., description="Type of assay (biochemical, hts, screening)")
    features: Optional[List[str]] = Field(default=None, description="Optional specific features to use")
    min_actives: int = Field(default=50, ge=10, le=1000, description="Minimum active compounds required")
    min_inactives: int = Field(default=50, ge=10, le=1000, description="Minimum inactive compounds required")


class TrainingDataStatsResponse(BaseModel):
    """Response schema for training data statistics."""
    
    total_compounds: int
    active_compounds: int
    inactive_compounds: int
    activity_rate: float = Field(..., ge=0.0, le=1.0, description="Fraction of active compounds")
    feature_count: int
    data_quality_score: float = Field(..., ge=0.0, le=1.0, description="Overall data quality score")


class AssayPredictorTrainResponse(BaseModel):
    """Response schema for assay predictor training."""
    
    model_id: UUID
    program_id: UUID
    assay_type: str
    model_performance: Dict[str, float]
    training_stats: TrainingDataStatsResponse
    feature_names: List[str]
    training_time_seconds: float
    success: bool
    error_message: Optional[str] = None


class AssayPredictionRequest(BaseModel):
    """Request schema for assay prediction."""
    
    compound_smiles: List[str] = Field(..., min_length=1, max_length=1000, description="SMILES strings to predict")


class PredictionResultResponse(BaseModel):
    """Response schema for a single prediction result."""
    
    compound_smiles: str
    prediction: str = Field(..., description="Predicted outcome: active or inactive")
    probability_active: float = Field(..., ge=0.0, le=1.0, description="Probability of being active")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Prediction confidence")
    feature_vector: Optional[List[float]] = Field(default=None, description="Optional feature vector")


class AssayPredictionResponse(BaseModel):
    """Response schema for assay prediction."""
    
    model_id: UUID
    predictions: List[PredictionResultResponse]
    total_predictions: int
    processing_time_seconds: float


class AssayModelResponse(BaseModel):
    """Response schema for assay model metadata."""
    
    model_id: UUID
    name: str
    version: str
    program_id: Optional[str] = None
    assay_type: Optional[str] = None
    created_at: datetime
    performance: Dict[str, float] = {}
    training_stats: Dict[str, Any] = {}


class AssayModelsListResponse(BaseModel):
    """Response schema for listing assay models."""
    
    models: List[AssayModelResponse]
    total_models: int


# ============================================================================
# Active Learning schemas
# ============================================================================


class ScreenedCompoundRequest(BaseModel):
    """Schema for a screened compound."""
    
    compound_id: UUID
    smiles: str
    activity: bool = Field(..., description="True if compound is active")


class CandidateCompoundRequest(BaseModel):
    """Schema for a candidate compound."""
    
    compound_id: UUID
    smiles: str


class ActiveLearningRequest(BaseModel):
    """Request schema for active learning suggestions."""
    
    screened: List[ScreenedCompoundRequest] = Field(..., description="Already screened compounds with results")
    candidates: List[CandidateCompoundRequest] = Field(..., min_length=1, description="Candidate compounds to choose from")
    strategy: str = Field(..., description="Active learning strategy: uncertainty or diversity")
    batch_size: int = Field(default=10, ge=1, le=100, description="Number of compounds to suggest")
    model_id: Optional[UUID] = Field(default=None, description="Model ID for uncertainty-based strategies")


class ProgramActiveLearningRequest(BaseModel):
    """Request schema for program-based active learning."""
    
    program_id: UUID
    strategy: str = Field(default="uncertainty", description="Active learning strategy: uncertainty or diversity")
    batch_size: int = Field(default=10, ge=1, le=100, description="Number of compounds to suggest")
    model_id: Optional[UUID] = Field(default=None, description="Model ID for uncertainty-based strategies")


class SuggestionResultResponse(BaseModel):
    """Response schema for a compound suggestion."""
    
    compound_id: UUID
    smiles: str
    acquisition_score: float = Field(..., ge=0.0, description="Acquisition score (higher = better)")
    strategy_used: str
    rank: int = Field(..., ge=1, description="Rank in suggestion list")
    explanation: str


class ActiveLearningResponse(BaseModel):
    """Response schema for active learning suggestions."""
    
    suggestions: List[SuggestionResultResponse]
    strategy_used: str
    total_candidates: int
    total_screened: int
    batch_size: int
    processing_time_seconds: float


class AcquisitionScoreResponse(BaseModel):
    """Response schema for acquisition scores."""
    
    compound_id: UUID
    score: float = Field(..., ge=0.0, description="Acquisition score")
    explanation: str


class StrategyPerformanceResponse(BaseModel):
    """Response schema for strategy performance evaluation."""
    
    hit_rate: float = Field(..., ge=0.0, le=1.0, description="Fraction of suggested compounds that are active")
    enrichment: float = Field(..., ge=0.0, description="Enrichment over random selection")
    coverage: float = Field(..., ge=0.0, le=1.0, description="Fraction of suggestions with known activity")
    total_suggested: int
    total_hits: int


# ============================================================================
# ADMET prediction schemas
# ============================================================================


class ADMETPredictRequest(BaseSchema):
    smiles: List[str] = Field(..., max_length=100, description="SMILES list (max 100)")
    endpoints: Optional[List[str]] = None
    include_uncertainty: bool = True


class ADMETEndpointPrediction(BaseSchema):
    mean: float
    std: Optional[float] = None
    ci_low: Optional[float] = None
    ci_high: Optional[float] = None
    in_domain: Optional[bool] = None
    similarity: Optional[float] = None
    calibrated: bool = False


class ADMETCompoundPrediction(BaseSchema):
    smiles: str
    predictions: Dict[str, ADMETEndpointPrediction]
    error: Optional[str] = None


class ADMETPredictResponse(BaseSchema):
    results: List[ADMETCompoundPrediction]
    model_info: Dict[str, Any]


class ADMETExplainRequest(BaseSchema):
    smiles: str
    endpoint: str = Field("herg", description="ADMET endpoint (herg|logs|logp)")
    top_k: int = Field(10, ge=1, le=100, description="Number of top SHAP features")


class SHAPFeature(BaseSchema):
    name: str
    value: float
    rank: int


class ADMETExplainResponse(BaseSchema):
    smiles: str
    endpoint: str
    prediction: Dict[str, Any]
    shap: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


class QSARPredictRequest(BaseSchema):
    smiles_list: List[str] = Field(..., description="SMILES list (max 100)")
    targets: List[str] = Field(..., description="Target names (e.g., EGFR, BRAF)")


class QSARTargetPrediction(BaseSchema):
    probability: float
    std: Optional[float] = None
    ci_low: Optional[float] = None
    ci_high: Optional[float] = None
    in_domain: Optional[bool] = None
    similarity: Optional[float] = None
    calibrated: bool
    active: bool


class QSARCompoundResult(BaseSchema):
    smiles: str
    predictions: Dict[str, QSARTargetPrediction]
    error: Optional[str] = None


class QSARPredictResponse(BaseSchema):
    results: List[QSARCompoundResult]


class QSARTargetInfo(BaseSchema):
    target: str
    model_name: str
    version: str
    metrics: Optional[Dict[str, Any]] = None


class ConformerRequest(BaseSchema):
    smiles: str
    n_conformers: int = Field(5, ge=1, le=50)
    optimize: bool = True


class ConformerResponse(BaseSchema):
    pdb_strings: List[str]
    energies: List[float]


class OverlayRequest(BaseSchema):
    smiles_list: List[str]
    reference_idx: int = Field(0, ge=0)


class OverlayResponse(BaseSchema):
    aligned_pdb_strings: List[str]


class NetworkNode(BaseSchema):
    data: Dict[str, Any]


class NetworkEdge(BaseSchema):
    data: Dict[str, Any]


class CompoundTargetNetworkRequest(BaseSchema):
    compound_ids: List[UUID] = Field(default_factory=list, description="Compound UUIDs (GraphEdge source_entity_id)")
    target_ids: List[UUID] = Field(default_factory=list, description="Target Feature UUIDs (GraphEdge target_entity_id)")
    filters: Dict[str, Any] = Field(default_factory=dict, description="Filters: ic50_range{min_nm,max_nm}, activity_type")


class CompoundTargetNetworkResponse(BaseSchema):
    nodes: List[NetworkNode] = Field(default_factory=list)
    edges: List[NetworkEdge] = Field(default_factory=list)
    meta: Dict[str, Any] = Field(default_factory=dict)


class DoseResponseFitRequest(BaseSchema):
    concentrations: List[float]
    responses: List[float]
    model: str = Field("4PL", description="3PL, 4PL, bayesian_4pl")
    concentration_unit: str = Field("nM", description="Unit label for concentrations (metadata)")
    response_unit: str = Field("percent_inhibition", description="Unit label for responses (metadata)")


class DoseResponseFitResponse(BaseSchema):
    ec50: float
    ec50_ci: Optional[Tuple[float, float]] = None
    hill_slope: float
    top: float
    bottom: Optional[float] = None
    r_squared: Optional[float] = None
    curve_x: List[float]
    curve_y: List[float]
    ci_lower: Optional[List[float]] = None
    ci_upper: Optional[List[float]] = None
    warnings: List[str] = Field(default_factory=list)


class DoseResponseCompareRequest(BaseSchema):
    fits: List[DoseResponseFitRequest]


class DoseResponseCompareResponse(BaseSchema):
    results: List[DoseResponseFitResponse]
    comparison_table: Dict[str, Any]


class TimeseriesAnalyzeRequest(BaseSchema):
    values: List[float]
    timepoints: List[float]
    smooth_method: Optional[str] = Field(None, description="savgol or lowess")
    detect_changepoints: bool = False
    changepoint_threshold: float = 2.0


class TimeseriesAnalyzeResponse(BaseSchema):
    slope: float
    pvalue: float
    direction: str
    smoothed_values: Optional[List[float]] = None
    changepoint_indices: Optional[List[int]] = None


class TrajectoryCompareRequest(BaseSchema):
    series: List[TimeseriesAnalyzeRequest]


class TrajectoryCompareResponse(BaseSchema):
    results: List[TimeseriesAnalyzeResponse]
    comparison_table: Dict[str, Any]
    cluster_labels: Optional[List[int]] = None


class BiomarkerFeature(BaseSchema):
    feature: str
    avg_rank: float
    methods: int


class BiomarkerDiscoverRequest(BaseSchema):
    experiment_id: UUID
    group1_samples: List[str] = Field(..., description="Sample IDs for group 1")
    group2_samples: List[str] = Field(..., description="Sample IDs for group 2")
    methods: List[str] = Field(
        default_factory=lambda: ["statistical", "stability", "importance"],
        description="Discovery methods to run",
    )
    fdr_threshold: float = Field(0.05, ge=0.0, le=1.0, description="FDR threshold for statistical filtering")


class BiomarkerDiscoverResponse(BaseSchema):
    consensus_ranking: List[BiomarkerFeature]
    method_results: Dict[str, Any]


class AlertCheckRequest(BaseSchema):
    smiles: str
    filters: Optional[List[str]] = None


class AlertBatchRequest(BaseSchema):
    smiles_list: List[str] = Field(..., description="SMILES list (max 100)")
    filters: Optional[List[str]] = None


class AlertResultSchema(BaseSchema):
    alert_type: str
    pattern_name: str
    description: str
    severity: str
    matched_smarts: str


class AlertCheckResponse(BaseSchema):
    smiles: str
    is_clean: bool
    traffic_light: str
    alert_count: int
    alerts: List[AlertResultSchema]
    summary: Dict[str, int]
    error: Optional[str] = None


class AlertBatchResponse(BaseSchema):
    results: List[AlertCheckResponse]
    total_checked: int
    clean_count: int
    flagged_count: int


class FilterInfo(BaseSchema):
    name: str
    pattern_count: int
    description: str


class CatalogDatasetsResponse(BaseSchema):
    """Response wrapper for catalog dataset listings."""

    datasets: List[CatalogDataset]
    total: int


# ============================================================================
# Subscription schemas
# ============================================================================


class SubscriptionCreate(BaseSchema):
    """Create a repository subscription."""

    name: str
    repository_source: Literal["GEO", "PRIDE", "MetaboLights", "MW", "all"]
    query_params: Optional[dict] = None
    notify_email: bool = False
    notify_in_app: bool = True


class SubscriptionUpdate(BaseSchema):
    """Update an existing subscription."""

    name: Optional[str] = None
    notify_email: Optional[bool] = None
    notify_in_app: Optional[bool] = None
    is_active: Optional[bool] = None
    query_params: Optional[dict] = None


class SubscriptionResponse(BaseSchema):
    """Subscription response."""

    id: UUID
    user_id: Optional[UUID] = None
    name: str
    repository_source: str
    query_params: Optional[dict] = None
    notify_email: bool
    notify_in_app: bool
    is_active: bool
    last_checked: Optional[datetime] = None
    created_at: datetime
    updated_at: datetime


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
    hit_flag: bool
    hit_category: Optional[str] = None


# ============================================================================
# SAR schemas
# ============================================================================


class TargetResponse(BaseSchema):
    target: str
    compound_count: int = 0


class CompoundActivityResponse(BaseSchema):
    compound_id: str
    smiles: str
    ic50: Optional[float] = None
    units: Optional[str] = None
    assay_name: Optional[str] = None
    result_id: Optional[str] = None


class ActivityCliffResponse(BaseSchema):
    compound_1: str
    smiles_1: str
    activity_1: float
    compound_2: str
    smiles_2: str
    activity_2: float
    similarity: float
    fold_change: float
    assay_id: Optional[str] = None


class ScaffoldSummary(BaseSchema):
    scaffold_smiles: str
    compound_count: int


class SARGridRequest(BaseSchema):
    compound_ids: List[str]
    core_smarts: Optional[str] = None  # Auto-detect if None
    x_axis: str = "R1"
    y_axis: str = "R2"


class SARGridResponse(BaseSchema):
    matrix: List[List[Optional[float]]]  # 2D array of activity values
    row_labels: List[str]  # R-group strings for rows
    col_labels: List[str]  # R-group strings for columns
    core_smiles: str
    total_compounds: int


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


class SignatureStatusUpdate(BaseSchema):
    """Schema for updating signature validation status."""

    status: Literal["pending", "approved", "rejected"]
    reviewer_notes: Optional[str] = None


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
    validation_status: Optional[str] = None


# ============================================================================
# Signature explainability schemas
# ============================================================================


class FeatureContribution(BaseModel):
    """Per-feature contribution to signature match score."""

    feature_name: str
    matched_to: Optional[str] = None
    direction_expected: Optional[str] = None
    direction_actual: Optional[str] = None
    weight: float
    contribution: float
    match_type: str
    direction_match: str


class SignatureExplanationResponse(BaseModel):
    """Explanation of why a signature matches a dataset."""

    signature_id: UUID
    signature_name: str
    dataset_id: UUID
    dataset_name: str
    total_score: float
    direction_concordance: float
    top_positive: List[FeatureContribution]
    top_negative: List[FeatureContribution]
    all_contributions: List[FeatureContribution]


# ============================================================================
# Reports schemas
# ============================================================================


class ReportRequest(BaseModel):
    """Request to generate a narrative report."""

    entity_type: str
    entity_id: str
    format: Literal["html", "pdf"] = "html"


class ReportResponse(BaseModel):
    """Response containing report location."""

    file_path: str
    download_url: Optional[str] = None


# ============================================================================
# Data Quality schemas
# ============================================================================


class DatasetQualityReport(BaseModel):
    """Quality report for a dataset."""

    dataset_id: UUID
    dataset_name: str
    score: float
    status: str
    issues: List[str]
    metrics: Dict[str, float]


class QualitySummary(BaseModel):
    """Summary of dataset quality distribution."""

    total: int
    high: int
    medium: int
    low: int


# ============================================================================
# Protocol diff/deviation schemas
# ============================================================================


class ProtocolDiff(BaseModel):
    protocol_id: UUID
    other_id: UUID
    added_steps: List[Dict[str, Any]]
    removed_steps: List[Dict[str, Any]]
    changed_steps: List[Dict[str, Any]]
    materials_added: List[Dict[str, Any]]
    materials_removed: List[Dict[str, Any]]
    parameters_changed: List[Dict[str, Any]]


class ProtocolHistoryItem(BaseModel):
    protocol_id: UUID
    version: int
    parent_id: Optional[UUID] = None
    name: str


class DeviationReport(BaseModel):
    experiment_id: UUID
    protocol_id: UUID
    protocol_name: str
    protocol_version: int
    deviations: List[Any]


# ============================================================================
# HTS QC schemas
# ============================================================================


class PlateQCSummary(BaseModel):
    campaign_id: UUID
    total_wells: int
    hit_rate: float
    z_prime: Optional[float]
    hits: int
    pos_controls: int
    neg_controls: int


class WellData(BaseModel):
    well_position: Optional[str] = None
    normalized_value: Optional[float] = None
    z_score: Optional[float] = None
    hit_flag: Optional[bool] = None
    compound_id: Optional[UUID] = None
    result_id: Optional[str] = None


class HitCompound(BaseModel):
    result_id: str
    compound_id: UUID
    well_position: Optional[str] = None
    normalized_value: Optional[float] = None
    z_score: Optional[float] = None


# ============================================================================
# Cross-omics pathway schemas
# ============================================================================


class OmicsFeatureSet(BaseModel):
    transcriptomics: List[str] = Field(default_factory=list)
    proteomics: List[str] = Field(default_factory=list)
    metabolomics: List[str] = Field(default_factory=list)
    lipidomics: List[str] = Field(default_factory=list)


class ConvergentPathway(BaseModel):
    pathway_id: str
    name: str
    source: str
    adjusted_p_value: float
    enrichment_ratio: float
    matched_by_omics: Dict[str, List[str]]


class CrossOmicsEnrichmentResult(BaseModel):
    program_id: UUID
    pathways: List[ConvergentPathway]
    features: OmicsFeatureSet


class PathwayFeatures(BaseModel):
    pathway_id: str
    pathway_name: str
    source: str
    features_by_omics: Dict[str, List[str]]


# ============================================================================
# Pathway maps (KEGG-only MVP) schemas
# ============================================================================


class PathwayNodeSchema(BaseModel):
    id: str
    name: str
    type: str
    x: float
    y: float
    kegg_ids: List[str] = Field(default_factory=list)


class PathwayEdgeSchema(BaseModel):
    source: str
    target: str
    type: str
    subtype: Optional[str] = None
    style: str = "solid"
    color: str = "#9E9E9E"


class PathwayStructureResponse(BaseModel):
    pathway_id: str
    name: str
    nodes: List[PathwayNodeSchema] = Field(default_factory=list)
    edges: List[PathwayEdgeSchema] = Field(default_factory=list)
    organism: str = ""


class PathwayOverlayRequest(BaseModel):
    pathway_id: str
    expression_data: Dict[str, float] = Field(default_factory=dict)
    colormap: str = "RdBu_r"
    vmin: float = -2.0
    vmax: float = 2.0


class PathwayOverlayNodeSchema(BaseModel):
    node_id: str
    gene_symbol: str
    value: float
    color: str
    label: str


class PathwayOverlayResponse(BaseModel):
    pathway_id: str
    overlays: List[PathwayOverlayNodeSchema] = Field(default_factory=list)


class PathwaySearchResult(BaseModel):
    pathway_id: str
    name: str
    organism: str
    gene_count: int = 0


class PathwayExpressionResponse(BaseModel):
    dataset_id: UUID
    gene_expression: Dict[str, float] = Field(default_factory=dict)




# ============================================================================
# MOA inference schemas
# ============================================================================


class EvidenceContribution(BaseModel):
    feature_name: str
    value: float
    weight: float


class MOACandidate(BaseModel):
    candidate_id: str
    type: str
    probability: float
    probability_ci: Optional[Tuple[float, float]] = None
    rank: int
    contributions: List[EvidenceContribution]


class MOAInferenceRequest(BaseModel):
    compound_id: UUID
    dataset_ids: List[UUID] = Field(default_factory=list)
    method: Literal["fusion", "bayesian"] = "fusion"


class MOAInferenceResult(BaseModel):
    compound_id: UUID
    candidates: List[MOACandidate]


# ============================================================================
# Analysis / Bayesian inference schemas
# ============================================================================


class BayesianDoseResponseRequest(BaseModel):
    concentrations: List[float]
    responses: List[float]
    prior_ec50: Optional[float] = None
    likelihood: str = "normal"
    include_diagnostics: bool = False


class BayesianDoseResponseResponse(BaseModel):
    ec50_mean: float
    ec50_ci: Tuple[float, float]
    hill_slope: float
    diagnostics: Optional[Dict[str, Any]] = None


# ============================================================================
# Multi-Omics Visualization schemas
# ============================================================================


class AlluvialRequest(BaseModel):
    """Request for alluvial diagram data."""
    
    dataset_ids: List[str]  # UUIDs as strings


class AlluvialNode(BaseModel):
    """Node in alluvial diagram."""
    
    id: str
    label: str
    color: str


class AlluvialLink(BaseModel):
    """Link in alluvial diagram."""
    
    source: int
    target: int
    value: int


class AlluvialResponse(BaseModel):
    """Response for alluvial diagram data."""
    
    nodes: List[AlluvialNode]
    links: List[AlluvialLink]


class UpSetRequest(BaseModel):
    """Request for UpSet plot data."""
    
    dataset_ids: List[str]


class UpSetSet(BaseModel):
    """Set definition for UpSet plot."""
    
    id: str
    name: str
    omics_type: str
    color: str
    size: int


class UpSetIntersection(BaseModel):
    """Intersection data for UpSet plot."""
    
    sets: List[str]  # dataset IDs in intersection
    count: int
    bitmask: int


class UpSetMatrixItem(BaseModel):
    """Matrix item for UpSet plot."""
    
    key: str
    label: str
    presence: List[int]
    dataset_count: int


class UpSetResponse(BaseModel):
    """Response for UpSet plot data."""
    
    sets: List[UpSetSet]
    intersections: List[UpSetIntersection]
    matrix: List[UpSetMatrixItem]  # List of matrix items with presence vectors


# ============================================================================
# Compound Ranking schemas
# ============================================================================


class RankingRequest(BaseModel):
    """Request for compound ranking."""
    
    compound_ids: List[str]  # UUIDs as strings
    weights: Optional[Dict[str, float]] = None  # Override preset
    preset: Optional[str] = None  # "balanced", "potency_first", etc.
    include_pareto: bool = True


class ObjectiveScoreSchema(BaseModel):
    """Individual objective score."""
    
    name: str
    raw_value: Optional[float]
    normalized: float
    weight: float
    confidence: Optional[float] = None


class CompoundRankingSchema(BaseModel):
    """Compound ranking result."""
    
    compound_id: str
    smiles: str
    objectives: List[ObjectiveScoreSchema]
    weighted_score: float
    pareto_rank: Optional[int]
    rank: int


class RankingResponse(BaseModel):
    """Response for compound ranking."""
    
    rankings: List[CompoundRankingSchema]
    pareto_front: List[CompoundRankingSchema]
    total_compounds: int
    skipped_compounds: int  # Missing potency or invalid SMILES


class ParetoRequest(BaseModel):
    """Request for Pareto plot data."""
    
    compound_ids: List[str]
    x_objective: str = "potency"
    y_objective: str = "liability_aggregate"  # or specific: "herg", "alerts"


class ParetoPoint(BaseModel):
    """Point in Pareto plot."""
    
    compound_id: str
    smiles: str
    x_value: float
    y_value: float
    pareto_rank: int
    is_frontier: bool


class ParetoResponse(BaseModel):
    """Response for Pareto plot data."""
    
    points: List[ParetoPoint]
    x_label: str
    y_label: str
    frontier_ids: List[str]


class RankingPreset(BaseModel):
    """Ranking weight preset."""
    
    name: str
    weights: Dict[str, float]
    description: str


# Model Monitoring Schemas
class MonitoringLogRequest(BaseModel):
    """Request to log a model prediction for monitoring."""
    
    model_id: UUID
    model_version: str
    prediction_id: UUID
    prediction: float
    input_hash: Optional[str] = None
    feature_summary: Dict[str, Any]


class MonitoringFeedbackRequest(BaseModel):
    """Request to provide ground truth feedback for a prediction."""
    
    prediction_id: UUID
    ground_truth: float


class DriftReportSchema(BaseModel):
    """Drift analysis report for a model."""
    
    model_id: UUID
    model_name: str
    status: str
    psi_scores: Dict[str, float]
    fp_aggregate_drift: Dict[str, float]
    window_hours: int
    n_predictions: int
    computed_at: datetime


class CalibrationReportSchema(BaseModel):
    """Calibration analysis report for a classification model."""
    
    model_id: UUID
    model_name: str
    status: str
    ece: Optional[float]
    n_predictions_with_truth: int
    computed_at: datetime


class HealthReportSchema(BaseModel):
    """Overall health report combining drift and calibration analysis."""
    
    model_id: UUID
    model_name: str
    overall_status: str
    drift_status: str
    calibration_status: str
    psi_max: Optional[float]
    ece: Optional[float]
    last_checked: datetime


# Activity & Notifications Schemas
class ActivityEventSchema(BaseModel):
    """Activity event schema for API responses."""
    
    id: UUID
    event_type: str
    actor_id: Optional[UUID] = None
    target_type: str
    target_id: UUID
    target_name: str
    program_id: Optional[UUID] = None
    metadata: Optional[Dict[str, Any]] = None
    created_at: datetime


class NotificationSchema(BaseModel):
    """Notification schema for API responses."""
    
    id: UUID
    user_id: UUID
    activity_event_id: Optional[UUID] = None
    notification_type: str
    title: Optional[str] = None
    message: Optional[str] = None
    is_read: bool
    created_at: datetime
    activity_event: Optional[ActivityEventSchema] = None


class NotificationCountSchema(BaseModel):
    """Notification count schema."""
    
    unread_count: int


# Comment schemas
class CommentCreate(BaseModel):
    """Schema for creating a comment."""
    
    entity_type: str
    entity_id: UUID
    content: str
    parent_id: Optional[UUID] = None


class CommentUpdate(BaseModel):
    """Schema for updating a comment."""
    
    content: str


class CommentResponse(BaseModel):
    """Schema for comment response."""
    
    id: UUID
    entity_type: str
    entity_id: UUID
    content: str
    author: str
    created_at: datetime
    updated_at: Optional[datetime] = None
    replies: List[dict] = []


# ============================================================================
# Entity Sharing schemas
# ============================================================================

class EntityShareCreate(BaseModel):
    """Schema for creating entity shares."""
    
    shared_with_user_id: Optional[UUID] = None
    shared_with_team_id: Optional[UUID] = None
    permission: Literal["view", "edit", "admin"] = "view"
    
    @field_validator("permission")
    @classmethod
    def validate_permission(cls, v):
        if v not in ["view", "edit", "admin"]:
            raise ValueError("Permission must be view, edit, or admin")
        return v


class EntityShareResponse(BaseModel):
    """Schema for entity share response."""
    
    id: UUID
    entity_type: str
    entity_id: UUID
    shared_with_user_id: Optional[UUID] = None
    shared_with_team_id: Optional[UUID] = None
    permission: str
    shared_by_id: UUID
    created_at: datetime


class EntityShareList(BaseModel):
    """Schema for list of entity shares."""
    
    shares: List[EntityShareResponse]


# ============================================================================
# Entity Review schemas
# ============================================================================

class EntityReviewCreate(BaseModel):
    """Schema for creating entity reviews."""
    
    pass  # No additional fields needed - entity info comes from URL


class EntityReviewResponse(BaseModel):
    """Schema for entity review response."""
    
    id: UUID
    entity_type: str
    entity_id: UUID
    reviewer_id: UUID
    status: Literal["draft", "submitted", "in_review", "approved", "rejected", "changes_requested"]
    comments: Optional[str] = None
    reviewed_at: datetime


class EntityReviewList(BaseModel):
    """Schema for list of entity reviews."""
    
    reviews: List[EntityReviewResponse]


class EntityReviewDecision(BaseModel):
    """Schema for review decision."""
    
    decision: Literal["approved", "rejected", "changes_requested"]
    comment: str = Field(..., min_length=1, description="Review comment/feedback")


class EntityReviewAssign(BaseModel):
    """Schema for reviewer assignment."""
    
    reviewer_id: UUID


class EntityReviewStatusInfo(BaseModel):
    """Schema for review status information."""
    
    valid_statuses: List[str]
    valid_transitions: Dict[str, List[str]]
    status_descriptions: Dict[str, str]


# ============================================================================
# Team Management schemas
# ============================================================================

class TeamCreate(BaseModel):
    """Schema for creating teams."""
    
    name: str = Field(..., min_length=1, max_length=255)
    description: Optional[str] = None


class TeamResponse(BaseModel):
    """Schema for team response."""
    
    id: UUID
    name: str
    description: Optional[str] = None
    created_at: datetime
    member_count: int
    user_role: str  # owner, admin, member, viewer


class TeamListResponse(BaseModel):
    """Schema for list of teams."""
    
    teams: List[TeamResponse]


class TeamMemberResponse(BaseModel):
    """Schema for team member response."""
    
    id: UUID
    user_id: UUID
    username: str
    email: str
    role: str
    joined_at: datetime


# ============================================================================
# Generative Chemistry schemas
# ============================================================================


class PropertyConstraintSchema(BaseModel):
    """Schema for property constraint in optimization."""
    
    name: str = Field(..., description="Property name (e.g., 'logP', 'herg', 'mw')")
    min_value: Optional[float] = Field(None, description="Minimum allowed value")
    max_value: Optional[float] = Field(None, description="Maximum allowed value")
    target_value: Optional[float] = Field(None, description="Target value for optimization")
    weight: float = Field(1.0, ge=0.0, description="Importance weight")
    
    @field_validator('name')
    @classmethod
    def validate_name(cls, v):
        if not v or not v.strip():
            raise ValueError("Property name cannot be empty")
        return v.strip().lower()


class MoleculeSchema(BaseModel):
    """Schema for generated molecule with properties."""
    
    smiles: str = Field(..., description="Generated SMILES string")
    properties: Dict[str, float] = Field(default_factory=dict, description="Predicted properties")
    score: Optional[float] = Field(None, description="Overall optimization score")
    step: Optional[int] = Field(None, description="Generation step (for interpolation)")
    iteration: Optional[int] = Field(None, description="Optimization iteration")


class SampleRequest(BaseModel):
    """Request schema for random sampling."""
    
    n_samples: int = Field(10, ge=1, le=100, description="Number of molecules to generate")
    temperature: float = Field(1.0, ge=0.1, le=2.0, description="Sampling temperature")
    max_length: int = Field(100, ge=20, le=200, description="Maximum SMILES length")


class SampleResponse(BaseModel):
    """Response schema for random sampling."""
    
    molecules: List[MoleculeSchema] = Field(..., description="Generated molecules")
    count: int = Field(..., description="Number of molecules generated")
    model_info: Dict = Field(..., description="Model information used")


class InterpolateRequest(BaseModel):
    """Request schema for molecular interpolation."""
    
    smiles_start: str = Field(..., description="Starting molecule SMILES")
    smiles_end: str = Field(..., description="Ending molecule SMILES")
    steps: int = Field(10, ge=2, le=50, description="Number of interpolation steps")
    interpolation_type: str = Field("linear", description="Interpolation method")
    
    @field_validator('smiles_start', 'smiles_end')
    @classmethod
    def validate_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError("SMILES cannot be empty")
        return v.strip()


class InterpolateResponse(BaseModel):
    """Response schema for molecular interpolation."""
    
    molecules: List[MoleculeSchema] = Field(..., description="Interpolated molecules")
    start_smiles: str = Field(..., description="Starting molecule")
    end_smiles: str = Field(..., description="Ending molecule")
    steps: int = Field(..., description="Number of steps")


class OptimizeRequest(BaseModel):
    """Request schema for property-guided optimization."""
    
    seed_smiles: str = Field(..., description="Starting molecule SMILES")
    constraints: List[PropertyConstraintSchema] = Field(..., description="Property constraints")
    n_iterations: int = Field(100, ge=1, le=1000, description="Number of optimization iterations")
    n_samples_per_iter: int = Field(10, ge=1, le=50, description="Samples per iteration")
    learning_rate: float = Field(0.1, ge=0.01, le=1.0, description="Optimization learning rate")
    temperature: float = Field(1.0, ge=0.1, le=2.0, description="Initial sampling temperature")
    
    @field_validator('seed_smiles')
    @classmethod
    def validate_seed_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError("Seed SMILES cannot be empty")
        return v.strip()
    
    @field_validator('constraints')
    @classmethod
    def validate_constraints(cls, v):
        if not v:
            raise ValueError("At least one constraint must be specified")
        return v


class OptimizeResponse(BaseModel):
    """Response schema for property-guided optimization."""
    
    optimized: List[MoleculeSchema] = Field(..., description="Optimized molecules")
    seed_smiles: str = Field(..., description="Starting molecule")
    seed_properties: Dict[str, float] = Field(..., description="Properties of seed molecule")
    best_score: float = Field(..., description="Best optimization score achieved")
    iterations_completed: int = Field(..., description="Number of iterations completed")


class ScaffoldHopRequest(BaseModel):
    """Request schema for scaffold hopping."""
    
    smiles: str = Field(..., description="Input molecule SMILES")
    n_analogs: int = Field(20, ge=1, le=50, description="Number of analogs to generate")
    preserve_scaffold: bool = Field(True, description="Whether to preserve the scaffold")
    similarity_threshold: float = Field(0.7, ge=0.1, le=0.9, description="Similarity threshold")
    
    @field_validator('smiles')
    @classmethod
    def validate_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError("SMILES cannot be empty")
        return v.strip()


class ScaffoldHopResponse(BaseModel):
    """Response schema for scaffold hopping."""
    
    scaffold: Optional[str] = Field(None, description="Extracted scaffold SMILES")
    analogs: List[MoleculeSchema] = Field(..., description="Generated analogs")
    input_smiles: str = Field(..., description="Original input molecule")
    n_generated: int = Field(..., description="Number of analogs generated")


class GenerativeModelInfo(BaseModel):
    """Schema for generative model information."""
    
    name: str = Field(..., description="Model name")
    version: str = Field(..., description="Model version")
    latent_dim: int = Field(..., description="Latent space dimensionality")
    vocab_size: int = Field(..., description="Vocabulary size")
    status: str = Field(..., description="Model status (loaded/not_found/error)")
    created_at: Optional[str] = Field(None, description="Model creation timestamp")
    description: Optional[str] = Field(None, description="Model description")


class GenerativeModelsResponse(BaseModel):
    """Response schema for listing available models."""
    
    models: List[GenerativeModelInfo] = Field(..., description="Available generative models")
    count: int = Field(..., description="Number of models")
    default_model: Optional[str] = Field(None, description="Default model name")


class ErrorResponse(BaseModel):
    """Schema for error responses."""
    
    error: str = Field(..., description="Error message")
    detail: Optional[str] = Field(None, description="Detailed error information")
    error_type: str = Field(..., description="Type of error")
    request_id: Optional[str] = Field(None, description="Request identifier for tracking")


# ============================================================================
# Imaging schemas
# ============================================================================


class ImageUploadRequest(BaseModel):
    """Request for uploading microscopy image."""
    
    well_id: Optional[UUID] = Field(None, description="Well ID to associate with image")
    channel: str = Field(..., description="Channel name (e.g., DAPI, GFP, RFP)")
    z_slice: int = Field(0, description="Z-stack slice index")
    timepoint: int = Field(0, description="Time series index")
    pixel_size_um: Optional[float] = Field(None, description="Pixel size in microns")


class ImageUploadResponse(BaseModel):
    """Response for image upload."""
    
    image_id: UUID = Field(description="Unique image identifier")
    image_path: str = Field(description="Storage path of uploaded image")
    width: int = Field(description="Image width in pixels")
    height: int = Field(description="Image height in pixels")
    channel: str = Field(description="Channel name")
    well_id: Optional[UUID] = Field(description="Associated well ID")
    message: str = Field(description="Success message")
    
    model_config = ConfigDict(from_attributes=True)


class SegmentationRequest(BaseModel):
    """Request for cell segmentation."""
    
    image_id: UUID = Field(description="Image to segment")
    model_name: Optional[str] = Field("cyto", description="CellPose model type")
    diameter: Optional[float] = Field(30.0, description="Expected cell diameter in pixels")
    channels: Optional[List[int]] = Field([0, 0], description="Channel configuration [cytoplasm, nucleus]")
    extract_features: bool = Field(True, description="Whether to extract morphological features")


class SegmentationResponse(BaseModel):
    """Response for cell segmentation."""
    
    segmentation_id: UUID = Field(description="Unique segmentation identifier")
    image_id: UUID = Field(description="Source image ID")
    cell_count: int = Field(description="Number of cells detected")
    model_name: str = Field(description="Model used for segmentation")
    mask_path: str = Field(description="Storage path of segmentation mask")
    features_extracted: bool = Field(description="Whether features were extracted")
    processing_time_seconds: float = Field(description="Processing time")
    
    model_config = ConfigDict(from_attributes=True)


class BatchSegmentationRequest(BaseModel):
    """Request for batch segmentation."""
    
    image_ids: List[UUID] = Field(description="List of images to segment")
    model_name: Optional[str] = Field("cyto", description="CellPose model type")
    diameter: Optional[float] = Field(30.0, description="Expected cell diameter in pixels")
    channels: Optional[List[int]] = Field([0, 0], description="Channel configuration")
    extract_features: bool = Field(True, description="Whether to extract features")


class BatchSegmentationResponse(BaseModel):
    """Response for batch segmentation."""
    
    task_id: str = Field(description="Celery task ID for tracking")
    image_count: int = Field(description="Number of images queued")
    status: str = Field(description="Task status")
    message: str = Field(description="Status message")


class ImageMetadataResponse(BaseModel):
    """Response for image metadata."""
    
    image_id: UUID = Field(description="Image identifier")
    well_id: Optional[UUID] = Field(description="Associated well ID")
    channel: str = Field(description="Channel name")
    z_slice: int = Field(description="Z-stack slice index")
    timepoint: int = Field(description="Time series index")
    width: int = Field(description="Image width in pixels")
    height: int = Field(description="Image height in pixels")
    bit_depth: int = Field(description="Image bit depth")
    pixel_size_um: Optional[float] = Field(description="Pixel size in microns")
    image_path: str = Field(description="Storage path")
    metadata: Optional[Dict[str, Any]] = Field(description="Additional metadata")
    acquired_at: Optional[datetime] = Field(description="Image acquisition time")
    created_at: datetime = Field(description="Record creation time")
    
    model_config = ConfigDict(from_attributes=True)


class SegmentationResultResponse(BaseModel):
    """Response for segmentation results."""
    
    segmentation_id: UUID = Field(description="Segmentation identifier")
    image_id: UUID = Field(description="Source image ID")
    model_name: str = Field(description="Model used")
    model_version: Optional[str] = Field(description="Model version")
    cell_count: int = Field(description="Number of cells detected")
    mask_path: str = Field(description="Segmentation mask path")
    parameters: Optional[Dict[str, Any]] = Field(description="Segmentation parameters")
    confidence_score: Optional[float] = Field(description="Overall confidence score")
    created_at: datetime = Field(description="Segmentation time")
    
    model_config = ConfigDict(from_attributes=True)


class CellFeaturesResponse(BaseModel):
    """Response for cell features."""
    
    segmentation_id: UUID = Field(description="Source segmentation ID")
    cell_count: int = Field(description="Number of cells")
    features: List[Dict[str, Any]] = Field(description="List of cell features")
    
    model_config = ConfigDict(from_attributes=True)


class WellSummaryResponse(BaseModel):
    """Response for well-level summary."""
    
    well_id: UUID = Field(description="Well identifier")
    image_count: int = Field(description="Number of images in well")
    total_cell_count: int = Field(description="Total cells across all images")
    channels: List[str] = Field(description="Available channels")
    aggregated_features: Dict[str, Any] = Field(description="Aggregated morphology features")
    summary_metrics: Dict[str, Any] = Field(description="High-level well metrics")
    
    model_config = ConfigDict(from_attributes=True)


# ============================================================================
# Flow Cytometry schemas
# ============================================================================

class BaseFlowSchema(BaseSchema):
    """Base schema with common configuration for flow cytometry."""
    pass


class FlowCytometryDatasetCreate(BaseFlowSchema):
    """Schema for creating a flow cytometry dataset."""
    
    dataset_id: Optional[UUID] = None
    sample_id: Optional[str] = None
    sample_volume_ul: Optional[float] = Field(None, gt=0)
    dilution_factor: Optional[float] = Field(None, gt=0)
    staining_protocol: Optional[str] = None


class FlowCytometryDatasetResponse(BaseFlowSchema):
    """Schema for flow cytometry dataset response."""
    
    id: UUID
    dataset_id: UUID
    events_parquet_path: str
    file_size_bytes: Optional[int] = None
    n_events: Optional[int] = None
    n_parameters: Optional[int] = None
    
    # Acquisition metadata
    acquisition_date: Optional[datetime] = None
    cytometer_model: Optional[str] = None
    cytometer_serial: Optional[str] = None
    acquisition_software: Optional[str] = None
    acquisition_settings: Optional[Dict[str, Any]] = None
    
    # Sample information
    sample_id: Optional[str] = None
    sample_volume_ul: Optional[float] = None
    dilution_factor: Optional[float] = None
    staining_protocol: Optional[str] = None
    
    # Processing status
    processing_status: str
    processing_log: Optional[str] = None
    
    # Timestamps
    ingested_at: datetime
    processed_at: Optional[datetime] = None


class FlowCytometryParameterResponse(BaseFlowSchema):
    """Schema for flow cytometry parameter response."""
    
    id: UUID
    flow_dataset_id: UUID
    parameter_index: int
    parameter_name: str
    parameter_short_name: Optional[str] = None
    
    # Channel configuration
    detector: Optional[str] = None
    excitation_wavelength: Optional[int] = None
    emission_wavelength: Optional[int] = None
    fluorophore: Optional[str] = None
    
    # Data range and scaling
    data_range: Optional[int] = None
    amplifier_gain: Optional[float] = None
    voltage: Optional[float] = None
    
    # Statistics
    min_value: Optional[float] = None
    max_value: Optional[float] = None


class GateCreate(BaseFlowSchema):
    """Schema for creating a gate."""
    
    gate_name: str = Field(..., min_length=1, max_length=200)
    gate_type: Literal["polygon", "rectangle", "quadrant"]
    gate_definition: Dict[str, Any]
    x_parameter_id: UUID
    y_parameter_id: Optional[UUID] = None
    parent_gate_id: Optional[UUID] = None
    boolean_operator: Optional[Literal["AND", "OR", "NOT"]] = None
    operand_gate_ids: Optional[List[UUID]] = None
    
    @field_validator('gate_definition')
    @classmethod
    def validate_gate_definition(cls, v: Dict[str, Any], info) -> Dict[str, Any]:
        """Validate gate definition based on gate type."""
        if 'gate_type' not in info.data:
            return v
            
        gate_type = info.data['gate_type']
        
        if gate_type == "polygon":
            if "vertices" not in v:
                raise ValueError("Polygon gate requires 'vertices' in gate_definition")
            vertices = v["vertices"]
            if not isinstance(vertices, list) or len(vertices) < 3:
                raise ValueError("Polygon gate requires at least 3 vertices")
            for vertex in vertices:
                if not isinstance(vertex, list) or len(vertex) != 2:
                    raise ValueError("Each vertex must be [x, y] coordinate pair")
                    
        elif gate_type == "rectangle":
            required_keys = {"x_min", "x_max", "y_min", "y_max"}
            if not required_keys.issubset(v.keys()):
                raise ValueError(f"Rectangle gate requires keys: {required_keys}")
            if v["x_min"] >= v["x_max"] or v["y_min"] >= v["y_max"]:
                raise ValueError("Rectangle bounds: min values must be less than max values")
                
        elif gate_type == "quadrant":
            required_keys = {"x_threshold", "y_threshold"}
            if not required_keys.issubset(v.keys()):
                raise ValueError(f"Quadrant gate requires keys: {required_keys}")
            if "quadrant" in v and v["quadrant"] not in ["Q1", "Q2", "Q3", "Q4"]:
                raise ValueError("Quadrant must be one of: Q1, Q2, Q3, Q4")
                
        return v


class GateUpdate(BaseFlowSchema):
    """Schema for updating a gate (all fields optional)."""
    
    gate_name: Optional[str] = Field(None, min_length=1, max_length=200)
    gate_definition: Optional[Dict[str, Any]] = None
    is_active: Optional[bool] = None
    parent_gate_id: Optional[UUID] = None
    boolean_operator: Optional[Literal["AND", "OR", "NOT"]] = None
    operand_gate_ids: Optional[List[UUID]] = None


class GateResponse(BaseFlowSchema):
    """Schema for gate response."""
    
    id: UUID
    flow_dataset_id: UUID
    gate_name: str
    gate_type: str
    gate_definition: Dict[str, Any]
    x_parameter_id: UUID
    y_parameter_id: Optional[UUID] = None
    parent_gate_id: Optional[UUID] = None
    boolean_operator: Optional[str] = None
    operand_gate_ids: Optional[List[UUID]] = None
    is_active: bool
    created_at: datetime


class PopulationResponse(BaseFlowSchema):
    """Schema for population statistics response."""
    
    id: UUID
    flow_dataset_id: UUID
    gate_id: UUID
    event_count: int
    parent_event_count: Optional[int] = None
    percentage_of_parent: Optional[float] = None
    percentage_of_total: Optional[float] = None
    population_name: Optional[str] = None
    cell_type: Optional[str] = None
    phenotype: Optional[str] = None
    parameter_statistics: Optional[Dict[str, Any]] = None
    analysis_software: Optional[str] = None
    analysis_version: Optional[str] = None
    analyzed_at: datetime


class EventsQueryParams(BaseModel):
    """Schema for querying flow cytometry events."""
    
    offset: int = Field(0, ge=0, description="Number of events to skip")
    limit: int = Field(10000, ge=1, le=100000, description="Maximum number of events to return")
    parameters: Optional[List[str]] = Field(None, description="Parameter names to include (all if None)")
    subsample: bool = Field(False, description="Apply random subsampling if dataset is large")
    
    model_config = ConfigDict(from_attributes=True)


class EventsResponse(BaseFlowSchema):
    """Schema for events data response."""
    
    events: List[List[float]]  # 2D array as list of lists
    parameter_names: List[str]
    total_events: int
    offset: int
    limit: int
    subsampled: bool


class UploadResponse(BaseFlowSchema):
    """Schema for FCS file upload response."""
    
    flow_dataset_id: UUID
    dataset_id: UUID
    filename: str
    file_size_bytes: int
    processing_status: str
    message: str


# ============================================================================
# Biophysical Assay Schemas
# ============================================================================


class BaseBiophysicalSchema(BaseSchema):
    """Base schema for biophysical assay schemas."""
    pass


# Upload request schemas
class SPRUploadRequest(BaseBiophysicalSchema):
    """Schema for SPR file upload request."""
    
    experiment_id: Optional[UUID] = None
    compound_id: Optional[UUID] = None
    target_id: Optional[UUID] = None
    target_name: Optional[str] = None


class MSTUploadRequest(BaseBiophysicalSchema):
    """Schema for MST file upload request."""
    
    experiment_id: Optional[UUID] = None
    compound_id: Optional[UUID] = None
    target_id: Optional[UUID] = None
    target_name: Optional[str] = None


class DSCUploadRequest(BaseBiophysicalSchema):
    """Schema for DSC file upload request."""
    
    experiment_id: Optional[UUID] = None
    compound_id: Optional[UUID] = None
    protein_id: Optional[UUID] = None
    protein_name: Optional[str] = None


# Response schemas
class SPRExperimentResponse(BaseBiophysicalSchema):
    """Schema for SPR experiment response."""
    
    id: UUID
    experiment_name: Optional[str]
    compound_id: Optional[UUID]
    target_id: Optional[UUID]
    target_name: Optional[str]
    instrument: Optional[str]
    chip_type: Optional[str]
    ka: Optional[float]
    kd_rate: Optional[float]
    kd_equilibrium: Optional[float]
    rmax: Optional[float]
    chi_squared: Optional[float]
    processing_status: str
    created_at: datetime


class SPRSensorgamResponse(BaseBiophysicalSchema):
    """Schema for SPR sensorgram response."""
    
    id: UUID
    cycle_number: int
    analyte_concentration: float
    time_seconds: List[float]
    response_values: List[float]


class MSTExperimentResponse(BaseBiophysicalSchema):
    """Schema for MST experiment response."""
    
    id: UUID
    experiment_name: Optional[str]
    compound_id: Optional[UUID]
    target_id: Optional[UUID]
    target_name: Optional[str]
    kd_value: Optional[float]
    kd_error: Optional[float]
    hill_coefficient: Optional[float]
    processing_status: str
    created_at: datetime


class MSTDoseResponsePointResponse(BaseBiophysicalSchema):
    """Schema for MST dose response point response."""
    
    id: UUID
    concentration: float
    fnorm: float
    fnorm_error: Optional[float]


class DSCExperimentResponse(BaseBiophysicalSchema):
    """Schema for DSC experiment response."""
    
    id: UUID
    experiment_name: Optional[str]
    compound_id: Optional[UUID]
    protein_id: Optional[UUID]
    protein_name: Optional[str]
    tm_value: Optional[float]
    tm_error: Optional[float]
    delta_h: Optional[float]
    processing_status: str
    created_at: datetime


class DSCScanResponse(BaseBiophysicalSchema):
    """Schema for DSC scan response."""
    
    id: UUID
    scan_number: int
    temperature_celsius: List[float]
    heat_capacity: List[float]


# Fit request schemas
class RefitRequest(BaseBiophysicalSchema):
    """Schema for refitting experiment data."""
    
    model: Optional[str] = Field(None, description="Fitting model: '1:1_langmuir', 'two_state', 'hill', etc.")


# Cross-assay comparison response
class BiophysicalCompareResponse(BaseBiophysicalSchema):
    """Schema for cross-assay biophysical comparison."""
    
    compound_id: UUID
    spr_results: List[SPRExperimentResponse]
    mst_results: List[MSTExperimentResponse]
    dsc_results: List[DSCExperimentResponse]


# Upload response schemas
class BiophysicalUploadResponse(BaseBiophysicalSchema):
    """Schema for biophysical file upload response."""
    
    experiment_id: UUID
    filename: str
    file_size_bytes: int
    processing_status: str
    message: str


# ============================================================================
# Inline Annotation schemas
# ============================================================================


class InlineAnnotationCreate(BaseSchema):
    """Schema for creating an inline annotation."""
    
    entity_type: str = Field(..., description="Type of entity being annotated")
    entity_id: UUID = Field(..., description="ID of the entity being annotated")
    position_type: str = Field(..., description="Type of position (cell, column, row, field, range)")
    position_data: Dict[str, Any] = Field(..., description="Position-specific data")
    content: str = Field(..., description="Annotation content")
    parent_id: Optional[UUID] = Field(None, description="Parent annotation ID for replies")


class InlineAnnotationResponse(BaseSchema):
    """Schema for inline annotation response."""
    
    id: UUID
    entity_type: str
    entity_id: UUID
    position_type: str
    position_data: Dict[str, Any]
    content: str
    status: str
    parent_id: Optional[UUID]
    created_by_id: Optional[UUID]
    resolved_by_id: Optional[UUID]
    created_at: datetime
    updated_at: Optional[datetime]
    resolved_at: Optional[datetime]


class InlineAnnotationReplyCreate(BaseSchema):
    """Schema for creating a reply to an annotation."""
    
    content: str = Field(..., description="Reply content")


class InlineAnnotationListResponse(BaseSchema):
    """Schema for inline annotation list response."""
    
    annotations: List[InlineAnnotationResponse]
    total: int
    open_count: int
    resolved_count: int
