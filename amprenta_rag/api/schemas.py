"""
Pydantic schemas for API request/response models.

These schemas define the structure of API requests and responses,
mapping between domain models and API representations.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Literal, Optional, Tuple
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field, field_validator

from amprenta_rag.models.domain import FeatureType, OmicsType, SignatureDirection


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


class AnnotationCreate(BaseSchema):
    """Schema for creating an annotation/note on an entity."""

    text: str
    annotation_type: Optional[str] = None


# Import PriorConfig from analysis layer to avoid circular imports
from amprenta_rag.analysis.models import PriorConfig


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
