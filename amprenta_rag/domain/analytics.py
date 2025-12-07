from datetime import datetime
from typing import Any, Dict, List, Optional, Set
from uuid import UUID

from pydantic import BaseModel


class DatasetFeatureProfile(BaseModel):
    dataset_id: UUID
    omics_type: str
    diseases: Optional[List[str]]
    matrix: Optional[str]
    features_by_type: Dict[str, Set[str]]
    directions: Optional[Dict[str, str]] = None


class DatasetComparisonResult(BaseModel):
    dataset_id_1: UUID
    dataset_id_2: UUID
    jaccard_by_type: Dict[str, float]  # e.g. {"gene": 0.6, ...}
    overall_score: float
    shared_features: Set[str]
    unique_to_1: Set[str]
    unique_to_2: Set[str]


class DatasetCluster(BaseModel):
    cluster_id: int
    member_dataset_ids: List[UUID]
    linkage_score: Optional[float] = None


class EvidenceReportRequest(BaseModel):
    entity_type: str  # 'program' | 'experiment' | 'dataset' | 'signature'
    id: UUID
    include_sections: Optional[dict[str, bool]] = None


class EvidenceSection(BaseModel):
    title: str
    summary_text: str
    supporting_datasets: Optional[List[UUID]] = None
    key_features: Optional[List[str]] = None
    signatures: Optional[List[UUID]] = None
    references: Optional[List[Any]] = None


class EvidenceReport(BaseModel):
    entity_id: UUID
    entity_type: str
    generated_at: datetime
    sections: List[EvidenceSection]


class SignatureValidationMetrics(BaseModel):
    signature_id: UUID
    num_matched_datasets: int
    num_total_datasets: int
    coverage: float
    specificity: Optional[float] = None
    reproducibility: Optional[float] = None
    mean_score: Optional[float] = None
    p_value: Optional[float] = None


class SignatureValidationResult(BaseModel):
    signature_id: UUID
    metrics: SignatureValidationMetrics
    matched_dataset_ids: List[UUID]
    summary: str
