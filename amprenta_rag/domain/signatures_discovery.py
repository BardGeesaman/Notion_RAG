from typing import List, Optional, Set
from uuid import UUID

from pydantic import BaseModel


class DiscoveryDatasetSummary(BaseModel):
    dataset_id: UUID
    omics_type: str
    disease: Optional[str] = None
    matrix: Optional[str] = None
    features: Set[str]
    directions: Optional[dict[str, str]] = None  # feature → "up"/"down"


class Component(BaseModel):
    feature: str
    weight: Optional[float] = None
    direction: Optional[str] = None


class DiscoveredSignature(BaseModel):
    name: str
    modality: str
    components: List[Component]
    support: int  # count of datasets
    provenance: dict
