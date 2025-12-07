from typing import Any, Optional

from pydantic import BaseModel


class ChemistryCampaignMetadata(BaseModel):
    campaign_id: str
    name: str
    description: Optional[str] = None
    metadata: Optional[dict] = None


class ChemistryHitRow(BaseModel):
    hit_id: str
    compound_id: str
    score: Optional[float] = None
    metadata: Optional[dict] = None


class BiochemicalResultRow(BaseModel):
    result_id: str
    compound_id: str
    campaign_id: str
    value: Any
    metadata: Optional[dict] = None
