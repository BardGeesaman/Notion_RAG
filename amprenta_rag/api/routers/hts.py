"""
HTS (high-throughput screening) API routes.

These endpoints expose read-only HTS quality control summaries and derived
analytics for a specific HTS campaign.
"""

from __future__ import annotations

from typing import List
from uuid import UUID

from fastapi import APIRouter, HTTPException

from amprenta_rag.analysis.hts_qc import (
    get_hit_compounds,
    get_plate_heatmap_data,
    get_plate_qc_summary,
)
from amprenta_rag.api import schemas
from amprenta_rag.database.models import HTSCampaign
from amprenta_rag.database.session import db_session

router = APIRouter()


def _campaign_exists(campaign_id: UUID) -> None:
    """Raise 404 if the requested HTS campaign does not exist."""
    with db_session() as db:
        exists = (
            db.query(HTSCampaign.id).filter(HTSCampaign.id == campaign_id).first() is not None
        )
    if not exists:
        raise HTTPException(status_code=404, detail="HTS campaign not found")


@router.get(
    "/hts/campaigns/{campaign_id}/qc",
    summary="QC summary for HTS campaign",
    response_model=schemas.PlateQCSummary,
)
def hts_qc_summary(campaign_id: UUID) -> schemas.PlateQCSummary:
    """Return a QC summary for a campaign (Z'-factor, controls, etc.)."""
    _campaign_exists(campaign_id)
    summary = get_plate_qc_summary(campaign_id)
    return schemas.PlateQCSummary(**summary.asdict())


@router.get(
    "/hts/campaigns/{campaign_id}/plate",
    summary="Plate heatmap data for HTS campaign",
    response_model=List[schemas.WellData],
)
def hts_plate_data(campaign_id: UUID) -> List[schemas.WellData]:
    """Return per-well heatmap data for a campaign plate view."""
    _campaign_exists(campaign_id)
    wells = get_plate_heatmap_data(campaign_id)
    return [schemas.WellData(**w.asdict()) for w in wells]


@router.get(
    "/hts/campaigns/{campaign_id}/hits",
    summary="Hit compounds for HTS campaign",
    response_model=List[schemas.HitCompound],
)
def hts_hits(campaign_id: UUID) -> List[schemas.HitCompound]:
    """Return hit compounds for a campaign based on QC/thresholding rules."""
    _campaign_exists(campaign_id)
    hits = get_hit_compounds(campaign_id)
    return [schemas.HitCompound(**h.asdict()) for h in hits]

