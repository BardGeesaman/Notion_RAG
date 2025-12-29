"""
Compound portfolio API endpoints.

Provides portfolio summary, ADMET rollup, SAR gaps, and recommendations.
"""

from __future__ import annotations

from typing import List

from fastapi import APIRouter, Depends
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.services.portfolio_service import (
    get_admet_summary,
    get_portfolio_summary,
    get_recommendations,
    get_sar_gaps,
)

router = APIRouter(prefix="/portfolio", tags=["Portfolio"])


class PortfolioSummaryResponse(BaseModel):
    """Portfolio summary response."""

    total_compounds: int
    scaffold_count: int
    date_from: str | None
    date_to: str | None


class ADMETSummaryResponse(BaseModel):
    """ADMET traffic light summary response."""

    green: int
    yellow: int
    red: int
    unknown: int


class SARGapResponse(BaseModel):
    """SAR gap response."""

    scaffold_id: str
    compound_count: int
    example_smiles: str | None


class RecommendationResponse(BaseModel):
    """Compound recommendation response."""

    compound_id: str
    smiles: str
    score: float
    reason: str


@router.get("/summary", response_model=PortfolioSummaryResponse)
async def get_portfolio_summary_endpoint(
    db: Session = Depends(get_database_session),
) -> PortfolioSummaryResponse:
    """Get compound portfolio summary."""
    summary = get_portfolio_summary(db)
    return PortfolioSummaryResponse(**summary)


@router.get("/admet", response_model=ADMETSummaryResponse)
async def get_admet_summary_endpoint(
    db: Session = Depends(get_database_session),
) -> ADMETSummaryResponse:
    """Get ADMET traffic light summary."""
    summary = get_admet_summary(db)
    return ADMETSummaryResponse(**summary)


@router.get("/gaps", response_model=List[SARGapResponse])
async def get_sar_gaps_endpoint(
    db: Session = Depends(get_database_session),
) -> List[SARGapResponse]:
    """Get SAR gaps (scaffolds needing expansion)."""
    gaps = get_sar_gaps(db)
    return [SARGapResponse(**gap) for gap in gaps]


@router.get("/recommendations", response_model=List[RecommendationResponse])
async def get_recommendations_endpoint(
    limit: int = 5,
    db: Session = Depends(get_database_session),
) -> List[RecommendationResponse]:
    """Get compound recommendations."""
    recommendations = get_recommendations(db, limit=limit)
    return [RecommendationResponse(**rec) for rec in recommendations]

