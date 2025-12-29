"""
Experiment planner API endpoints.

Provides power analysis, effect size estimation, plate layout, and cost estimation.
"""

from __future__ import annotations

from typing import List

from fastapi import APIRouter
from pydantic import BaseModel, Field

from amprenta_rag.analysis.power_analysis import (
    calculate_plate_layout,
    calculate_sample_size,
    estimate_effect_size_from_data,
    estimate_experiment_cost,
)

router = APIRouter(prefix="/planner", tags=["Planner"])


class PowerRequest(BaseModel):
    """Request for power analysis."""

    effect_size: float
    alpha: float = 0.05
    power: float = 0.80
    test_type: str = "t-test"


class PowerResponse(BaseModel):
    """Response for power analysis."""

    n_per_group: int


class EffectSizeRequest(BaseModel):
    """Request for effect size estimation."""

    group1: List[float]
    group2: List[float]


class EffectSizeResponse(BaseModel):
    """Response for effect size estimation."""

    effect_size: float


class PlateLayoutRequest(BaseModel):
    """Request for plate layout calculation."""

    n: int = Field(..., ge=1)
    plate_format: int = Field(96, description="96, 384, or 1536")


class PlateLayoutResponse(BaseModel):
    """Response for plate layout."""

    plates_needed: int
    wells_used: int
    empty_wells: int


class CostRequest(BaseModel):
    """Request for cost estimation."""

    n: int = Field(..., ge=1)
    cost_per_sample: float = Field(..., gt=0)
    overhead_pct: float = Field(0.1, ge=0, le=1.0)


class CostResponse(BaseModel):
    """Response for cost estimation."""

    sample_cost: float
    overhead: float
    total: float


@router.post("/power", response_model=PowerResponse)
async def calculate_power_endpoint(request: PowerRequest) -> PowerResponse:
    """Calculate required sample size for desired power."""
    n = calculate_sample_size(
        effect_size=request.effect_size,
        alpha=request.alpha,
        power=request.power,
        test_type=request.test_type,  # type: ignore
    )
    return PowerResponse(n_per_group=n)


@router.post("/effect-size", response_model=EffectSizeResponse)
async def estimate_effect_size_endpoint(request: EffectSizeRequest) -> EffectSizeResponse:
    """Estimate effect size from pilot data."""
    effect_size = estimate_effect_size_from_data(request.group1, request.group2)
    return EffectSizeResponse(effect_size=effect_size)


@router.post("/plates", response_model=PlateLayoutResponse)
async def calculate_plates_endpoint(request: PlateLayoutRequest) -> PlateLayoutResponse:
    """Calculate plate layout for n samples."""
    layout = calculate_plate_layout(request.n, request.plate_format)
    return PlateLayoutResponse(**layout)


@router.post("/cost", response_model=CostResponse)
async def estimate_cost_endpoint(request: CostRequest) -> CostResponse:
    """Estimate experiment cost."""
    cost = estimate_experiment_cost(
        request.n,
        request.cost_per_sample,
        request.overhead_pct,
    )
    return CostResponse(**cost)

