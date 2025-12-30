"""
Experiment planner API endpoints.

Provides power analysis, effect size estimation, plate layout, and cost estimation.
Also provides LLM-based experiment planning, critique, refinement, and execution guidance.
"""

from __future__ import annotations

import asyncio
from typing import List

from fastapi import APIRouter
from pydantic import BaseModel, Field

from amprenta_rag.analysis.power_analysis import (
    calculate_plate_layout,
    calculate_sample_size,
    estimate_effect_size_from_data,
    estimate_experiment_cost,
)
from amprenta_rag.api.schemas import (
    PlanRequest,
    CritiqueRequest,
    RefineRequest,
    ExecuteRequest,
    PlanResult,
    CritiqueResult,
    RefinementResult,
    ExecutionGuidance,
)

router = APIRouter(prefix="/planner", tags=["Planner"])


# Sync helper functions for LLM-based planning
def _sync_create_plan(goal: str, context: dict, constraints: list):
    """Sync helper for LLM-based experiment planning."""
    from amprenta_rag.analysis.planner import create_plan
    return create_plan(goal, context, constraints)


def _sync_critique_plan(plan: str, criteria: list, context: dict):
    """Sync helper for LLM-based plan critique."""
    from amprenta_rag.analysis.planner import critique_plan
    return critique_plan(plan, criteria, context)


def _sync_refine_plan(original_plan: str, critique: str, additional_requirements: str):
    """Sync helper for LLM-based plan refinement."""
    from amprenta_rag.analysis.planner import refine_plan
    return refine_plan(original_plan, critique, additional_requirements)


def _sync_execute_plan(plan: str, current_step: int, issues: list):
    """Sync helper for LLM-based execution guidance."""
    from amprenta_rag.analysis.planner import execute_plan
    return execute_plan(plan, current_step, issues)


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


@router.post("/plan", response_model=PlanResult)
async def plan_endpoint(request: PlanRequest) -> PlanResult:
    """
    Create an experimental plan using LLM analysis.
    
    Uses LLM to generate a detailed experimental plan based on research goals,
    context, and constraints. Provides step-by-step guidance with resource
    estimates and risk assessment.
    """
    try:
        # Convert request to internal format
        context = {}
        if request.context:
            context = {
                "diseases": request.context.diseases or [],
                "targets": request.context.targets or [],
                "species": request.context.species or [],
                "assay_types": request.context.assay_types or [],
                "min_sample_size": request.context.min_sample_size,
            }
        
        constraints = request.constraints or []
        
        # Create plan using async thread pool
        result = await asyncio.to_thread(
            _sync_create_plan,
            request.goal,
            context,
            constraints
        )
        
        return result
        
    except Exception as e:
        from fastapi import HTTPException
        raise HTTPException(
            status_code=500,
            detail=f"LLM-based planning failed: {str(e)}"
        )


@router.post("/critique", response_model=CritiqueResult)
async def critique_endpoint(request: CritiqueRequest) -> CritiqueResult:
    """
    Critique an experimental plan using LLM analysis.
    
    Uses LLM to evaluate plan quality, identify strengths and weaknesses,
    and provide recommendations for improvement. Assesses feasibility
    and scientific rigor.
    """
    try:
        # Convert request to internal format
        context = {}
        if request.context:
            context = {
                "diseases": request.context.diseases or [],
                "targets": request.context.targets or [],
                "species": request.context.species or [],
                "assay_types": request.context.assay_types or [],
                "min_sample_size": request.context.min_sample_size,
            }
        
        criteria = request.criteria or []
        
        # Critique plan using async thread pool
        result = await asyncio.to_thread(
            _sync_critique_plan,
            request.plan,
            criteria,
            context
        )
        
        return result
        
    except Exception as e:
        from fastapi import HTTPException
        raise HTTPException(
            status_code=500,
            detail=f"LLM-based plan critique failed: {str(e)}"
        )


@router.post("/refine", response_model=RefinementResult)
async def refine_endpoint(request: RefineRequest) -> RefinementResult:
    """
    Refine an experimental plan using LLM analysis.
    
    Uses LLM to improve a plan based on critique feedback and additional
    requirements. Provides a refined plan with detailed rationale for
    changes made.
    """
    try:
        # Refine plan using async thread pool
        result = await asyncio.to_thread(
            _sync_refine_plan,
            request.original_plan,
            request.critique,
            request.additional_requirements or ""
        )
        
        return result
        
    except Exception as e:
        from fastapi import HTTPException
        raise HTTPException(
            status_code=500,
            detail=f"LLM-based plan refinement failed: {str(e)}"
        )


@router.post("/execute", response_model=ExecutionGuidance)
async def execute_endpoint(request: ExecuteRequest) -> ExecutionGuidance:
    """
    Get execution guidance for an experimental plan using LLM analysis.
    
    Uses LLM to provide step-by-step execution guidance, identify potential
    issues, and suggest troubleshooting approaches. Helps navigate complex
    experimental protocols.
    """
    try:
        # Get execution guidance using async thread pool
        result = await asyncio.to_thread(
            _sync_execute_plan,
            request.plan,
            request.current_step or 0,
            request.issues or []
        )
        
        return result
        
    except Exception as e:
        from fastapi import HTTPException
        raise HTTPException(
            status_code=500,
            detail=f"LLM-based execution guidance failed: {str(e)}"
        )

