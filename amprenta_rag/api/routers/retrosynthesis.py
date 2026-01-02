"""Retrosynthesis API endpoints."""

from __future__ import annotations

import logging
from typing import List
from uuid import uuid4

from fastapi import APIRouter, Depends, HTTPException, Query

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.api.schemas import (
    RetrosynthesisRequest,
    RetrosynthesisResponse,
    RouteScoreSchema,
    BuildingBlockResultSchema,
    SynthesisRouteSchema,
)
from amprenta_rag.chemistry.retrosynthesis import (
    RetrosynthesisPlanner,
    score_route,
    check_building_blocks,
)
from amprenta_rag.database.models import User

logger = logging.getLogger(__name__)

router = APIRouter()

# In-memory cache for MVP (replace with DB in Phase 2)
_analysis_cache: dict = {}


@router.post("/analyze", response_model=RetrosynthesisResponse)
def analyze_target(
    request: RetrosynthesisRequest,
    current_user: User = Depends(get_current_user),
):
    """Analyze target molecule and return synthesis routes."""
    logger.info(f"Retrosynthesis analysis requested for {request.smiles[:20]}...")
    
    planner = RetrosynthesisPlanner(backend="mock")
    tree = planner.analyze_target(request.smiles, request.max_depth)
    
    analysis_id = str(uuid4())
    _analysis_cache[analysis_id] = tree
    
    return RetrosynthesisResponse(analysis_id=analysis_id, tree=tree)


@router.get("/routes/{analysis_id}")
def get_routes(
    analysis_id: str,
    current_user: User = Depends(get_current_user),
):
    """Retrieve cached analysis by ID."""
    if analysis_id not in _analysis_cache:
        raise HTTPException(status_code=404, detail="Analysis not found")
    
    tree = _analysis_cache[analysis_id]
    return {"analysis_id": analysis_id, "tree": tree}


@router.post("/score", response_model=RouteScoreSchema)
def score_synthesis_route(
    route: SynthesisRouteSchema,
    current_user: User = Depends(get_current_user),
):
    """Score a synthesis route."""
    # Convert schema to dataclass
    from amprenta_rag.chemistry.retrosynthesis import SynthesisRoute, SynthesisStep
    
    steps = [SynthesisStep(
        reactants=s.reactants,
        product=s.product,
        reaction_type=s.reaction_type,
        conditions=s.conditions,
        confidence=s.confidence
    ) for s in route.steps]
    
    route_dc = SynthesisRoute(
        id=route.id,
        steps=steps,
        total_steps=route.total_steps,
        confidence=route.confidence
    )
    
    score = score_route(route_dc)
    return score


@router.get("/building-blocks", response_model=List[BuildingBlockResultSchema])
def get_building_blocks(
    smiles: List[str] = Query(..., description="SMILES strings to check"),
    current_user: User = Depends(get_current_user),
):
    """Check building block availability."""
    results = check_building_blocks(smiles)
    return results
