"""Compound ranking API endpoints."""

from __future__ import annotations

import asyncio
from typing import List
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from amprenta_rag.api.schemas import (
    CompoundRankingSchema,
    ObjectiveScoreSchema,
    ParetoPoint,
    ParetoRequest,
    ParetoResponse,
    RankingPreset,
    RankingRequest,
    RankingResponse,
    RankRequest,
    RerankRequest,
    RankingResult,
)
from amprenta_rag.database.base import get_db
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.ranking.scorer import MAX_COMPOUNDS, PRESETS
from amprenta_rag.ml.ranking.service import get_pareto_front, get_presets, rank_compounds

logger = get_logger(__name__)

router = APIRouter(prefix="/ranking", tags=["ranking"])


# Sync helper functions for LLM-based ranking
def _sync_rank_items(items: List[dict], criteria: str, context: dict):
    """Sync helper for LLM-based item ranking."""
    from amprenta_rag.analysis.ranking import rank_items
    return rank_items(items, criteria, context)


def _sync_rerank_items(items: List[dict], previous_ranking: List[str], criteria: str, context: dict):
    """Sync helper for LLM-based item re-ranking."""
    from amprenta_rag.analysis.ranking import rerank_items
    return rerank_items(items, previous_ranking, criteria, context)


@router.post("/score", response_model=RankingResponse)
async def score_compounds(
    request: RankingRequest,
    db: Session = Depends(get_db),
) -> RankingResponse:
    """Score and rank compounds using multi-objective optimization."""
    # Validate compound count
    if len(request.compound_ids) > MAX_COMPOUNDS:
        raise HTTPException(
            status_code=400,
            detail=f"Too many compounds ({len(request.compound_ids)}). Maximum allowed: {MAX_COMPOUNDS}"
        )
    
    if not request.compound_ids:
        raise HTTPException(status_code=400, detail="No compound IDs provided")
    
    # Convert string UUIDs to UUID objects
    try:
        compound_uuids = [UUID(id_str) for id_str in request.compound_ids]
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid UUID format in compound_ids")
    
    # Resolve weights from preset or use provided
    if request.weights:
        weights = request.weights
    elif request.preset and request.preset in PRESETS:
        weights = PRESETS[request.preset]
    else:
        # Default to balanced preset
        weights = PRESETS["balanced"]
    
    try:
        # Perform ranking
        rankings = rank_compounds(
            compound_ids=compound_uuids,
            weights=weights,
            db=db,
            include_pareto=request.include_pareto
        )
        
        # Convert to response schema
        ranking_schemas = []
        for ranking in rankings:
            objectives = [
                ObjectiveScoreSchema(
                    name=obj.name,
                    raw_value=obj.raw_value,
                    normalized=obj.normalized,
                    weight=obj.weight,
                    confidence=obj.confidence
                )
                for obj in ranking.objectives
            ]
            
            ranking_schemas.append(CompoundRankingSchema(
                compound_id=ranking.compound_id,
                smiles=ranking.smiles,
                objectives=objectives,
                weighted_score=ranking.weighted_score,
                pareto_rank=ranking.pareto_rank,
                rank=ranking.rank
            ))
        
        # Get Pareto front
        pareto_front_rankings = get_pareto_front(rankings)
        pareto_front_schemas = [
            schema for schema in ranking_schemas 
            if schema.compound_id in [r.compound_id for r in pareto_front_rankings]
        ]
        
        # Calculate skipped compounds
        skipped_compounds = len(request.compound_ids) - len(rankings)
        
        return RankingResponse(
            rankings=ranking_schemas,
            pareto_front=pareto_front_schemas,
            total_compounds=len(rankings),
            skipped_compounds=skipped_compounds
        )
        
    except Exception as e:
        logger.error("Compound ranking failed: %s", e)
        raise HTTPException(status_code=500, detail="Ranking computation failed")


@router.post("/pareto", response_model=ParetoResponse)
async def get_pareto_plot_data(
    request: ParetoRequest,
    db: Session = Depends(get_db),
) -> ParetoResponse:
    """Get data for Pareto plot visualization."""
    # Validate compound count
    if len(request.compound_ids) > MAX_COMPOUNDS:
        raise HTTPException(
            status_code=400,
            detail=f"Too many compounds ({len(request.compound_ids)}). Maximum allowed: {MAX_COMPOUNDS}"
        )
    
    # Convert string UUIDs to UUID objects
    try:
        compound_uuids = [UUID(id_str) for id_str in request.compound_ids]
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid UUID format in compound_ids")
    
    try:
        # Get rankings with balanced weights for Pareto analysis
        rankings = rank_compounds(
            compound_ids=compound_uuids,
            weights=PRESETS["balanced"],
            db=db,
            include_pareto=True
        )
        
        if not rankings:
            return ParetoResponse(
                points=[],
                x_label=request.x_objective,
                y_label=request.y_objective,
                frontier_ids=[]
            )
        
        # Build Pareto points
        points = []
        frontier_ids = []
        
        for ranking in rankings:
            # Extract objective values
            obj_dict = {obj.name: obj.normalized for obj in ranking.objectives}
            
            # Get x-axis value
            x_value = obj_dict.get(request.x_objective, 0.0)
            
            # Get y-axis value
            if request.y_objective == "liability_aggregate":
                # Aggregate liability score (inverted for plotting - higher = worse)
                herg_liability = 1.0 - obj_dict.get("herg", 0.0)
                alerts_liability = 1.0 - obj_dict.get("alerts", 0.0)
                logs_liability = 1.0 - obj_dict.get("logs", 0.0)
                y_value = (herg_liability + alerts_liability + logs_liability) / 3.0
            else:
                y_value = obj_dict.get(request.y_objective, 0.0)
            
            is_frontier = ranking.pareto_rank == 1
            if is_frontier:
                frontier_ids.append(ranking.compound_id)
            
            points.append(ParetoPoint(
                compound_id=ranking.compound_id,
                smiles=ranking.smiles,
                x_value=x_value,
                y_value=y_value,
                pareto_rank=ranking.pareto_rank,
                is_frontier=is_frontier
            ))
        
        # Generate labels
        x_label = request.x_objective.replace("_", " ").title()
        y_label = request.y_objective.replace("_", " ").title()
        if request.y_objective == "liability_aggregate":
            y_label = "Aggregate Liability Score"
        
        return ParetoResponse(
            points=points,
            x_label=x_label,
            y_label=y_label,
            frontier_ids=frontier_ids
        )
        
    except Exception as e:
        logger.error("Pareto plot data generation failed: %s", e)
        raise HTTPException(status_code=500, detail="Pareto plot data generation failed")


@router.get("/presets", response_model=List[RankingPreset])
async def get_ranking_presets() -> List[RankingPreset]:
    """Get available ranking weight presets."""
    try:
        presets = get_presets()
        return [
            RankingPreset(
                name=preset["name"],
                weights=preset["weights"],
                description=preset["description"]
            )
            for preset in presets
        ]
    except Exception as e:
        logger.error("Failed to get ranking presets: %s", e)
        raise HTTPException(status_code=500, detail="Failed to get ranking presets")


@router.post("/rank", response_model=RankingResult)
async def rank_items_endpoint(request: RankRequest) -> RankingResult:
    """
    Rank items using LLM-based analysis.
    
    Uses LLM to intelligently rank items based on provided criteria and context.
    Useful for ranking datasets, papers, or other research items by relevance,
    quality, or custom criteria.
    """
    try:
        # Convert request to internal format
        items = []
        for item in request.items:
            items.append({
                "id": item.id,
                "title": item.title,
                "description": item.description,
                "species": item.species,
                "assay_type": item.assay_type,
                "sample_count": item.sample_count,
            })
        
        context = {}
        if request.context:
            context = {
                "diseases": request.context.diseases or [],
                "targets": request.context.targets or [],
                "species": request.context.species or [],
                "assay_types": request.context.assay_types or [],
                "min_sample_size": request.context.min_sample_size,
            }
        
        # Rank items using async thread pool
        result = await asyncio.to_thread(
            _sync_rank_items,
            items,
            request.criteria or "",
            context
        )
        
        return RankingResult(
            ranked_items=result.ranked_items,
            criteria_used=result.criteria_used,
            explanation=result.explanation,
            processing_time_seconds=result.processing_time_seconds,
            cached=result.cached,
        )
        
    except Exception as e:
        logger.error("LLM-based ranking failed: %s", e)
        raise HTTPException(
            status_code=500,
            detail=f"LLM-based ranking failed: {str(e)}"
        )


@router.post("/rerank", response_model=RankingResult)
async def rerank_items_endpoint(request: RerankRequest) -> RankingResult:
    """
    Re-rank items using LLM-based analysis.
    
    Takes a previous ranking and re-ranks items based on new criteria or context.
    Useful for refining rankings with additional information or different priorities.
    """
    try:
        # Convert request to internal format
        items = []
        for item in request.items:
            items.append({
                "id": item.id,
                "title": item.title,
                "description": item.description,
                "species": item.species,
                "assay_type": item.assay_type,
                "sample_count": item.sample_count,
            })
        
        context = {}
        if request.context:
            context = {
                "diseases": request.context.diseases or [],
                "targets": request.context.targets or [],
                "species": request.context.species or [],
                "assay_types": request.context.assay_types or [],
                "min_sample_size": request.context.min_sample_size,
            }
        
        # Re-rank items using async thread pool
        result = await asyncio.to_thread(
            _sync_rerank_items,
            items,
            request.previous_ranking,
            request.criteria or "",
            context
        )
        
        return RankingResult(
            ranked_items=result.ranked_items,
            criteria_used=result.criteria_used,
            explanation=result.explanation,
            processing_time_seconds=result.processing_time_seconds,
            cached=result.cached,
        )
        
    except Exception as e:
        logger.error("LLM-based re-ranking failed: %s", e)
        raise HTTPException(
            status_code=500,
            detail=f"LLM-based re-ranking failed: {str(e)}"
        )
