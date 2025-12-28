from __future__ import annotations

import time
from typing import List
from uuid import UUID

from fastapi import APIRouter, HTTPException

from amprenta_rag.api import schemas
from amprenta_rag.api.services import screening as service

router = APIRouter()


@router.get(
    "/campaigns",
    summary="List HTS campaigns",
    response_model=List[schemas.CampaignResponse],
)
def list_campaigns():
    return service.list_campaigns()


@router.get(
    "/campaigns/{campaign_id}",
    summary="Get HTS campaign by ID",
    response_model=schemas.CampaignResponse,
)
def get_campaign(campaign_id: str):
    campaign = service.get_campaign(campaign_id)
    if not campaign:
        raise HTTPException(status_code=404, detail="Campaign not found")
    return campaign


@router.get(
    "/campaigns/{campaign_id}/hits",
    summary="Get HTS hits for campaign",
    response_model=List[schemas.HTSHitResponse],
)
def get_campaign_hits(campaign_id: str):
    return service.get_campaign_hits(campaign_id)


@router.post(
    "/suggest",
    summary="Get active learning suggestions for screening",
    response_model=schemas.ActiveLearningResponse,
)
async def suggest_compounds_for_screening(
    request: schemas.ActiveLearningRequest,
) -> schemas.ActiveLearningResponse:
    """
    Suggest next compounds to screen using active learning strategies.
    
    Uses acquisition functions to prioritize compounds based on uncertainty
    or diversity criteria. Supports both uncertainty sampling (requires model)
    and diversity-based selection using chemical fingerprints.
    """
    try:
        from amprenta_rag.analysis.active_learning import suggest_next_compounds
        
        start_time = time.time()
        
        # Convert request to internal format
        screened = [
            {
                "compound_id": s.compound_id,
                "smiles": s.smiles,
                "activity": s.activity,
            }
            for s in request.screened
        ]
        
        candidates = [
            {
                "compound_id": c.compound_id,
                "smiles": c.smiles,
            }
            for c in request.candidates
        ]
        
        # Get suggestions
        suggestions = suggest_next_compounds(
            screened=screened,
            candidates=candidates,
            strategy=request.strategy,
            batch_size=request.batch_size,
            model_id=request.model_id,
        )
        
        # Convert to response format
        suggestion_responses = [
            schemas.SuggestionResultResponse(
                compound_id=s.compound_id,
                smiles=s.smiles,
                acquisition_score=s.acquisition_score,
                strategy_used=s.strategy_used,
                rank=s.rank,
                explanation=s.explanation,
            )
            for s in suggestions
        ]
        
        processing_time = time.time() - start_time
        
        return schemas.ActiveLearningResponse(
            suggestions=suggestion_responses,
            strategy_used=request.strategy,
            total_candidates=len(candidates),
            total_screened=len(screened),
            batch_size=request.batch_size,
            processing_time_seconds=processing_time,
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Active learning suggestion failed: {str(e)}"
        )


@router.post(
    "/suggest/program/{program_id}",
    summary="Get active learning suggestions for a program",
    response_model=schemas.ActiveLearningResponse,
)
async def suggest_compounds_for_program(
    program_id: UUID,
    request: schemas.ProgramActiveLearningRequest,
) -> schemas.ActiveLearningResponse:
    """
    Suggest next compounds to screen for a specific program.
    
    Automatically fetches screened and candidate compounds from the database
    for the specified program and applies active learning strategies.
    """
    try:
        from amprenta_rag.analysis.active_learning import get_compound_suggestions_for_program
        
        start_time = time.time()
        
        # Get suggestions for the program
        suggestions = get_compound_suggestions_for_program(
            program_id=program_id,
            strategy=request.strategy,
            batch_size=request.batch_size,
            model_id=request.model_id,
        )
        
        # Convert to response format
        suggestion_responses = [
            schemas.SuggestionResultResponse(
                compound_id=s.compound_id,
                smiles=s.smiles,
                acquisition_score=s.acquisition_score,
                strategy_used=s.strategy_used,
                rank=s.rank,
                explanation=s.explanation,
            )
            for s in suggestions
        ]
        
        processing_time = time.time() - start_time
        
        return schemas.ActiveLearningResponse(
            suggestions=suggestion_responses,
            strategy_used=request.strategy,
            total_candidates=0,  # Not calculated in program-based approach
            total_screened=0,    # Not calculated in program-based approach
            batch_size=request.batch_size,
            processing_time_seconds=processing_time,
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Program-based active learning suggestion failed: {str(e)}"
        )

