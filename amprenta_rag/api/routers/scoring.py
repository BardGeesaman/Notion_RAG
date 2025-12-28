"""API router for relevance and novelty scoring."""
from __future__ import annotations

import time
from typing import List

from fastapi import APIRouter, HTTPException

from amprenta_rag.api.schemas import (
    BatchScoreRequest,
    BatchScoreResponse,
    NoveltyScoreRequest,
    NoveltyScoreResponse,
    RelevanceScoreRequest,
    RelevanceScoreResponse,
    ScoredItemResponse,
)

router = APIRouter()


@router.post("/relevance", response_model=RelevanceScoreResponse)
async def score_relevance_endpoint(
    request: RelevanceScoreRequest,
) -> RelevanceScoreResponse:
    """
    Score relevance of an item using LLM with structured output.
    
    Evaluates how relevant an item (dataset, paper, etc.) is to a given
    research context based on disease match, target overlap, and data quality.
    """
    try:
        from amprenta_rag.analysis.relevance_scoring import score_relevance
        
        # Convert request to internal format
        item = {
            "id": request.item.id,
            "title": request.item.title,
            "description": request.item.description,
            "species": request.item.species,
            "assay_type": request.item.assay_type,
            "sample_count": request.item.sample_count,
        }
        
        context = {
            "diseases": request.context.diseases or [],
            "targets": request.context.targets or [],
            "species": request.context.species or [],
            "assay_types": request.context.assay_types or [],
            "min_sample_size": request.context.min_sample_size,
        }
        
        # Score relevance
        result = score_relevance(item, context, request.criteria)
        
        return RelevanceScoreResponse(
            item_id=result.item_id,
            overall_score=result.overall_score,
            disease_match=result.disease_match,
            target_overlap=result.target_overlap,
            data_quality=result.data_quality,
            explanation=result.explanation,
            processing_time_seconds=result.processing_time_seconds,
            cached=result.cached,
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Relevance scoring failed: {str(e)}"
        )


@router.post("/novelty", response_model=NoveltyScoreResponse)
async def score_novelty_endpoint(
    request: NoveltyScoreRequest,
) -> NoveltyScoreResponse:
    """
    Score novelty using hybrid approach: embeddings for similarity + LLM for explanation.
    
    Computes how novel an item is compared to existing items using semantic
    similarity and provides LLM-generated explanation.
    """
    try:
        from amprenta_rag.analysis.relevance_scoring import score_novelty
        
        # Convert request to internal format
        item = {
            "id": request.item.id,
            "title": request.item.title,
            "description": request.item.description,
            "species": request.item.species,
            "assay_type": request.item.assay_type,
            "sample_count": request.item.sample_count,
        }
        
        existing_items = []
        for existing_item in request.existing_items:
            existing_items.append({
                "id": existing_item.id,
                "title": existing_item.title,
                "description": existing_item.description,
                "species": existing_item.species,
                "assay_type": existing_item.assay_type,
                "sample_count": existing_item.sample_count,
            })
        
        # Score novelty
        result = score_novelty(item, existing_items)
        
        return NoveltyScoreResponse(
            item_id=result.item_id,
            novelty_score=result.novelty_score,
            max_similarity=result.max_similarity,
            most_similar_item_id=result.most_similar_item_id,
            explanation=result.explanation,
            processing_time_seconds=result.processing_time_seconds,
            cached=result.cached,
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Novelty scoring failed: {str(e)}"
        )


@router.post("/batch", response_model=BatchScoreResponse)
async def batch_score_endpoint(
    request: BatchScoreRequest,
) -> BatchScoreResponse:
    """
    Score multiple items in batch for relevance and/or novelty.
    
    Efficiently processes multiple items with optional relevance and novelty scoring.
    Useful for ranking large sets of datasets or papers.
    """
    try:
        from amprenta_rag.analysis.relevance_scoring import batch_score
        
        start_time = time.time()
        
        # Convert request to internal format
        items = []
        for item_req in request.items:
            items.append({
                "id": item_req.id,
                "title": item_req.title,
                "description": item_req.description,
                "species": item_req.species,
                "assay_type": item_req.assay_type,
                "sample_count": item_req.sample_count,
            })
        
        context = {
            "diseases": request.context.diseases or [],
            "targets": request.context.targets or [],
            "species": request.context.species or [],
            "assay_types": request.context.assay_types or [],
            "min_sample_size": request.context.min_sample_size,
        }
        
        # Batch score
        results = batch_score(
            items, 
            context, 
            score_relevance_flag=request.score_relevance,
            score_novelty_flag=request.score_novelty
        )
        
        # Convert results to response format
        response_items = []
        for result in results:
            item_response = ScoredItemResponse(item_id=result.item_id)
            
            if result.relevance_score:
                item_response.relevance_score = RelevanceScoreResponse(
                    item_id=result.relevance_score.item_id,
                    overall_score=result.relevance_score.overall_score,
                    disease_match=result.relevance_score.disease_match,
                    target_overlap=result.relevance_score.target_overlap,
                    data_quality=result.relevance_score.data_quality,
                    explanation=result.relevance_score.explanation,
                    processing_time_seconds=result.relevance_score.processing_time_seconds,
                    cached=result.relevance_score.cached,
                )
            
            if result.novelty_score:
                item_response.novelty_score = NoveltyScoreResponse(
                    item_id=result.novelty_score.item_id,
                    novelty_score=result.novelty_score.novelty_score,
                    max_similarity=result.novelty_score.max_similarity,
                    most_similar_item_id=result.novelty_score.most_similar_item_id,
                    explanation=result.novelty_score.explanation,
                    processing_time_seconds=result.novelty_score.processing_time_seconds,
                    cached=result.novelty_score.cached,
                )
            
            response_items.append(item_response)
        
        processing_time = time.time() - start_time
        
        return BatchScoreResponse(
            items=response_items,
            total_items=len(response_items),
            processing_time_seconds=processing_time,
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Batch scoring failed: {str(e)}"
        )
