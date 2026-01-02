"""Active Learning API endpoints."""
from uuid import UUID
from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException, Query, status
from datetime import datetime, timezone

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.api.schemas import (
    SampleSelectionRequest,
    LabelSubmitRequest,
    LabelQueueItemResponse,
    LabelQueueItemWithCompound,
    ActiveLearningCycleResponse,
    ActiveLearningStatsResponse,
    CreateCycleRequest,
)
from amprenta_rag.ml.active_learning import ActiveLearningService
from amprenta_rag.models.auth import User
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import LabelQueueItem, ActiveLearningCycle, Compound
import logging

logger = logging.getLogger(__name__)

router = APIRouter()


@router.post(
    "/models/{model_id}/select",
    response_model=List[LabelQueueItemResponse],
    summary="Select samples for labeling",
    description="Use uncertainty sampling to select informative compounds for labeling."
)
def select_samples(
    model_id: UUID,
    request: SampleSelectionRequest,
    current_user: User = Depends(get_current_user)
):
    """Select high-uncertainty samples for labeling.
    
    Uses various sampling strategies to identify the most informative compounds
    for manual labeling to improve model performance.
    """
    service = ActiveLearningService(model_id)
    try:
        items = service.select_samples(
            compound_ids=request.compound_ids,
            strategy=request.strategy,
            batch_size=request.batch_size
        )
        return items
    except ValueError as e:
        logger.error(f"Sample selection error for model {model_id}: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        logger.error(f"Unexpected error in sample selection: {e}")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Sample selection failed")


@router.get(
    "/models/{model_id}/queue",
    response_model=List[LabelQueueItemWithCompound],
    summary="Get labeling queue",
    description="Get pending items in the labeling queue for a model."
)
def get_queue(
    model_id: UUID,
    status: str = Query("pending", description="Filter by status"),
    limit: int = Query(50, le=200, description="Maximum items to return"),
    current_user: User = Depends(get_current_user)
):
    """Get items in the labeling queue.
    
    Returns queue items with compound details for the labeling interface.
    """
    service = ActiveLearningService(model_id)
    try:
        return service.get_labeling_queue(status=status, limit=limit)
    except Exception as e:
        logger.error(f"Error retrieving queue for model {model_id}: {e}")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Failed to retrieve queue")


@router.get(
    "/queue/{item_id}",
    response_model=LabelQueueItemWithCompound,
    summary="Get queue item details",
    description="Get detailed information about a specific queue item."
)
def get_queue_item(
    item_id: UUID,
    current_user: User = Depends(get_current_user)
):
    """Get a specific queue item with compound details.
    
    Returns detailed information about a queue item including
    compound structure and prediction metadata.
    """
    with db_session() as db:
        item = db.query(LabelQueueItem).filter(LabelQueueItem.id == item_id).first()
        if not item:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Queue item not found")
        
        compound = db.query(Compound).filter(Compound.id == item.compound_id).first()
        
        # Convert to response format
        result_data = {
            "id": item.id,
            "compound_id": item.compound_id,
            "model_id": item.model_id,
            "prediction": item.prediction,
            "uncertainty": item.uncertainty,
            "applicability_distance": item.applicability_distance,
            "selection_strategy": item.selection_strategy,
            "selection_batch": item.selection_batch,
            "priority_score": item.priority_score,
            "status": item.status,
            "assigned_to": item.assigned_to,
            "assigned_at": item.assigned_at,
            "label": item.label,
            "label_source": item.label_source,
            "label_confidence": item.label_confidence,
            "labeled_at": item.labeled_at,
            "labeled_by": item.labeled_by,
            "notes": item.notes,
            "created_at": item.created_at,
            "updated_at": item.updated_at,
            "compound_smiles": compound.smiles if compound else None,
            "compound_name": compound.compound_id if compound else None
        }
        
        return LabelQueueItemWithCompound(**result_data)


@router.post(
    "/queue/{item_id}/label",
    response_model=LabelQueueItemResponse,
    summary="Submit label",
    description="Submit a label for a queue item."
)
def submit_label(
    item_id: UUID,
    request: LabelSubmitRequest,
    current_user: User = Depends(get_current_user)
):
    """Submit a label for a queue item.
    
    Records the experimental or expert-estimated label value
    for a compound in the active learning queue.
    """
    with db_session() as db:
        item = db.query(LabelQueueItem).filter(LabelQueueItem.id == item_id).first()
        if not item:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Queue item not found")
        
        if item.status not in ["pending", "in_progress"]:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Item is not available for labeling")
        
        service = ActiveLearningService(item.model_id)
        try:
            return service.submit_label(
                item_id=item_id,
                label=request.label,
                source=request.source,
                confidence=request.confidence,
                user_id=current_user.id,
                notes=request.notes
            )
        except ValueError as e:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
        except Exception as e:
            logger.error(f"Error submitting label for item {item_id}: {e}")
            raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Failed to submit label")


@router.post(
    "/queue/{item_id}/skip",
    response_model=LabelQueueItemResponse,
    summary="Skip item",
    description="Mark a queue item as skipped."
)
def skip_item(
    item_id: UUID,
    current_user: User = Depends(get_current_user)
):
    """Skip a queue item.
    
    Marks a compound as skipped when it cannot be labeled
    (e.g., insufficient data, synthesis issues).
    """
    with db_session() as db:
        item = db.query(LabelQueueItem).filter(LabelQueueItem.id == item_id).first()
        if not item:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Queue item not found")
        
        if item.status not in ["pending", "in_progress"]:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Item cannot be skipped")
        
        service = ActiveLearningService(item.model_id)
        try:
            return service.skip_item(item_id=item_id, user_id=current_user.id)
        except ValueError as e:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
        except Exception as e:
            logger.error(f"Error skipping item {item_id}: {e}")
            raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Failed to skip item")


@router.post(
    "/queue/{item_id}/assign",
    response_model=LabelQueueItemResponse,
    summary="Assign item",
    description="Assign a queue item to a user."
)
def assign_item(
    item_id: UUID,
    assignee_id: UUID = Query(..., description="User ID to assign to"),
    current_user: User = Depends(get_current_user)
):
    """Assign a queue item to a user.
    
    Assigns a labeling task to a specific user and marks it as in progress.
    """
    with db_session() as db:
        item = db.query(LabelQueueItem).filter(LabelQueueItem.id == item_id).first()
        if not item:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Queue item not found")
        
        if item.status != "pending":
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Item is not available for assignment")
        
        item.assigned_to = assignee_id
        item.assigned_at = datetime.now(timezone.utc)
        item.status = "in_progress"
        db.commit()
        db.expunge(item)
        
        logger.info(f"Assigned queue item {item_id} to user {assignee_id}")
        return item


@router.get(
    "/models/{model_id}/cycles",
    response_model=List[ActiveLearningCycleResponse],
    summary="List AL cycles",
    description="List all active learning cycles for a model."
)
def list_cycles(
    model_id: UUID,
    current_user: User = Depends(get_current_user)
):
    """List active learning cycles for a model.
    
    Returns all cycles for a model in reverse chronological order
    with progress statistics and performance metrics.
    """
    with db_session() as db:
        cycles = db.query(ActiveLearningCycle).filter(
            ActiveLearningCycle.model_id == model_id
        ).order_by(ActiveLearningCycle.cycle_number.desc()).all()
        
        for c in cycles:
            db.expunge(c)
        return cycles


@router.post(
    "/models/{model_id}/cycles",
    response_model=ActiveLearningCycleResponse,
    summary="Create AL cycle",
    description="Create a new active learning cycle.",
    status_code=status.HTTP_201_CREATED
)
def create_cycle(
    model_id: UUID,
    request: CreateCycleRequest,
    current_user: User = Depends(get_current_user)
):
    """Create a new active learning cycle.
    
    Starts a new iteration of the active learning process
    with the specified strategy and batch size.
    """
    service = ActiveLearningService(model_id)
    try:
        return service.create_cycle(
            selection_strategy=request.selection_strategy,
            batch_size=request.batch_size
        )
    except ValueError as e:
        logger.error(f"Cycle creation error for model {model_id}: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        logger.error(f"Unexpected error creating cycle: {e}")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Failed to create cycle")


@router.get(
    "/models/{model_id}/stats",
    response_model=ActiveLearningStatsResponse,
    summary="Get AL statistics",
    description="Get active learning statistics for a model."
)
def get_stats(
    model_id: UUID,
    current_user: User = Depends(get_current_user)
):
    """Get active learning statistics for a model.
    
    Returns comprehensive statistics including cycle history,
    labeling progress, and performance metrics.
    """
    service = ActiveLearningService(model_id)
    try:
        stats = service.get_cycle_stats()
        return ActiveLearningStatsResponse(**stats)
    except Exception as e:
        logger.error(f"Error retrieving stats for model {model_id}: {e}")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Failed to retrieve statistics")
