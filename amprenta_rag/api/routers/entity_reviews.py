"""API router for entity review workflow functionality."""

from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.api.schemas import (
    EntityReviewResponse,
    EntityReviewList,
    EntityReviewDecision,
    EntityReviewAssign,
    EntityReviewStatusInfo,
)
from amprenta_rag.models.auth import User
from amprenta_rag.services.entity_reviews import (
    create_review,
    submit_for_review,
    assign_reviewer,
    decide_review,
    cancel_review,
    list_reviews_for_entity,
    get_pending_reviews,
    get_my_reviews,
    get_review_status_info,
)

router = APIRouter(prefix="/reviews", tags=["entity-reviews"])


@router.post("/entities/{entity_type}/{entity_id}/reviews", status_code=status.HTTP_201_CREATED)
def create_entity_review(
    entity_type: str,
    entity_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityReviewResponse:
    """
    Create a new review for an entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity to review
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        EntityReview details
    """
    # Validate entity type
    valid_types = ["dataset", "experiment", "compound", "signature"]
    if entity_type not in valid_types:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type. Must be one of: {', '.join(valid_types)}"
        )
    
    # TODO: Add authorization check - user must have permission to create reviews for this entity
    
    review = create_review(
        entity_type=entity_type,
        entity_id=entity_id,
        submitted_by_id=current_user.id,
        db=db,
    )
    
    if not review:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Failed to create review. An active review may already exist for this entity."
        )
    
    return EntityReviewResponse(
        id=review.id,
        entity_type=review.entity_type,
        entity_id=review.entity_id,
        reviewer_id=review.reviewer_id,
        status=review.status,
        comments=review.comments,
        reviewed_at=review.reviewed_at,
    )


@router.get("/entities/{entity_type}/{entity_id}/reviews")
def list_entity_reviews(
    entity_type: str,
    entity_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityReviewList:
    """
    List all reviews for a specific entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        List of reviews for the entity
    """
    # Validate entity type
    valid_types = ["dataset", "experiment", "compound", "signature"]
    if entity_type not in valid_types:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type. Must be one of: {', '.join(valid_types)}"
        )
    
    # TODO: Add authorization check - user must have permission to view reviews for this entity
    
    reviews = list_reviews_for_entity(
        entity_type=entity_type,
        entity_id=entity_id,
        db=db,
    )
    
    review_responses = [
        EntityReviewResponse(
            id=review.id,
            entity_type=review.entity_type,
            entity_id=review.entity_id,
            reviewer_id=review.reviewer_id,
            status=review.status,
            comments=review.comments,
            reviewed_at=review.reviewed_at,
        )
        for review in reviews
    ]
    
    return EntityReviewList(reviews=review_responses)


@router.post("/reviews/{review_id}/submit")
def submit_review(
    review_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityReviewResponse:
    """
    Submit a review for approval (transition from draft to submitted).
    
    Args:
        review_id: UUID of the review
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        Updated review details
    """
    # TODO: Add authorization check - user must own the review or have admin permissions
    
    review = submit_for_review(review_id=review_id, db=db)
    
    if not review:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Failed to submit review. Check review status and permissions."
        )
    
    return EntityReviewResponse(
        id=review.id,
        entity_type=review.entity_type,
        entity_id=review.entity_id,
        reviewer_id=review.reviewer_id,
        status=review.status,
        comments=review.comments,
        reviewed_at=review.reviewed_at,
    )


@router.post("/reviews/{review_id}/assign")
def assign_review(
    review_id: UUID,
    assignment: EntityReviewAssign,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityReviewResponse:
    """
    Assign a reviewer to a review (transition to in_review status).
    
    Args:
        review_id: UUID of the review
        assignment: Reviewer assignment details
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        Updated review details
    """
    # TODO: Add authorization check - user must have admin permissions or own the review
    
    review = assign_reviewer(
        review_id=review_id,
        reviewer_id=assignment.reviewer_id,
        db=db,
    )
    
    if not review:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Failed to assign reviewer. Check review status, reviewer exists, and permissions."
        )
    
    return EntityReviewResponse(
        id=review.id,
        entity_type=review.entity_type,
        entity_id=review.entity_id,
        reviewer_id=review.reviewer_id,
        status=review.status,
        comments=review.comments,
        reviewed_at=review.reviewed_at,
    )


@router.post("/reviews/{review_id}/decide")
def decide_review_endpoint(
    review_id: UUID,
    decision: EntityReviewDecision,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityReviewResponse:
    """
    Make a decision on a review (approve, reject, or request changes).
    
    Args:
        review_id: UUID of the review
        decision: Review decision details
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        Updated review details
    """
    # TODO: Add authorization check - user must be the assigned reviewer or admin
    
    review = decide_review(
        review_id=review_id,
        decision=decision.decision,
        comment=decision.comment,
        db=db,
    )
    
    if not review:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Failed to make review decision. Check review status and permissions."
        )
    
    return EntityReviewResponse(
        id=review.id,
        entity_type=review.entity_type,
        entity_id=review.entity_id,
        reviewer_id=review.reviewer_id,
        status=review.status,
        comments=review.comments,
        reviewed_at=review.reviewed_at,
    )


@router.delete("/reviews/{review_id}")
def cancel_review_endpoint(
    review_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> dict:
    """
    Cancel a review (transition to draft status).
    
    Args:
        review_id: UUID of the review
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        Success message
    """
    # TODO: Add authorization check - user must own the review or have admin permissions
    
    review = cancel_review(review_id=review_id, db=db)
    
    if not review:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Failed to cancel review. Check review status and permissions."
        )
    
    return {"message": "Review cancelled successfully"}


@router.get("/my-reviews")
def get_my_reviews_endpoint(
    status: Optional[str] = None,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityReviewList:
    """
    Get all reviews created by or assigned to the current user.
    
    Args:
        status: Optional status filter (query param)
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        List of user's reviews
    """
    if status:
        valid_statuses = ["draft", "submitted", "in_review", "approved", "rejected", "changes_requested"]
        if status not in valid_statuses:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Invalid status. Must be one of: {', '.join(valid_statuses)}"
            )
    
    reviews = get_my_reviews(
        user_id=current_user.id,
        status_filter=status,
        db=db,
    )
    
    review_responses = [
        EntityReviewResponse(
            id=review.id,
            entity_type=review.entity_type,
            entity_id=review.entity_id,
            reviewer_id=review.reviewer_id,
            status=review.status,
            comments=review.comments,
            reviewed_at=review.reviewed_at,
        )
        for review in reviews
    ]
    
    return EntityReviewList(reviews=review_responses)


@router.get("/pending")
def get_pending_reviews_endpoint(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityReviewList:
    """
    Get all pending reviews assigned to the current user.
    
    Args:
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        List of pending reviews
    """
    reviews = get_pending_reviews(
        reviewer_id=current_user.id,
        db=db,
    )
    
    review_responses = [
        EntityReviewResponse(
            id=review.id,
            entity_type=review.entity_type,
            entity_id=review.entity_id,
            reviewer_id=review.reviewer_id,
            status=review.status,
            comments=review.comments,
            reviewed_at=review.reviewed_at,
        )
        for review in reviews
    ]
    
    return EntityReviewList(reviews=review_responses)


@router.get("/status-info")
def get_status_info() -> EntityReviewStatusInfo:
    """
    Get information about valid review statuses and transitions.
    
    Returns:
        Status and transition information
    """
    info = get_review_status_info()
    
    return EntityReviewStatusInfo(
        valid_statuses=info["valid_statuses"],
        valid_transitions=info["valid_transitions"],
        status_descriptions=info["status_descriptions"],
    )
