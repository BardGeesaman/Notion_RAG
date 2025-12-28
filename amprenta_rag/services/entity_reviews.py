"""Entity review workflow service with state machine."""

import logging
from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session
from sqlalchemy import and_

from amprenta_rag.models.auth import EntityReview, User
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import ActivityEventType
from amprenta_rag.services.activity import log_activity

logger = logging.getLogger(__name__)

# State machine constants
VALID_STATUSES = ["draft", "submitted", "in_review", "approved", "rejected", "changes_requested"]

VALID_TRANSITIONS = {
    "draft": ["submitted"],
    "submitted": ["in_review", "draft"],  # draft = cancel
    "in_review": ["approved", "rejected", "changes_requested"],
    "changes_requested": ["draft", "submitted"],
}


def create_review(
    entity_type: str,
    entity_id: UUID,
    submitted_by_id: UUID,
    db: Optional[Session] = None
) -> Optional[EntityReview]:
    """
    Create a new review for an entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity to review
        submitted_by_id: UUID of user creating the review
        db: Database session (optional)
    
    Returns:
        EntityReview object if successful, None if failed
    """
    if entity_type not in ["dataset", "experiment", "compound", "signature"]:
        logger.error(f"Invalid entity_type: {entity_type}")
        return None
    
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            # Check if there's already an active review for this entity
            existing_review = session.query(EntityReview).filter(
                and_(
                    EntityReview.entity_type == entity_type,
                    EntityReview.entity_id == entity_id,
                    EntityReview.status.in_(["draft", "submitted", "in_review", "changes_requested"])
                )
            ).first()
            
            if existing_review:
                logger.warning(
                    f"Active review already exists for {entity_type}:{entity_id} "
                    f"(status: {existing_review.status})"
                )
                return None
            
            # Create new review in draft state
            review = EntityReview(
                entity_type=entity_type,
                entity_id=entity_id,
                reviewer_id=submitted_by_id,  # Initially self-assigned
                status="draft",
                comments="Review created"
            )
            
            session.add(review)
            session.commit()
            session.expunge(review)
            
            logger.info(f"Review created: {review.id} for {entity_type}:{entity_id}")
            
            return review
            
    except Exception as e:
        logger.error(f"Failed to create review for {entity_type}:{entity_id}: {e}")
        return None


def submit_for_review(
    review_id: UUID,
    db: Optional[Session] = None
) -> Optional[EntityReview]:
    """
    Submit a review from draft to submitted status.
    
    Args:
        review_id: UUID of the review
        db: Database session (optional)
    
    Returns:
        Updated EntityReview object if successful, None if failed
    """
    review = _transition_review(
        review_id=review_id,
        new_status="submitted",
        comment="Review submitted for approval",
        db=db
    )
    
    if review:
        # Log activity event
        try:
            log_activity(
                event_type=ActivityEventType.REVIEW_SUBMITTED,
                target_type=review.entity_type,
                target_id=review.entity_id,
                target_name=f"{review.entity_type}:{review.entity_id}",
                actor_id=review.reviewer_id,  # The person who submitted
                program_id=None,  # Reviews are not program-specific
                metadata={
                    "entity_type": review.entity_type,
                    "entity_id": str(review.entity_id),
                    "review_id": str(review.id)
                }
            )
        except Exception as e:
            logger.error(f"Failed to log REVIEW_SUBMITTED activity: {e}")
    
    return review


def assign_reviewer(
    review_id: UUID,
    reviewer_id: UUID,
    db: Optional[Session] = None
) -> Optional[EntityReview]:
    """
    Assign a reviewer to a submitted review and move to in_review status.
    
    Args:
        review_id: UUID of the review
        reviewer_id: UUID of the reviewer to assign
        db: Database session (optional)
    
    Returns:
        Updated EntityReview object if successful, None if failed
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            review = session.query(EntityReview).filter(EntityReview.id == review_id).first()
            
            if not review:
                logger.error(f"Review not found: {review_id}")
                return None
            
            # Validate transition
            if not _is_valid_transition(review.status, "in_review"):
                logger.error(
                    f"Invalid transition from {review.status} to in_review "
                    f"for review {review_id}"
                )
                return None
            
            # Verify reviewer exists
            reviewer = session.query(User).filter(User.id == reviewer_id).first()
            if not reviewer:
                logger.error(f"Reviewer not found: {reviewer_id}")
                return None
            
            # Update review
            review.reviewer_id = reviewer_id
            review.status = "in_review"
            review.comments = f"Review assigned to {reviewer.username}"
            
            session.commit()
            session.expunge(review)
            
            logger.info(f"Review {review_id} assigned to {reviewer_id}")
            
            # Log activity event
            try:
                log_activity(
                    event_type=ActivityEventType.REVIEW_ASSIGNED,
                    target_type=review.entity_type,
                    target_id=review.entity_id,
                    target_name=f"{review.entity_type}:{review.entity_id}",
                    actor_id=reviewer_id,  # The assigned reviewer
                    program_id=None,  # Reviews are not program-specific
                    metadata={
                        "entity_type": review.entity_type,
                        "entity_id": str(review.entity_id),
                        "review_id": str(review.id),
                        "reviewer_id": str(reviewer_id)
                    }
                )
            except Exception as e:
                logger.error(f"Failed to log REVIEW_ASSIGNED activity: {e}")
            
            return review
            
    except Exception as e:
        logger.error(f"Failed to assign reviewer for review {review_id}: {e}")
        return None


def start_review(
    review_id: UUID,
    db: Optional[Session] = None
) -> Optional[EntityReview]:
    """
    Start a review (transition from submitted to in_review).
    
    Args:
        review_id: UUID of the review
        db: Database session (optional)
    
    Returns:
        Updated EntityReview object if successful, None if failed
    """
    return _transition_review(
        review_id=review_id,
        new_status="in_review",
        comment="Review started",
        db=db
    )


def decide_review(
    review_id: UUID,
    decision: str,
    comment: str,
    db: Optional[Session] = None
) -> Optional[EntityReview]:
    """
    Make a decision on a review (approve, reject, or request changes).
    
    Args:
        review_id: UUID of the review
        decision: Decision to make (approved, rejected, changes_requested)
        comment: Review comment/feedback
        db: Database session (optional)
    
    Returns:
        Updated EntityReview object if successful, None if failed
    """
    if decision not in ["approved", "rejected", "changes_requested"]:
        logger.error(f"Invalid decision: {decision}")
        return None
    
    review = _transition_review(
        review_id=review_id,
        new_status=decision,
        comment=comment,
        db=db
    )
    
    if review:
        # Log activity event
        try:
            log_activity(
                event_type=ActivityEventType.REVIEW_DECIDED,
                target_type=review.entity_type,
                target_id=review.entity_id,
                target_name=f"{review.entity_type}:{review.entity_id}",
                actor_id=review.reviewer_id,  # The reviewer who made the decision
                program_id=None,  # Reviews are not program-specific
                metadata={
                    "entity_type": review.entity_type,
                    "entity_id": str(review.entity_id),
                    "review_id": str(review.id),
                    "decision": decision
                }
            )
        except Exception as e:
            logger.error(f"Failed to log REVIEW_DECIDED activity: {e}")
    
    return review


def cancel_review(
    review_id: UUID,
    db: Optional[Session] = None
) -> Optional[EntityReview]:
    """
    Cancel a review (transition to draft status).
    
    Args:
        review_id: UUID of the review
        db: Database session (optional)
    
    Returns:
        Updated EntityReview object if successful, None if failed
    """
    return _transition_review(
        review_id=review_id,
        new_status="draft",
        comment="Review cancelled",
        db=db
    )


def list_reviews_for_entity(
    entity_type: str,
    entity_id: UUID,
    db: Optional[Session] = None
) -> List[EntityReview]:
    """
    List all reviews for a specific entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        db: Database session (optional)
    
    Returns:
        List of EntityReview objects
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            reviews = session.query(EntityReview).filter(
                and_(
                    EntityReview.entity_type == entity_type,
                    EntityReview.entity_id == entity_id
                )
            ).order_by(EntityReview.reviewed_at.desc()).all()
            
            # Detach from session
            for review in reviews:
                session.expunge(review)
            
            return reviews
            
    except Exception as e:
        logger.error(f"Failed to list reviews for {entity_type}:{entity_id}: {e}")
        return []


def get_pending_reviews(
    reviewer_id: UUID,
    db: Optional[Session] = None
) -> List[EntityReview]:
    """
    Get all pending reviews assigned to a reviewer.
    
    Args:
        reviewer_id: UUID of the reviewer
        db: Database session (optional)
    
    Returns:
        List of EntityReview objects with status in_review
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            reviews = session.query(EntityReview).filter(
                and_(
                    EntityReview.reviewer_id == reviewer_id,
                    EntityReview.status == "in_review"
                )
            ).order_by(EntityReview.reviewed_at.asc()).all()
            
            # Detach from session
            for review in reviews:
                session.expunge(review)
            
            return reviews
            
    except Exception as e:
        logger.error(f"Failed to get pending reviews for reviewer {reviewer_id}: {e}")
        return []


def get_my_reviews(
    user_id: UUID,
    status_filter: Optional[str] = None,
    db: Optional[Session] = None
) -> List[EntityReview]:
    """
    Get all reviews created by or assigned to a user.
    
    Args:
        user_id: UUID of the user
        status_filter: Optional status filter
        db: Database session (optional)
    
    Returns:
        List of EntityReview objects
    """
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            query = session.query(EntityReview).filter(
                EntityReview.reviewer_id == user_id
            )
            
            if status_filter and status_filter in VALID_STATUSES:
                query = query.filter(EntityReview.status == status_filter)
            
            reviews = query.order_by(EntityReview.reviewed_at.desc()).all()
            
            # Detach from session
            for review in reviews:
                session.expunge(review)
            
            return reviews
            
    except Exception as e:
        logger.error(f"Failed to get reviews for user {user_id}: {e}")
        return []


def _transition_review(
    review_id: UUID,
    new_status: str,
    comment: str,
    db: Optional[Session] = None
) -> Optional[EntityReview]:
    """
    Helper function to transition a review to a new status.
    
    Args:
        review_id: UUID of the review
        new_status: New status to transition to
        comment: Comment for the transition
        db: Database session (optional)
    
    Returns:
        Updated EntityReview object if successful, None if failed
    """
    if new_status not in VALID_STATUSES:
        logger.error(f"Invalid status: {new_status}")
        return None
    
    try:
        session_provided = db is not None
        
        with db if session_provided else db_session() as session:
            review = session.query(EntityReview).filter(EntityReview.id == review_id).first()
            
            if not review:
                logger.error(f"Review not found: {review_id}")
                return None
            
            # Validate transition
            if not _is_valid_transition(review.status, new_status):
                logger.error(
                    f"Invalid transition from {review.status} to {new_status} "
                    f"for review {review_id}"
                )
                return None
            
            # Update review
            old_status = review.status
            review.status = new_status
            review.comments = comment
            
            session.commit()
            session.expunge(review)
            
            logger.info(
                f"Review {review_id} transitioned from {old_status} to {new_status}"
            )
            
            return review
            
    except Exception as e:
        logger.error(f"Failed to transition review {review_id}: {e}")
        return None


def _is_valid_transition(current_status: str, new_status: str) -> bool:
    """
    Check if a status transition is valid.
    
    Args:
        current_status: Current review status
        new_status: Proposed new status
    
    Returns:
        True if transition is valid, False otherwise
    """
    if current_status not in VALID_TRANSITIONS:
        return False
    
    return new_status in VALID_TRANSITIONS[current_status]


def get_review_status_info() -> dict:
    """
    Get information about valid statuses and transitions.
    
    Returns:
        Dictionary with status and transition information
    """
    return {
        "valid_statuses": VALID_STATUSES,
        "valid_transitions": VALID_TRANSITIONS,
        "status_descriptions": {
            "draft": "Review is being prepared",
            "submitted": "Review submitted, waiting for reviewer assignment",
            "in_review": "Review is being conducted by assigned reviewer",
            "approved": "Review approved, entity can be published/used",
            "rejected": "Review rejected, entity needs major changes",
            "changes_requested": "Minor changes requested before approval"
        }
    }
