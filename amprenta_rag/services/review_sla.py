"""SLA enforcement service for review management."""

from __future__ import annotations

from datetime import datetime, timezone, timedelta
from typing import Dict, List, Optional

from sqlalchemy.orm import Session

from amprenta_rag.database.models import ReviewSLA
from amprenta_rag.models.auth import EntityReview
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_default_sla(entity_type: str, db: Session) -> Optional[ReviewSLA]:
    """Find default SLA for entity type or global default.
    
    Args:
        entity_type: Entity type to find SLA for (dataset, experiment, etc.)
        db: Database session
        
    Returns:
        Default SLA for the entity type, or None if not found
    """
    try:
        # First try to find entity-specific default SLA
        sla = db.query(ReviewSLA).filter(
            ReviewSLA.entity_type == entity_type,
            ReviewSLA.is_default.is_(True),
            ReviewSLA.is_active.is_(True)
        ).first()
        
        if sla:
            logger.debug(f"Found entity-specific default SLA for {entity_type}: {sla.name}")
            return sla
        
        # Fall back to global default SLA (entity_type is None)
        sla = db.query(ReviewSLA).filter(
            ReviewSLA.entity_type.is_(None),
            ReviewSLA.is_default.is_(True),
            ReviewSLA.is_active.is_(True)
        ).first()
        
        if sla:
            logger.debug(f"Found global default SLA for {entity_type}: {sla.name}")
            return sla
        
        logger.warning(f"No default SLA found for entity type {entity_type}")
        return None
        
    except Exception as e:
        logger.error(f"Failed to get default SLA for {entity_type}: {e}")
        return None


def apply_sla(review: EntityReview, sla: ReviewSLA, db: Session) -> EntityReview:
    """Apply SLA rules to a review.
    
    Args:
        review: EntityReview to apply SLA to
        sla: ReviewSLA to apply
        db: Database session
        
    Returns:
        Updated EntityReview with SLA applied
    """
    try:
        # Calculate due date based on max_review_hours
        now = datetime.now(timezone.utc)
        due_at = now + timedelta(hours=sla.max_review_hours)
        
        # Apply SLA to review
        review.sla_id = sla.id
        review.due_at = due_at
        review.escalation_level = 0
        review.reminder_sent_at = None
        review.escalated_at = None
        
        db.commit()
        
        logger.info(f"Applied SLA {sla.name} to review {review.id}, due at {due_at}")
        return review
        
    except Exception as e:
        logger.error(f"Failed to apply SLA {sla.id} to review {review.id}: {e}")
        db.rollback()
        raise


def check_sla_status(review: EntityReview) -> Dict[str, any]:
    """Check SLA status for a review.
    
    Args:
        review: EntityReview to check status for
        
    Returns:
        Dict with status, hours_remaining, pct_elapsed
    """
    if not review.due_at or not review.sla:
        return {
            "status": "no_sla",
            "hours_remaining": None,
            "pct_elapsed": None
        }
    
    now = datetime.now(timezone.utc)
    
    # Calculate time metrics
    total_hours = review.sla.max_review_hours
    elapsed_hours = (now - review.reviewed_at).total_seconds() / 3600 if review.reviewed_at else 0
    hours_remaining = (review.due_at - now).total_seconds() / 3600
    pct_elapsed = min(100, (elapsed_hours / total_hours) * 100) if total_hours > 0 else 0
    
    # Determine status
    if hours_remaining <= 0:
        if review.status in ("pending", "in_review"):
            status = "breached"
        else:
            status = "overdue"
    elif pct_elapsed >= review.sla.warning_threshold_pct:
        status = "warning"
    else:
        status = "on_track"
    
    return {
        "status": status,
        "hours_remaining": hours_remaining,
        "pct_elapsed": pct_elapsed,
        "due_at": review.due_at,
        "warning_threshold": review.sla.warning_threshold_pct
    }


def get_overdue_reviews(db: Session) -> List[EntityReview]:
    """Find reviews that are overdue.
    
    Args:
        db: Database session
        
    Returns:
        List of overdue EntityReview objects
    """
    try:
        now = datetime.now(timezone.utc)
        
        overdue_reviews = db.query(EntityReview).filter(
            EntityReview.due_at < now,
            EntityReview.status.in_(["pending", "in_review"]),
            EntityReview.due_at.isnot(None)
        ).all()
        
        logger.info(f"Found {len(overdue_reviews)} overdue reviews")
        return overdue_reviews
        
    except Exception as e:
        logger.error(f"Failed to get overdue reviews: {e}")
        return []


def get_at_risk_reviews(db: Session) -> List[EntityReview]:
    """Find reviews that are at risk of breaching SLA.
    
    Args:
        db: Database session
        
    Returns:
        List of at-risk EntityReview objects
    """
    try:
        at_risk_reviews = []
        
        # Get active reviews with SLAs
        active_reviews = db.query(EntityReview).filter(
            EntityReview.status.in_(["pending", "in_review"]),
            EntityReview.sla_id.isnot(None),
            EntityReview.due_at.isnot(None)
        ).all()
        
        for review in active_reviews:
            status_info = check_sla_status(review)
            if status_info["status"] in ("warning", "overdue", "breached"):
                at_risk_reviews.append(review)
        
        logger.info(f"Found {len(at_risk_reviews)} at-risk reviews")
        return at_risk_reviews
        
    except Exception as e:
        logger.error(f"Failed to get at-risk reviews: {e}")
        return []


def send_reminder(review: EntityReview, db: Session) -> bool:
    """Send reminder notification for review.
    
    Args:
        review: EntityReview to send reminder for
        db: Database session
        
    Returns:
        True if reminder sent successfully
    """
    try:
        # TODO: Integrate with notification system when available
        # For now, just update the reminder timestamp
        
        review.reminder_sent_at = datetime.now(timezone.utc)
        db.commit()
        
        logger.info(f"Reminder sent for review {review.id} to reviewer {review.reviewer_id}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to send reminder for review {review.id}: {e}")
        db.rollback()
        return False


def escalate_review(review: EntityReview, db: Session) -> bool:
    """Escalate review to next level in escalation chain.
    
    Args:
        review: EntityReview to escalate
        db: Database session
        
    Returns:
        True if escalation successful
    """
    try:
        if not review.sla or not review.sla.escalation_chain:
            logger.warning(f"No escalation chain for review {review.id}")
            return False
        
        escalation_chain = review.sla.escalation_chain
        current_level = review.escalation_level
        
        if current_level >= len(escalation_chain):
            logger.warning(f"Review {review.id} already at max escalation level")
            return False
        
        # Get next escalation target
        next_reviewer_id = escalation_chain[current_level]
        
        # TODO: Integrate with notification system when available
        # For now, just update escalation tracking
        
        review.escalation_level = current_level + 1
        review.escalated_at = datetime.now(timezone.utc)
        # Optionally reassign reviewer
        # review.reviewer_id = UUID(next_reviewer_id)
        
        db.commit()
        
        logger.info(f"Escalated review {review.id} to level {review.escalation_level}, target: {next_reviewer_id}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to escalate review {review.id}: {e}")
        db.rollback()
        return False
