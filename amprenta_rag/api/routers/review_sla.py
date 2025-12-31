"""API endpoints for SLA management and review cycles."""

from datetime import datetime
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, status
from pydantic import BaseModel, ConfigDict
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.database.models import ReviewSLA, ReviewCycle
from amprenta_rag.models.auth import EntityReview, User
from amprenta_rag.services.review_sla import (
    check_sla_status,
    get_overdue_reviews,
)
from amprenta_rag.services.review_cycles import (
    create_cycle,
    get_cycles,
    create_cycle_reviews,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)
router = APIRouter(prefix="/sla", tags=["SLA"])


def require_admin_role(user: User):
    """Helper function to check admin role."""
    if user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Admin role required"
        )


# Pydantic Schemas
class SLARuleCreate(BaseModel):
    name: str
    entity_type: Optional[str] = None
    max_review_hours: int = 120
    warning_threshold_pct: int = 75
    escalation_chain: Optional[List[str]] = None
    is_default: bool = False
    is_active: bool = True


class SLARuleUpdate(BaseModel):
    name: Optional[str] = None
    entity_type: Optional[str] = None
    max_review_hours: Optional[int] = None
    warning_threshold_pct: Optional[int] = None
    escalation_chain: Optional[List[str]] = None
    is_default: Optional[bool] = None
    is_active: Optional[bool] = None


class SLARuleResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)
    
    id: UUID
    name: str
    entity_type: Optional[str]
    max_review_hours: int
    warning_threshold_pct: int
    escalation_chain: Optional[List[str]]
    is_default: bool
    is_active: bool
    created_at: datetime


class ReviewCycleCreate(BaseModel):
    name: str
    entity_type: str
    frequency: str
    reviewer_pool: List[str]
    program_id: Optional[UUID] = None
    day_of_week: Optional[int] = None
    day_of_month: Optional[int] = None


class ReviewCycleUpdate(BaseModel):
    name: Optional[str] = None
    entity_type: Optional[str] = None
    frequency: Optional[str] = None
    reviewer_pool: Optional[List[str]] = None
    program_id: Optional[UUID] = None
    day_of_week: Optional[int] = None
    day_of_month: Optional[int] = None
    is_active: Optional[bool] = None


class ReviewCycleResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)
    
    id: UUID
    name: str
    entity_type: str
    frequency: str
    reviewer_pool: List[str]
    program_id: Optional[UUID]
    day_of_week: Optional[int]
    day_of_month: Optional[int]
    next_run_at: Optional[datetime]
    is_active: bool
    created_at: datetime


class SLAStatusResponse(BaseModel):
    on_track: int
    warning: int
    overdue: int
    breached: int
    total: int


class ReviewSLAStatusResponse(BaseModel):
    review_id: UUID
    status: str
    hours_remaining: Optional[float]
    pct_elapsed: Optional[float]
    due_at: Optional[datetime]
    warning_threshold: Optional[int]
    sla_name: Optional[str]


# SLA Rules Endpoints (Admin Only)
@router.get("/rules", response_model=List[SLARuleResponse])
def list_sla_rules(
    entity_type: Optional[str] = Query(None),
    is_active: Optional[bool] = Query(None),
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """List all SLA rules with optional filtering."""
    require_admin_role(current_user)
    query = db.query(ReviewSLA)
    
    if entity_type:
        query = query.filter(ReviewSLA.entity_type == entity_type)
    
    if is_active is not None:
        query = query.filter(ReviewSLA.is_active == is_active)
    
    rules = query.order_by(ReviewSLA.created_at.desc()).all()
    
    logger.info(f"Admin {current_user.id} listed {len(rules)} SLA rules")
    return rules


@router.post("/rules", response_model=SLARuleResponse, status_code=status.HTTP_201_CREATED)
def create_sla_rule(
    rule_data: SLARuleCreate,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Create a new SLA rule."""
    require_admin_role(current_user)
    try:
        # Create SLA rule
        sla_rule = ReviewSLA(
            name=rule_data.name,
            entity_type=rule_data.entity_type,
            max_review_hours=rule_data.max_review_hours,
            warning_threshold_pct=rule_data.warning_threshold_pct,
            escalation_chain=rule_data.escalation_chain,
            is_default=rule_data.is_default,
            is_active=rule_data.is_active,
        )
        
        db.add(sla_rule)
        db.commit()
        db.refresh(sla_rule)
        
        logger.info(f"Admin {current_user.id} created SLA rule {sla_rule.id}: {sla_rule.name}")
        return sla_rule
        
    except Exception as e:
        logger.error(f"Failed to create SLA rule: {e}")
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to create SLA rule"
        )


@router.patch("/rules/{rule_id}", response_model=SLARuleResponse)
def update_sla_rule(
    rule_id: UUID,
    rule_data: SLARuleUpdate,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Update an existing SLA rule."""
    require_admin_role(current_user)
    rule = db.query(ReviewSLA).filter(ReviewSLA.id == rule_id).first()
    
    if not rule:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"SLA rule {rule_id} not found"
        )
    
    try:
        # Update fields
        update_data = rule_data.model_dump(exclude_unset=True)
        for field, value in update_data.items():
            setattr(rule, field, value)
        
        db.commit()
        db.refresh(rule)
        
        logger.info(f"Admin {current_user.id} updated SLA rule {rule_id}")
        return rule
        
    except Exception as e:
        logger.error(f"Failed to update SLA rule {rule_id}: {e}")
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to update SLA rule"
        )


@router.delete("/rules/{rule_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_sla_rule(
    rule_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Delete an SLA rule."""
    require_admin_role(current_user)
    rule = db.query(ReviewSLA).filter(ReviewSLA.id == rule_id).first()
    
    if not rule:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"SLA rule {rule_id} not found"
        )
    
    try:
        db.delete(rule)
        db.commit()
        
        logger.info(f"Admin {current_user.id} deleted SLA rule {rule_id}")
        
    except Exception as e:
        logger.error(f"Failed to delete SLA rule {rule_id}: {e}")
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to delete SLA rule"
        )


# Review Cycles Endpoints (Admin Only)
@router.get("/review-cycles", response_model=List[ReviewCycleResponse])
def list_review_cycles(
    program_id: Optional[UUID] = Query(None),
    entity_type: Optional[str] = Query(None),
    is_active: Optional[bool] = Query(None),
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """List review cycles with optional filtering."""
    require_admin_role(current_user)
    cycles = get_cycles(db, program_id=program_id, entity_type=entity_type)
    
    if is_active is not None:
        cycles = [cycle for cycle in cycles if cycle.is_active == is_active]
    
    logger.info(f"Admin {current_user.id} listed {len(cycles)} review cycles")
    return cycles


@router.post("/review-cycles", response_model=ReviewCycleResponse, status_code=status.HTTP_201_CREATED)
def create_review_cycle(
    cycle_data: ReviewCycleCreate,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Create a new review cycle."""
    require_admin_role(current_user)
    try:
        cycle = create_cycle(
            name=cycle_data.name,
            entity_type=cycle_data.entity_type,
            frequency=cycle_data.frequency,
            reviewer_pool=cycle_data.reviewer_pool,
            created_by_id=current_user.id,
            db=db,
            program_id=cycle_data.program_id,
            day_of_week=cycle_data.day_of_week,
            day_of_month=cycle_data.day_of_month,
        )
        
        logger.info(f"Admin {current_user.id} created review cycle {cycle.id}: {cycle.name}")
        return cycle
        
    except Exception as e:
        logger.error(f"Failed to create review cycle: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e)
        )


@router.patch("/review-cycles/{cycle_id}", response_model=ReviewCycleResponse)
def update_review_cycle(
    cycle_id: UUID,
    cycle_data: ReviewCycleUpdate,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Update an existing review cycle."""
    require_admin_role(current_user)
    cycle = db.query(ReviewCycle).filter(ReviewCycle.id == cycle_id).first()
    
    if not cycle:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Review cycle {cycle_id} not found"
        )
    
    try:
        # Update fields
        update_data = cycle_data.model_dump(exclude_unset=True)
        for field, value in update_data.items():
            setattr(cycle, field, value)
        
        db.commit()
        db.refresh(cycle)
        
        logger.info(f"Admin {current_user.id} updated review cycle {cycle_id}")
        return cycle
        
    except Exception as e:
        logger.error(f"Failed to update review cycle {cycle_id}: {e}")
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to update review cycle"
        )


@router.delete("/review-cycles/{cycle_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_review_cycle(
    cycle_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Delete a review cycle."""
    require_admin_role(current_user)
    cycle = db.query(ReviewCycle).filter(ReviewCycle.id == cycle_id).first()
    
    if not cycle:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Review cycle {cycle_id} not found"
        )
    
    try:
        db.delete(cycle)
        db.commit()
        
        logger.info(f"Admin {current_user.id} deleted review cycle {cycle_id}")
        
    except Exception as e:
        logger.error(f"Failed to delete review cycle {cycle_id}: {e}")
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to delete review cycle"
        )


@router.post("/review-cycles/{cycle_id}/run-now", response_model=dict)
def run_cycle_now(
    cycle_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Manually trigger a review cycle."""
    require_admin_role(current_user)
    cycle = db.query(ReviewCycle).filter(ReviewCycle.id == cycle_id).first()
    
    if not cycle:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Review cycle {cycle_id} not found"
        )
    
    try:
        created_reviews = create_cycle_reviews(cycle, db)
        
        logger.info(f"Admin {current_user.id} manually triggered cycle {cycle_id}")
        return {
            "cycle_id": cycle_id,
            "reviews_created": len(created_reviews),
            "message": f"Created {len(created_reviews)} reviews for cycle {cycle.name}"
        }
        
    except Exception as e:
        logger.error(f"Failed to run cycle {cycle_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e)
        )


# SLA Status Endpoints (Authenticated Users)
@router.get("/status", response_model=SLAStatusResponse)
def get_sla_status_summary(
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Get SLA status summary dashboard."""
    try:
        # Get all active reviews with SLAs
        active_reviews = db.query(EntityReview).filter(
            EntityReview.status.in_(["pending", "in_review"]),
            EntityReview.sla_id.isnot(None),
        ).all()
        
        counts = {"on_track": 0, "warning": 0, "overdue": 0, "breached": 0}
        
        for review in active_reviews:
            status_info = check_sla_status(review)
            status_key = status_info["status"]
            
            if status_key in counts:
                counts[status_key] += 1
        
        counts["total"] = sum(counts.values())
        
        logger.info(f"User {current_user.id} requested SLA status summary")
        return SLAStatusResponse(**counts)
        
    except Exception as e:
        logger.error(f"Failed to get SLA status summary: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get SLA status"
        )


@router.get("/overdue", response_model=List[dict])
def get_overdue_reviews_list(
    limit: int = Query(50, le=200),
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Get list of overdue reviews."""
    try:
        overdue_reviews = get_overdue_reviews(db)
        
        # Limit results and add SLA status
        limited_reviews = overdue_reviews[:limit]
        
        result = []
        for review in limited_reviews:
            status_info = check_sla_status(review)
            result.append({
                "review_id": review.id,
                "entity_type": review.entity_type,
                "entity_id": review.entity_id,
                "reviewer_id": review.reviewer_id,
                "status": review.status,
                "sla_status": status_info["status"],
                "hours_remaining": status_info["hours_remaining"],
                "due_at": review.due_at,
            })
        
        logger.info(f"User {current_user.id} requested {len(result)} overdue reviews")
        return result
        
    except Exception as e:
        logger.error(f"Failed to get overdue reviews: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get overdue reviews"
        )


@router.get("/reviews/{review_id}/sla", response_model=ReviewSLAStatusResponse)
def get_review_sla_status(
    review_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Get SLA status for a specific review."""
    review = db.query(EntityReview).filter(EntityReview.id == review_id).first()
    
    if not review:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Review {review_id} not found"
        )
    
    try:
        status_info = check_sla_status(review)
        
        return ReviewSLAStatusResponse(
            review_id=review.id,
            status=status_info["status"],
            hours_remaining=status_info["hours_remaining"],
            pct_elapsed=status_info["pct_elapsed"],
            due_at=status_info.get("due_at"),
            warning_threshold=status_info.get("warning_threshold"),
            sla_name=review.sla.name if review.sla else None,
        )
        
    except Exception as e:
        logger.error(f"Failed to get SLA status for review {review_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get review SLA status"
        )
