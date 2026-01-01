"""Review cycle management service for scheduled reviews."""

from __future__ import annotations

from datetime import datetime, timezone, timedelta
from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import ReviewCycle
from amprenta_rag.models.auth import EntityReview
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_cycle(
    name: str,
    entity_type: str,
    frequency: str,
    reviewer_pool: List[str],
    created_by_id: UUID,
    db: Session,
    program_id: Optional[UUID] = None,
    day_of_week: Optional[int] = None,
    day_of_month: Optional[int] = None,
    **kwargs
) -> ReviewCycle:
    """Create a new review cycle with computed next_run_at.
    
    Args:
        name: Human-readable cycle name
        entity_type: Entity type to review (dataset, experiment, etc.)
        frequency: Review frequency (weekly, monthly, quarterly, yearly)
        reviewer_pool: List of reviewer UUIDs
        created_by_id: Creator user ID
        db: Database session
        program_id: Optional program scoping
        day_of_week: Day of week for weekly cycles (0=Mon, 6=Sun)
        day_of_month: Day of month for monthly/quarterly cycles (1-28)
        
    Returns:
        Created ReviewCycle
    """
    try:
        # Compute initial next_run_at
        next_run_at = _compute_next_run(frequency, day_of_week, day_of_month)
        
        cycle = ReviewCycle(
            name=name,
            entity_type=entity_type,
            frequency=frequency,
            day_of_week=day_of_week,
            day_of_month=day_of_month,
            next_run_at=next_run_at,
            reviewer_pool=reviewer_pool,
            is_active=True,
            program_id=program_id,
            created_by_id=created_by_id,
            created_at=datetime.now(timezone.utc),
        )
        
        db.add(cycle)
        db.commit()
        db.refresh(cycle)
        
        logger.info(f"Created review cycle {cycle.id}: {name}, next run at {next_run_at}")
        return cycle
        
    except Exception as e:
        logger.error(f"Failed to create review cycle: {e}")
        db.rollback()
        raise


def get_cycles(
    db: Session,
    program_id: Optional[UUID] = None,
    entity_type: Optional[str] = None
) -> List[ReviewCycle]:
    """Get review cycles with optional filtering.
    
    Args:
        db: Database session
        program_id: Optional program filter
        entity_type: Optional entity type filter
        
    Returns:
        List of matching ReviewCycle objects
    """
    try:
        query = db.query(ReviewCycle)
        
        if program_id:
            query = query.filter(ReviewCycle.program_id == program_id)
        
        if entity_type:
            query = query.filter(ReviewCycle.entity_type == entity_type)
        
        cycles = query.filter(ReviewCycle.is_active.is_(True)).all()
        
        logger.debug(f"Found {len(cycles)} active cycles (program={program_id}, entity_type={entity_type})")
        return cycles
        
    except Exception as e:
        logger.error(f"Failed to get cycles: {e}")
        return []


def get_due_cycles(db: Session) -> List[ReviewCycle]:
    """Find cycles that are due to run.
    
    Args:
        db: Database session
        
    Returns:
        List of due ReviewCycle objects
    """
    try:
        now = datetime.now(timezone.utc)
        
        due_cycles = db.query(ReviewCycle).filter(
            ReviewCycle.next_run_at <= now,
            ReviewCycle.is_active.is_(True)
        ).all()
        
        logger.info(f"Found {len(due_cycles)} due cycles")
        return due_cycles
        
    except Exception as e:
        logger.error(f"Failed to get due cycles: {e}")
        return []


def create_cycle_reviews(cycle: ReviewCycle, db: Session) -> List[EntityReview]:
    """Create reviews for entities matching the cycle criteria.
    
    Args:
        cycle: ReviewCycle to create reviews for
        db: Database session
        
    Returns:
        List of created EntityReview objects
    """
    try:
        # NOTE: Simplified implementation - entity querying enhancement tracked in ROADMAP
        # Future: query actual entities (datasets, experiments, etc.) based on cycle.entity_type and cycle.program_id
        
        # For now, we'll simulate creating reviews for existing entities
        # that don't already have active reviews
        
        # Find entities without active reviews
        # This would need to be implemented based on your entity models
        
        created_reviews = []
        
        # Simulate finding entities that need review
        # In practice, this would query your Dataset, Experiment, etc. tables
        mock_entity_ids = [
            "550e8400-e29b-41d4-a716-446655440001",
            "550e8400-e29b-41d4-a716-446655440002",
            "550e8400-e29b-41d4-a716-446655440003",
        ]
        
        for i, entity_id in enumerate(mock_entity_ids):
            # Check if entity already has active review
            existing_review = db.query(EntityReview).filter(
                EntityReview.entity_type == cycle.entity_type,
                EntityReview.entity_id == UUID(entity_id),
                EntityReview.status.in_(["pending", "in_review"])
            ).first()
            
            if existing_review:
                logger.debug(f"Entity {entity_id} already has active review, skipping")
                continue
            
            # Assign reviewer from pool (rotate)
            reviewer_id = rotate_reviewer(cycle, i)
            
            # Create review
            review = EntityReview(
                entity_type=cycle.entity_type,
                entity_id=UUID(entity_id),
                reviewer_id=reviewer_id,
                status="pending",
                cycle_id=cycle.id,
                reviewed_at=datetime.now(timezone.utc),
                escalation_level=0,
            )
            
            db.add(review)
            created_reviews.append(review)
        
        db.commit()
        
        logger.info(f"Created {len(created_reviews)} reviews for cycle {cycle.id}")
        return created_reviews
        
    except Exception as e:
        logger.error(f"Failed to create cycle reviews for {cycle.id}: {e}")
        db.rollback()
        return []


def advance_cycle(cycle: ReviewCycle, db: Session) -> ReviewCycle:
    """Advance cycle to next scheduled run.
    
    Args:
        cycle: ReviewCycle to advance
        db: Database session
        
    Returns:
        Updated ReviewCycle
    """
    try:
        # Compute next run time
        next_run_at = _compute_next_run(
            cycle.frequency,
            cycle.day_of_week,
            cycle.day_of_month,
            from_date=cycle.next_run_at
        )
        
        cycle.next_run_at = next_run_at
        db.commit()
        
        logger.info(f"Advanced cycle {cycle.id} to next run at {next_run_at}")
        return cycle
        
    except Exception as e:
        logger.error(f"Failed to advance cycle {cycle.id}: {e}")
        db.rollback()
        raise


def rotate_reviewer(cycle: ReviewCycle, index: int = 0) -> UUID:
    """Get next reviewer from pool using round-robin.
    
    Args:
        cycle: ReviewCycle with reviewer pool
        index: Index for round-robin selection
        
    Returns:
        UUID of selected reviewer
    """
    if not cycle.reviewer_pool:
        raise ValueError(f"Cycle {cycle.id} has no reviewer pool")
    
    # Use modulo for round-robin selection
    reviewer_index = index % len(cycle.reviewer_pool)
    reviewer_id = cycle.reviewer_pool[reviewer_index]
    
    logger.debug(f"Selected reviewer {reviewer_id} (index {reviewer_index}) for cycle {cycle.id}")
    return UUID(reviewer_id)


def _compute_next_run(
    frequency: str,
    day_of_week: Optional[int] = None,
    day_of_month: Optional[int] = None,
    from_date: Optional[datetime] = None
) -> datetime:
    """Compute next run time based on frequency and constraints.
    
    Args:
        frequency: Review frequency (weekly, monthly, quarterly, yearly)
        day_of_week: Day of week for weekly cycles (0=Mon, 6=Sun)
        day_of_month: Day of month for monthly/quarterly cycles (1-28)
        from_date: Base date to compute from (default: now)
        
    Returns:
        Next scheduled run datetime
    """
    base_date = from_date or datetime.now(timezone.utc)
    
    if frequency == "weekly":
        # Find next occurrence of day_of_week
        if day_of_week is None:
            day_of_week = 1  # Default to Tuesday
        
        days_ahead = day_of_week - base_date.weekday()
        if days_ahead <= 0:  # Target day already happened this week
            days_ahead += 7
        
        next_run = base_date + timedelta(days=days_ahead)
        
    elif frequency == "monthly":
        # Find next occurrence of day_of_month
        if day_of_month is None:
            day_of_month = 1  # Default to 1st of month
        
        # Start with next month
        if base_date.month == 12:
            next_month = base_date.replace(year=base_date.year + 1, month=1, day=1)
        else:
            next_month = base_date.replace(month=base_date.month + 1, day=1)
        
        # Set to desired day of month
        try:
            next_run = next_month.replace(day=min(day_of_month, 28))  # Cap at 28 to avoid month issues
        except ValueError:
            next_run = next_month.replace(day=28)
        
    elif frequency == "quarterly":
        # Find next quarter
        if day_of_month is None:
            day_of_month = 1
        
        current_quarter = (base_date.month - 1) // 3
        next_quarter_month = ((current_quarter + 1) % 4) * 3 + 1
        
        if next_quarter_month <= base_date.month:
            # Next year
            next_run = base_date.replace(year=base_date.year + 1, month=next_quarter_month, day=min(day_of_month, 28))
        else:
            # Same year
            next_run = base_date.replace(month=next_quarter_month, day=min(day_of_month, 28))
        
    elif frequency == "yearly":
        # Next year, same month/day
        next_run = base_date.replace(year=base_date.year + 1)
        
    else:
        raise ValueError(f"Unsupported frequency: {frequency}")
    
    # Ensure time is set to business hours (9 AM UTC)
    next_run = next_run.replace(hour=9, minute=0, second=0, microsecond=0)
    
    return next_run
