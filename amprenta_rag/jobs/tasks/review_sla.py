"""Celery tasks for SLA monitoring and review cycle processing."""

from datetime import datetime, timezone
from typing import Dict

from celery import shared_task

from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.services.review_sla import (
    get_overdue_reviews,
    send_reminder,
    escalate_review,
    check_sla_status,
)
from amprenta_rag.services.review_cycles import (
    get_due_cycles,
    create_cycle_reviews,
    advance_cycle,
)

logger = get_logger(__name__)


@shared_task(bind=True)
def check_overdue_reviews_task(self) -> Dict[str, int]:
    """
    Hourly task to check overdue reviews and send reminders/escalations.
    
    Returns:
        Dict with counts of checked, reminded, and escalated reviews
    """
    try:
        with db_session() as db:
            # Get all overdue reviews
            overdue_reviews = get_overdue_reviews(db)
            
            counts = {
                "checked": len(overdue_reviews),
                "reminded": 0,
                "escalated": 0,
                "errors": 0,
            }
            
            now = datetime.now(timezone.utc)
            
            for review in overdue_reviews:
                try:
                    # Check SLA status for severity
                    sla_status = check_sla_status(review)
                    hours_remaining = sla_status.get("hours_remaining", 0)
                    
                    # Severely overdue (>1.5x SLA time past due)
                    is_severely_overdue = hours_remaining < -review.sla.max_review_hours * 0.5 if review.sla else False
                    
                    # Check if reminder needed (last reminder > 24h ago or no reminder sent)
                    needs_reminder = (
                        review.reminder_sent_at is None or 
                        (now - review.reminder_sent_at).total_seconds() > 24 * 3600
                    )
                    
                    if is_severely_overdue and review.sla and review.sla.escalation_chain:
                        # Escalate severely overdue reviews
                        if escalate_review(review, db):
                            counts["escalated"] += 1
                            logger.info(f"Escalated severely overdue review {review.id}")
                        else:
                            counts["errors"] += 1
                    elif needs_reminder:
                        # Send reminder for overdue reviews
                        if send_reminder(review, db):
                            counts["reminded"] += 1
                            logger.info(f"Sent reminder for overdue review {review.id}")
                        else:
                            counts["errors"] += 1
                            
                except Exception as e:
                    logger.error(f"Error processing overdue review {review.id}: {e}")
                    counts["errors"] += 1
            
            logger.info(
                f"SLA check completed: {counts['checked']} checked, "
                f"{counts['reminded']} reminded, {counts['escalated']} escalated, "
                f"{counts['errors']} errors"
            )
            
            return counts
            
    except Exception as e:
        logger.error(f"Failed to run overdue reviews check: {e}")
        self.retry(exc=e, countdown=300, max_retries=3)  # Retry in 5 minutes
        return {"checked": 0, "reminded": 0, "escalated": 0, "errors": 1}


@shared_task(bind=True)
def process_due_cycles_task(self) -> Dict[str, int]:
    """
    Daily task to process due review cycles.
    
    Returns:
        Dict with counts of cycles processed and reviews created
    """
    try:
        with db_session() as db:
            # Get all due cycles
            due_cycles = get_due_cycles(db)
            
            counts = {
                "cycles_processed": 0,
                "reviews_created": 0,
                "errors": 0,
            }
            
            for cycle in due_cycles:
                try:
                    # Create reviews for the cycle
                    created_reviews = create_cycle_reviews(cycle, db)
                    counts["reviews_created"] += len(created_reviews)
                    
                    # Advance the cycle to next run
                    advance_cycle(cycle, db)
                    counts["cycles_processed"] += 1
                    
                    logger.info(
                        f"Processed cycle {cycle.id} ({cycle.name}): "
                        f"created {len(created_reviews)} reviews, "
                        f"next run at {cycle.next_run_at}"
                    )
                    
                except Exception as e:
                    logger.error(f"Error processing cycle {cycle.id}: {e}")
                    counts["errors"] += 1
            
            logger.info(
                f"Cycle processing completed: {counts['cycles_processed']} cycles, "
                f"{counts['reviews_created']} reviews created, "
                f"{counts['errors']} errors"
            )
            
            return counts
            
    except Exception as e:
        logger.error(f"Failed to run cycle processing: {e}")
        self.retry(exc=e, countdown=300, max_retries=3)  # Retry in 5 minutes
        return {"cycles_processed": 0, "reviews_created": 0, "errors": 1}
