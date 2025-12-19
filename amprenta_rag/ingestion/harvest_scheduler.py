"""Harvest scheduler service for automated repository scanning."""
from __future__ import annotations

import logging
from datetime import datetime, timedelta
from typing import Optional

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import HarvestSchedule
from amprenta_rag.ingestion.discovery_service import run_discovery_job

logger = logging.getLogger(__name__)

try:
    from apscheduler.schedulers.background import BackgroundScheduler
    from apscheduler.triggers.interval import IntervalTrigger
    APSCHEDULER_AVAILABLE = True
except ImportError:
    APSCHEDULER_AVAILABLE = False
    BackgroundScheduler = None  # type: ignore
    IntervalTrigger = None  # type: ignore


_scheduler_instance: Optional[BackgroundScheduler] = None


def get_scheduler() -> BackgroundScheduler:
    """Get or create the scheduler singleton."""
    global _scheduler_instance
    if not APSCHEDULER_AVAILABLE:
        raise ImportError("APScheduler is not installed. Install with: pip install apscheduler")

    if _scheduler_instance is None:
        _scheduler_instance = BackgroundScheduler()
        logger.info("[HARVEST] Created scheduler instance")

    return _scheduler_instance


def start_scheduler() -> None:
    """Start the scheduler if not already running."""
    scheduler = get_scheduler()
    if not scheduler.running:
        scheduler.start()
        logger.info("[HARVEST] Scheduler started")
    else:
        logger.debug("[HARVEST] Scheduler already running")


def stop_scheduler() -> None:
    """Shutdown the scheduler."""
    global _scheduler_instance
    if _scheduler_instance is not None and _scheduler_instance.running:
        _scheduler_instance.shutdown()
        logger.info("[HARVEST] Scheduler stopped")
    _scheduler_instance = None


def run_harvest(schedule_id: str) -> None:
    """
    Run a harvest job for a specific schedule.

    Args:
        schedule_id: UUID of the harvest schedule
    """
    db_gen = get_db()
    db = next(db_gen)
    try:
        schedule = db.query(HarvestSchedule).filter(HarvestSchedule.id == schedule_id).first()
        if not schedule:
            logger.warning("[HARVEST] Schedule not found: %s", schedule_id)
            return

        if not schedule.is_active:
            logger.debug("[HARVEST] Schedule is inactive, skipping: %s", schedule.name)
            return

        logger.info("[HARVEST] Running harvest: %s (%s)", schedule.name, schedule.repository)

        # Update last_run
        schedule.last_run = datetime.utcnow()

        # Run discovery job
        try:
            run_discovery_job(
                repository=schedule.repository,
                query=schedule.query,
                max_results=50,
                user_id=str(schedule.created_by_id) if schedule.created_by_id else None,
            )

            # Calculate next_run
            schedule.next_run = datetime.utcnow() + timedelta(hours=schedule.interval_hours)
            db.commit()
            logger.info("[HARVEST] Harvest completed, next run: %s", schedule.next_run)
        except Exception as e:
            logger.error("[HARVEST] Harvest failed for %s: %r", schedule.name, e)
            db.rollback()
    finally:
        db_gen.close()


def schedule_harvest(schedule_id: str) -> None:
    """
    Add a harvest schedule to the scheduler.

    Args:
        schedule_id: UUID of the harvest schedule
    """
    db_gen = get_db()
    db = next(db_gen)
    try:
        schedule = db.query(HarvestSchedule).filter(HarvestSchedule.id == schedule_id).first()
        if not schedule:
            logger.warning("[HARVEST] Schedule not found: %s", schedule_id)
            return

        scheduler = get_scheduler()
        job_id = f"harvest_{schedule_id}"

        # Remove existing job if present
        if scheduler.get_job(job_id):
            scheduler.remove_job(job_id)

        # Calculate next run time
        if schedule.next_run:
            start_date = schedule.next_run
        elif schedule.last_run:
            start_date = schedule.last_run + timedelta(hours=schedule.interval_hours)
        else:
            start_date = datetime.utcnow()

        # Add job
        scheduler.add_job(
            run_harvest,
            trigger=IntervalTrigger(hours=schedule.interval_hours),
            args=[str(schedule.id)],
            id=job_id,
            name=f"Harvest: {schedule.name}",
            start_date=start_date,
            replace_existing=True,
        )
        logger.info("[HARVEST] Scheduled harvest: %s (every %d hours)", schedule.name, schedule.interval_hours)
    finally:
        db_gen.close()


def load_active_schedules() -> None:
    """Load all active harvest schedules from DB and add them to the scheduler."""
    db_gen = get_db()
    db = next(db_gen)
    try:
        schedules = db.query(HarvestSchedule).filter(HarvestSchedule.is_active).all()
        logger.info("[HARVEST] Loading %d active schedules", len(schedules))

        for schedule in schedules:
            try:
                schedule_harvest(str(schedule.id))
            except Exception as e:
                logger.error("[HARVEST] Failed to schedule %s: %r", schedule.name, e)
    finally:
        db_gen.close()
