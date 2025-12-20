"""APScheduler integration for report generation schedules."""

from __future__ import annotations

import json
import logging
import hashlib
from datetime import datetime
from typing import Optional, cast
from uuid import UUID

from apscheduler.triggers.cron import CronTrigger

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import ReportSchedule
from amprenta_rag.ingestion.harvest_scheduler import get_scheduler, start_scheduler
from amprenta_rag.reports.generator import generate_report
from amprenta_rag.reports.artifact_registry import save_artifact

logger = logging.getLogger(__name__)


def _params_hash(entity_type: str, entity_id: Optional[UUID], fmt: str) -> str:
    payload = {"entity_type": entity_type, "entity_id": str(entity_id) if entity_id else None, "format": fmt}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def run_scheduled_report(schedule_id: UUID) -> None:
    """Execute a scheduled report and persist the artifact."""
    with db_session() as db:
        schedule = db.query(ReportSchedule).filter(ReportSchedule.id == schedule_id).first()
        if not schedule:
            logger.warning("[REPORT-SCHED] Schedule not found: %s", schedule_id)
            return
        if not schedule.enabled:
            logger.info("[REPORT-SCHED] Schedule disabled, skipping: %s", schedule.name)
            return

        entity_type = schedule.entity_type or ""
        entity_id = cast(Optional[UUID], schedule.entity_id)
        fmt = schedule.format or ""
        created_by = schedule.created_by_id
        schedule.last_run_at = datetime.utcnow()
        db.commit()

    params = {
        "entity_type": entity_type,
        "entity_id": str(entity_id) if entity_id else "",
        "format": fmt,
    }

    try:
        report_path = generate_report(
            template_name="narrative_report.ipynb",
            params=params,
            format=fmt,
        )
    except Exception as exc:
        logger.error("[REPORT-SCHED] Report generation failed for %s: %r", schedule_id, exc)
        return

    phash = _params_hash(entity_type, entity_id, fmt)
    if entity_id:
        save_artifact(
            entity_type=entity_type,
            entity_id=entity_id,
            format=fmt,
            file_path=report_path,
            params_hash=phash,
            user_id=cast(Optional[UUID], created_by) if created_by is not None else None,  # type: ignore[arg-type]
        )
    logger.info("[REPORT-SCHED] Generated report artifact for schedule %s", schedule_id)


def add_report_schedule(schedule_id: UUID, cron_expression: Optional[str] = None) -> None:
    """Add a single report schedule to the shared scheduler."""
    with db_session() as db:
        schedule = db.query(ReportSchedule).filter(ReportSchedule.id == schedule_id).first()
        if not schedule:
            logger.warning("[REPORT-SCHED] Schedule not found: %s", schedule_id)
            return
        if not schedule.enabled:
            logger.info("[REPORT-SCHED] Schedule disabled, skipping add: %s", schedule.name)
            return
        cron_expr = cron_expression or schedule.cron_expression
        name = schedule.name

    scheduler = get_scheduler()
    start_scheduler()

    job_id = f"report_{schedule_id}"
    if scheduler.get_job(job_id):
        scheduler.remove_job(job_id)

    trigger = CronTrigger.from_crontab(cron_expr)
    scheduler.add_job(
        run_scheduled_report,
        trigger=trigger,
        args=[schedule_id],
        id=job_id,
        name=f"Report: {name}",
        replace_existing=True,
    )
    logger.info("[REPORT-SCHED] Scheduled report %s with cron '%s'", name, cron_expr)


def remove_report_schedule(schedule_id: UUID) -> None:
    """Remove a scheduled report job."""
    scheduler = get_scheduler()
    job_id = f"report_{schedule_id}"
    if scheduler.get_job(job_id):
        scheduler.remove_job(job_id)
        logger.info("[REPORT-SCHED] Removed schedule job %s", job_id)


def load_report_schedules() -> None:
    """Load all enabled report schedules and add to scheduler."""
    with db_session() as db:
        schedules = db.query(ReportSchedule).filter(ReportSchedule.enabled == True).all()  # noqa: E712
        logger.info("[REPORT-SCHED] Loading %d enabled report schedules", len(schedules))
        for schedule in schedules:
            try:
                add_report_schedule(cast(UUID, schedule.id), schedule.cron_expression)  # type: ignore[arg-type]
            except Exception as exc:
                logger.error("[REPORT-SCHED] Failed to schedule %s: %r", schedule.name, exc)
