"""Scheduled weekly executive digests (Papermill execution + HTML rendering).

This module mirrors the APScheduler pattern used by harvest/report scheduling.
"""

from __future__ import annotations

import logging
import os
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional
from uuid import UUID

try:
    from apscheduler.triggers.cron import CronTrigger

    APSCHEDULER_AVAILABLE = True
except ImportError:  # pragma: no cover
    CronTrigger = None  # type: ignore
    APSCHEDULER_AVAILABLE = False

from amprenta_rag.automation.notifications import send_digest_notification
from amprenta_rag.database.models import DigestSchedule, Program
from amprenta_rag.database.session import db_session
from amprenta_rag.ingestion.harvest_scheduler import get_scheduler, start_scheduler
from amprenta_rag.notebooks.registry import load_registry, resolve_notebook_path

logger = logging.getLogger(__name__)

try:
    import papermill as pm  # type: ignore

    PAPERMILL_AVAILABLE = True
except ImportError:
    pm = None  # type: ignore
    PAPERMILL_AVAILABLE = False

try:
    import nbformat  # type: ignore
    from nbconvert import HTMLExporter  # type: ignore

    NBCONVERT_AVAILABLE = True
except ImportError:
    nbformat = None  # type: ignore
    HTMLExporter = None  # type: ignore
    NBCONVERT_AVAILABLE = False


def _validate_notebook_path(nb_path: str) -> Path:
    p = (nb_path or "").strip()
    if not p:
        raise ValueError("notebook_path is required")
    for tpl in load_registry():
        if not isinstance(tpl, dict):
            continue
        if str(tpl.get("notebook_path") or "").strip() != p:
            continue
        nb = resolve_notebook_path(tpl)
        if not nb.exists():
            raise FileNotFoundError(str(nb))
        if nb.suffix != ".ipynb":
            raise ValueError("notebook must be .ipynb")
        return nb
    raise FileNotFoundError("notebook not found in registry")


def _cron_for_weekly(day_of_week: str, hour: int) -> str:
    # Cron: minute hour day-of-month month day-of-week
    # Use 0 minute for top of hour.
    dow = (day_of_week or "").strip().lower()
    if dow.isdigit():
        # 0-6 where 0=Sunday (crontab)
        return f"0 {int(hour)} * * {int(dow)}"
    # mon,tue,wed,thu,fri,sat,sun
    if dow not in {"mon", "tue", "wed", "thu", "fri", "sat", "sun"}:
        raise ValueError("day_of_week must be mon..sun or 0..6")
    return f"0 {int(hour)} * * {dow}"


@dataclass
class DigestRunResult:
    schedule_id: UUID
    output_ipynb: Path
    output_html: Optional[Path]
    status: str


class DigestScheduler:
    """Scheduler facade for weekly executive digests."""

    def schedule_weekly_digest(
        self,
        program_id: UUID,
        notebook_path: str,
        recipients: List[str],
        day_of_week: str,
        hour: int,
    ) -> DigestSchedule:
        """Create a DigestSchedule row and register it in APScheduler."""
        _validate_notebook_path(notebook_path)
        cron = _cron_for_weekly(day_of_week, hour)

        with db_session() as db:
            schedule_id = uuid.uuid4()
            sched = DigestSchedule(
                id=schedule_id,
                program_id=program_id,
                notebook_path=notebook_path.strip(),
                schedule_cron=cron,
                recipients=list(recipients or []),
                enabled=True,
            )
            db.add(sched)
            db.commit()
            db.refresh(sched)

        self.add_digest_schedule(schedule_id)
        return sched

    def add_digest_schedule(self, schedule_id: UUID) -> None:
        """Add a schedule to APScheduler (or replace existing)."""
        import os
        USE_CELERY = os.environ.get("USE_CELERY", "true").lower() == "true"
        
        if USE_CELERY:
            # Celery Beat handles digest schedules via periodic checks - just log
            logger.info("[DIGEST] Schedule %s managed by Celery Beat", schedule_id)
            return
        
        if not APSCHEDULER_AVAILABLE:
            raise ImportError("APScheduler is not installed")

        with db_session() as db:
            sched = db.query(DigestSchedule).filter(DigestSchedule.id == schedule_id).first()
            if not sched:
                logger.warning("[DIGEST] Schedule not found: %s", schedule_id)
                return
            if not sched.enabled:
                logger.info("[DIGEST] Schedule disabled, skipping add: %s", schedule_id)
                return
            cron_expr = str(sched.schedule_cron or "").strip()
            prog = db.query(Program).filter(Program.id == sched.program_id).first()
            prog_name = (prog.name if prog else "Program")

        scheduler = get_scheduler()
        start_scheduler()

        job_id = f"digest_{schedule_id}"
        if scheduler.get_job(job_id):
            scheduler.remove_job(job_id)

        trigger = CronTrigger.from_crontab(cron_expr)
        scheduler.add_job(
            self.run_digest,
            trigger=trigger,
            args=[schedule_id],
            id=job_id,
            name=f"Digest: {prog_name}",
            replace_existing=True,
        )
        logger.info("[DIGEST] Scheduled digest %s with cron '%s'", schedule_id, cron_expr)

    def remove_digest_schedule(self, schedule_id: UUID) -> None:
        if not APSCHEDULER_AVAILABLE:
            raise ImportError("APScheduler is not installed")
        scheduler = get_scheduler()
        job_id = f"digest_{schedule_id}"
        if scheduler.get_job(job_id):
            scheduler.remove_job(job_id)

    def run_digest(self, schedule_id: UUID) -> DigestRunResult:
        """Execute a digest immediately (scheduler entrypoint)."""
        with db_session() as db:
            sched = db.query(DigestSchedule).filter(DigestSchedule.id == schedule_id).first()
            if not sched:
                logger.warning("[DIGEST] Schedule not found: %s", schedule_id)
                return DigestRunResult(
                    schedule_id=schedule_id,
                    output_ipynb=Path(""),
                    output_html=None,
                    status="failed",
                )
            if not sched.enabled:
                logger.info("[DIGEST] Schedule disabled, skipping run: %s", schedule_id)
                return DigestRunResult(
                    schedule_id=schedule_id,
                    output_ipynb=Path(""),
                    output_html=None,
                    status="failed",
                )

            prog = db.query(Program).filter(Program.id == sched.program_id).first()
            program_name = (prog.name if prog else "Program")

            try:
                res = self._execute(schedule_id=schedule_id, program_id=sched.program_id, notebook_path=sched.notebook_path)
                sched.last_run_at = datetime.now(timezone.utc)
                sched.last_status = "success" if res.output_html else "failed"
                db.add(sched)
                db.commit()
            except Exception as e:  # noqa: BLE001
                logger.error("[DIGEST] Run failed for %s: %r", schedule_id, e)
                sched.last_run_at = datetime.now(timezone.utc)
                sched.last_status = "failed"
                db.add(sched)
                db.commit()
                return DigestRunResult(schedule_id=schedule_id, output_ipynb=Path(""), output_html=None, status="failed")

            # Notify on success
            if res.output_html:
                base_url = (  # full URL if configured, else relative
                    (os.environ.get("APP_BASE_URL") or os.environ.get("API_URL") or "").rstrip("/")
                )
                digest_url = f"{base_url}/api/digests/{schedule_id}/output" if base_url else f"/api/digests/{schedule_id}/output"
                try:
                    send_digest_notification(
                        recipients=list(sched.recipients or []),
                        digest_url=digest_url,
                        program_name=program_name,
                    )
                except Exception as e:  # noqa: BLE001
                    logger.warning("[DIGEST] Notification failed: %r", e)

            return DigestRunResult(
                schedule_id=schedule_id,
                output_ipynb=res.output_ipynb,
                output_html=res.output_html,
                status="success" if res.output_html else "failed",
            )

    def _execute(self, *, schedule_id: UUID, program_id: UUID, notebook_path: str) -> DigestRunResult:
        """Run papermill and render HTML. Writes latest.ipynb/html under data/digests/{schedule_id}/."""
        if not PAPERMILL_AVAILABLE:
            raise ImportError("papermill is not installed")
        if not NBCONVERT_AVAILABLE:
            raise ImportError("nbconvert is not installed")

        input_nb = _validate_notebook_path(notebook_path)
        root = Path(__file__).resolve().parents[2]  # .../RAG
        out_dir = root / "data" / "digests" / str(schedule_id)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_ipynb = out_dir / "latest.ipynb"
        out_html = out_dir / "latest.html"

        params = {"program_id": str(program_id), "executed_at": datetime.now(timezone.utc).isoformat()}
        pm.execute_notebook(str(input_nb), str(out_ipynb), parameters=params)  # type: ignore[union-attr]

        nb = nbformat.read(str(out_ipynb), as_version=4)  # type: ignore[union-attr]
        exporter = HTMLExporter()  # type: ignore[operator]
        body, _resources = exporter.from_notebook_node(nb)
        out_html.write_text(body, encoding="utf-8")

        return DigestRunResult(schedule_id=schedule_id, output_ipynb=out_ipynb, output_html=out_html, status="success")


__all__ = ["DigestScheduler", "DigestRunResult"]


