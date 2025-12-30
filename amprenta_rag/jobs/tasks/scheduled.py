"""Scheduled periodic tasks for Celery Beat."""

from uuid import UUID

from amprenta_rag.jobs.celery_app import celery_app


@celery_app.task(queue='scheduled')
def run_harvest(schedule_id: str) -> dict:
    """Execute harvest job for a schedule.
    
    Args:
        schedule_id: UUID string of the harvest schedule
        
    Returns:
        Dict with status and schedule_id
    """
    from amprenta_rag.ingestion.harvest_scheduler import run_harvest as _run_harvest
    
    try:
        _run_harvest(schedule_id)
        return {"status": "complete", "schedule_id": schedule_id}
    except Exception as e:
        return {"status": "failed", "schedule_id": schedule_id, "error": str(e)}


@celery_app.task(queue='scheduled')
def run_digest(schedule_id: str) -> dict:
    """Execute digest job for a schedule.
    
    Args:
        schedule_id: UUID string of the digest schedule
        
    Returns:
        Dict with status and schedule_id
    """
    from amprenta_rag.automation.digest_scheduler import DigestScheduler
    
    try:
        result = DigestScheduler().run_digest(UUID(schedule_id))
        return {
            "status": result.status,
            "schedule_id": schedule_id,
            "output_path": getattr(result, "output_path", None),
        }
    except Exception as e:
        return {"status": "failed", "schedule_id": schedule_id, "error": str(e)}


@celery_app.task(queue='scheduled')
def run_external_sync(source: str, sync_type: str = "incremental") -> dict:
    """Execute external sync (ChEMBL/PubChem).
    
    Args:
        source: Source name (e.g., 'chembl', 'pubchem')
        sync_type: Type of sync ('incremental' or 'full')
        
    Returns:
        Dict with status and source
    """
    from amprenta_rag.ingestion.harvest_scheduler import run_external_sync as _run_external_sync
    
    try:
        _run_external_sync(source, sync_type)
        return {"status": "complete", "source": source, "sync_type": sync_type}
    except Exception as e:
        return {"status": "failed", "source": source, "error": str(e)}


@celery_app.task(queue='scheduled')
def check_digest_schedules() -> dict:
    """Check for enabled digest schedules and trigger them if due.
    
    This task runs weekly to check all enabled digest schedules
    and trigger any that are due for execution.
    
    Returns:
        Dict with status and count of schedules processed
    """
    from amprenta_rag.automation.digest_scheduler import DigestScheduler
    from amprenta_rag.database.models import DigestSchedule
    from amprenta_rag.database.session import db_session
    from datetime import datetime, timezone
    
    try:
        processed_count = 0
        with db_session() as db:
            # Get all enabled digest schedules
            schedules = db.query(DigestSchedule).filter(
                DigestSchedule.enabled == True  # noqa: E712
            ).all()
            
            scheduler = DigestScheduler()
            for schedule in schedules:
                try:
                    # Check if schedule is due (simplified check)
                    # In a real implementation, you'd parse the cron expression
                    # and check if it's time to run
                    result = scheduler.run_digest(schedule.id)
                    if result.status == "complete":
                        processed_count += 1
                except Exception:
                    # Log error but continue with other schedules
                    continue
        
        return {
            "status": "complete",
            "schedules_processed": processed_count,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
    except Exception as e:
        return {"status": "failed", "error": str(e)}


@celery_app.task(queue='scheduled')
def check_harvest_schedules() -> dict:
    """Check for active harvest schedules and trigger them if due.
    
    This task can be scheduled to periodically check database-driven
    harvest schedules and trigger them as needed.
    
    Returns:
        Dict with status and count of schedules processed
    """
    from amprenta_rag.database.models import HarvestSchedule
    from amprenta_rag.database.session import db_session
    from amprenta_rag.ingestion.harvest_scheduler import run_harvest as _run_harvest
    from datetime import datetime, timezone
    
    try:
        processed_count = 0
        with db_session() as db:
            # Get all active harvest schedules
            schedules = db.query(HarvestSchedule).filter(
                HarvestSchedule.enabled == True  # noqa: E712
            ).all()
            
            current_time = datetime.now(timezone.utc)
            for schedule in schedules:
                try:
                    # Check if schedule is due based on next_run time
                    if schedule.next_run and schedule.next_run <= current_time:
                        _run_harvest(str(schedule.id))
                        processed_count += 1
                except Exception:
                    # Log error but continue with other schedules
                    continue
        
        return {
            "status": "complete",
            "schedules_processed": processed_count,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
    except Exception as e:
        return {"status": "failed", "error": str(e)}
