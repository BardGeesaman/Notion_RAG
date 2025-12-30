"""Celery tasks for document extraction processing."""

from uuid import UUID

from amprenta_rag.jobs.celery_app import celery_app


@celery_app.task(bind=True, max_retries=3, default_retry_delay=60, queue='high')
def process_extraction_job(self, job_id: str) -> dict:
    """Process batch document extraction."""
    from amprenta_rag.extraction.batch_service import ExtractionBatchService
    from amprenta_rag.database.session import db_session
    
    job_uuid = UUID(job_id)
    
    try:
        # Create service and process job
        svc = ExtractionBatchService(db_session)
        job = svc.process_job(job_uuid)
        
        return {
            "status": job.status,
            "job_id": job_id,
            "completed_count": job.completed_count,
            "file_count": job.file_count
        }
    
    except Exception as exc:
        # Update job status to failed if possible
        try:
            from amprenta_rag.database.models import ExtractionJob
            from datetime import datetime, timezone
            
            with db_session() as db:
                job = db.query(ExtractionJob).filter(ExtractionJob.id == job_uuid).first()
                if job is not None:
                    job.status = "failed"
                    job.completed_at = datetime.now(timezone.utc)
                    db.add(job)
                    db.commit()
        except Exception:
            pass  # Don't let DB errors prevent retry logic
        
        # Retry with exponential backoff
        if self.request.retries >= self.max_retries:
            return {"status": "failed", "error": str(exc), "job_id": job_id}
        
        self.retry(exc=exc, countdown=60 * (2 ** self.request.retries))
