"""Celery tasks for external source synchronization."""

from uuid import UUID

from amprenta_rag.jobs.celery_app import celery_app


@celery_app.task(bind=True, max_retries=3, default_retry_delay=300, queue='low')
def run_sync_job(self, job_id: str) -> dict:
    """Execute external source synchronization (ChEMBL/PubChem)."""
    from amprenta_rag.sync.manager import SyncManager
    from amprenta_rag.database.session import db_session
    from amprenta_rag.sync.adapters.chembl import ChEMBLAdapter
    from amprenta_rag.sync.adapters.pubchem import PubChemAdapter
    
    job_uuid = UUID(job_id)
    
    try:
        # Create manager and register adapters
        mgr = SyncManager(db_session)
        mgr.register_adapter(ChEMBLAdapter())
        mgr.register_adapter(PubChemAdapter(db_session))
        
        # Run sync job
        job = mgr.run_sync(job_uuid)
        
        return {
            "status": job.status,
            "job_id": job_id,
            "records_synced": job.records_synced,
            "records_new": job.records_new,
            "records_updated": job.records_updated,
            "conflicts_detected": job.conflicts_detected
        }
    
    except Exception as exc:
        # Update job status to failed if possible
        try:
            from amprenta_rag.database.models import SyncJob
            from datetime import datetime, timezone
            
            with db_session() as db:
                job = db.query(SyncJob).filter(SyncJob.id == job_uuid).first()
                if job is not None:
                    job.status = "failed"
                    job.completed_at = datetime.now(timezone.utc)
                    job.error_log = str(exc)
                    db.add(job)
                    db.commit()
        except Exception:
            pass  # Don't let DB errors prevent retry logic
        
        # Retry with exponential backoff
        if self.request.retries >= self.max_retries:
            return {"status": "failed", "error": str(exc), "job_id": job_id}
        
        self.retry(exc=exc, countdown=300 * (2 ** self.request.retries))
