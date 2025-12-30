"""Celery tasks for genomics pipeline processing."""

from datetime import datetime
from pathlib import Path
from uuid import UUID

from amprenta_rag.jobs.celery_app import celery_app


@celery_app.task(bind=True, max_retries=3, default_retry_delay=60, queue='default')
def run_genomics_pipeline(self, job_id: str) -> dict:
    """Execute Salmon/Kallisto quantification pipeline."""
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import PipelineJob, GenomicsIndex
    from amprenta_rag.ingestion.genomics import pipeline
    from amprenta_rag.config import get_config
    
    job_uuid = UUID(job_id)
    JOB_ROOT = Path(get_config().genomics.job_root)
    
    try:
        # Mark running
        with db_session() as db:
            job = db.query(PipelineJob).filter(PipelineJob.id == job_uuid).first()
            if job is None:
                return {"status": "failed", "error": "Job not found", "job_id": job_id}
            job.status = "running"
            job.started_at = datetime.utcnow()
            job.progress_percent = 5
            db.add(job)
            db.commit()

        # Get job details and run pipeline
        with db_session() as db:
            job = db.query(PipelineJob).filter(PipelineJob.id == job_uuid).first()
            if job is None:
                return {"status": "failed", "error": "Job not found", "job_id": job_id}
            idx = db.query(GenomicsIndex).filter(GenomicsIndex.id == job.index_id).first()
            if idx is None:
                raise ValueError(f"Index not found for job: {job.index_id}")

            tool = (job.tool or "").lower()
            fastq = Path(job.input_fastq_path or "")
            index_path = Path(idx.file_path)
            out_dir = Path(job.output_dir or (JOB_ROOT / str(job_uuid)))
            out_dir.mkdir(parents=True, exist_ok=True)

            # Progress: about to start quantification
            job.progress_percent = 50
            db.add(job)
            db.commit()

        # Run the appropriate pipeline
        if tool == "salmon":
            result = pipeline.run_salmon_quant(fastq_path=fastq, index_path=index_path, output_dir=out_dir)
        elif tool == "kallisto":
            result = pipeline.run_kallisto_quant(fastq_path=fastq, index_path=index_path, output_dir=out_dir)
        else:
            raise ValueError(f"Unknown tool: {tool}")

        if result is None:
            raise RuntimeError("Pipeline returned no result file.")

        # Mark complete
        with db_session() as db:
            job = db.query(PipelineJob).filter(PipelineJob.id == job_uuid).first()
            if job is None:
                return {"status": "failed", "error": "Job not found", "job_id": job_id}
            job.status = "complete"
            job.progress_percent = 100
            job.completed_at = datetime.utcnow()
            job.result_file = str(result)
            job.error_message = None
            db.add(job)
            db.commit()

        return {"status": "complete", "job_id": job_id, "result_file": str(result)}

    except Exception as exc:
        # Update job status to failed
        try:
            with db_session() as db:
                job = db.query(PipelineJob).filter(PipelineJob.id == job_uuid).first()
                if job is not None:
                    job.status = "failed"
                    job.completed_at = datetime.utcnow()
                    job.error_message = str(exc)
                    job.progress_percent = min(int(job.progress_percent or 0), 99)
                    db.add(job)
                    db.commit()
        except Exception:
            pass  # Don't let DB errors prevent retry logic
        
        # Retry with exponential backoff
        if self.request.retries >= self.max_retries:
            return {"status": "failed", "error": str(exc), "job_id": job_id}
        
        self.retry(exc=exc, countdown=60 * (2 ** self.request.retries))
