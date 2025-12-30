"""Async job runner for genomics pipeline executions.

Persists job state in the database and runs Salmon/Kallisto quantification in a
background daemon thread.
"""

from __future__ import annotations

import threading
import uuid
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from amprenta_rag.database.models import GenomicsIndex, PipelineJob
from amprenta_rag.database.session import db_session
from amprenta_rag.ingestion.genomics import pipeline


JOB_ROOT = Path("/data/genomics/jobs")


def submit_job(tool: str, fastq_path: str, index_id: str, user_id: str) -> uuid.UUID:
    """Create a new job and kick off async execution. Returns job_id."""
    job_id = uuid.uuid4()
    out_dir = JOB_ROOT / str(job_id)
    out_dir_str = str(out_dir)

    with db_session() as db:
        job = PipelineJob(
            id=job_id,
            status="pending",
            tool=tool,
            input_fastq_path=fastq_path,
            index_id=index_id,
            output_dir=out_dir_str,
            progress_percent=0,
            created_by=user_id,
            created_at=datetime.utcnow(),
        )
        db.add(job)

    # Submit job via Celery or fallback to threading
    import os
    if os.environ.get("USE_CELERY", "true").lower() == "true":
        from amprenta_rag.jobs.tasks.genomics import run_genomics_pipeline
        run_genomics_pipeline.delay(str(job_id))
    else:
        # Fallback to threading for gradual rollout
        threading.Thread(target=_run_job_async, args=(job_id,), daemon=True).start()
    
    return job_id


def get_job_status(job_id: str) -> PipelineJob:
    """Fetch a job by id; raises if not found."""
    with db_session() as db:
        job = db.query(PipelineJob).filter(PipelineJob.id == job_id).first()
        if job is None:
            raise ValueError(f"Job not found: {job_id}")
        return job


def list_jobs(user_id: Optional[str] = None, limit: int = 20) -> List[PipelineJob]:
    """List recent jobs (optionally filtered by creator)."""
    with db_session() as db:
        q = db.query(PipelineJob)
        if user_id:
            q = q.filter(PipelineJob.created_by == user_id)
        return q.order_by(PipelineJob.created_at.desc()).limit(limit).all()


def _run_job_async(job_id: uuid.UUID) -> None:
    """Background runner: updates job state and executes the pipeline."""
    # Mark running
    with db_session() as db:
        job = db.query(PipelineJob).filter(PipelineJob.id == job_id).first()
        if job is None:
            return
        job.status = "running"
        job.started_at = datetime.utcnow()
        job.progress_percent = 5
        db.add(job)

    try:
        with db_session() as db:
            job = db.query(PipelineJob).filter(PipelineJob.id == job_id).first()
            if job is None:
                return
            idx = db.query(GenomicsIndex).filter(GenomicsIndex.id == job.index_id).first()
            if idx is None:
                raise ValueError(f"Index not found for job: {job.index_id}")

            tool = (job.tool or "").lower()
            fastq = Path(job.input_fastq_path or "")
            index_path = Path(idx.file_path)
            out_dir = Path(job.output_dir or (JOB_ROOT / str(job_id)))
            out_dir.mkdir(parents=True, exist_ok=True)

            # Progress: about to start quantification.
            job.progress_percent = 50
            db.add(job)

        if tool == "salmon":
            result = pipeline.run_salmon_quant(fastq_path=fastq, index_path=index_path, output_dir=out_dir)
        elif tool == "kallisto":
            result = pipeline.run_kallisto_quant(fastq_path=fastq, index_path=index_path, output_dir=out_dir)
        else:
            raise ValueError(f"Unknown tool: {tool}")

        if result is None:
            raise RuntimeError("Pipeline returned no result file.")

        with db_session() as db:
            job = db.query(PipelineJob).filter(PipelineJob.id == job_id).first()
            if job is None:
                return
            job.status = "complete"
            job.progress_percent = 100
            job.completed_at = datetime.utcnow()
            job.result_file = str(result)
            job.error_message = None
            db.add(job)

    except Exception as e:  # noqa: BLE001
        with db_session() as db:
            job = db.query(PipelineJob).filter(PipelineJob.id == job_id).first()
            if job is None:
                return
            job.status = "failed"
            job.completed_at = datetime.utcnow()
            job.error_message = str(e)
            job.progress_percent = min(int(job.progress_percent or 0), 99)
            db.add(job)


__all__ = [
    "submit_job",
    "get_job_status",
    "list_jobs",
]


