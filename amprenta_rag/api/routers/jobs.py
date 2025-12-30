"""Job queue management API endpoints."""

from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.database.models import (
    PipelineJob,
    DockingRun,
    ExtractionJob,
    SyncJob,
    SingleCellDataset,
)
from amprenta_rag.database.session import db_session

router = APIRouter(prefix="/jobs", tags=["Jobs"])

# Map job types to model classes and task names
JOB_TYPE_MAP = {
    "genomics": {
        "model": PipelineJob,
        "task": "amprenta_rag.jobs.tasks.genomics.run_genomics_pipeline",
        "id_field": "id",
        "status_field": "status",
    },
    "docking": {
        "model": DockingRun,
        "task": "amprenta_rag.jobs.tasks.docking.run_docking",
        "id_field": "id",
        "status_field": "status",
    },
    "extraction": {
        "model": ExtractionJob,
        "task": "amprenta_rag.jobs.tasks.extraction.process_extraction_job",
        "id_field": "id",
        "status_field": "status",
    },
    "sync": {
        "model": SyncJob,
        "task": "amprenta_rag.jobs.tasks.sync.run_sync_job",
        "id_field": "id",
        "status_field": "status",
    },
    "single_cell": {
        "model": SingleCellDataset,
        "task": "amprenta_rag.jobs.tasks.single_cell.process_single_cell",
        "id_field": "id",
        "status_field": "processing_status",
    },
}


def _get_job_config(job_type: str) -> Dict[str, Any]:
    """Get job configuration or raise 400 error."""
    if job_type not in JOB_TYPE_MAP:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown job type: {job_type}. Valid types: {list(JOB_TYPE_MAP.keys())}"
        )
    return JOB_TYPE_MAP[job_type]


def _serialize_job(job: Any, job_type: str) -> Dict[str, Any]:
    """Serialize job model to dict."""
    config = JOB_TYPE_MAP[job_type]
    status_field = config["status_field"]
    
    base = {
        "id": str(job.id),
        "type": job_type,
        "status": getattr(job, status_field, None),
        "created_at": getattr(job, "created_at", None),
        "started_at": getattr(job, "started_at", None),
        "completed_at": getattr(job, "completed_at", None),
    }
    
    # Add type-specific fields
    if job_type == "genomics":
        base.update({
            "tool": getattr(job, "tool", None),
            "progress_percent": getattr(job, "progress_percent", None),
            "error_message": getattr(job, "error_message", None),
        })
    elif job_type == "docking":
        base.update({
            "total_compounds": getattr(job, "total_compounds", None),
            "completed_compounds": getattr(job, "completed_compounds", None),
            "error_log": getattr(job, "error_log", None),
        })
    elif job_type == "extraction":
        base.update({
            "file_count": getattr(job, "file_count", None),
            "completed_count": getattr(job, "completed_count", None),
        })
    elif job_type == "sync":
        base.update({
            "source": getattr(job, "source", None),
            "sync_type": getattr(job, "sync_type", None),
            "records_synced": getattr(job, "records_synced", None),
            "records_new": getattr(job, "records_new", None),
            "records_updated": getattr(job, "records_updated", None),
        })
    elif job_type == "single_cell":
        base.update({
            "n_cells": getattr(job, "n_cells", None),
            "n_genes": getattr(job, "n_genes", None),
            "processing_log": getattr(job, "processing_log", None),
        })
    
    return base


@router.get("")
def list_jobs(
    job_type: Optional[str] = Query(None, description="Filter by job type"),
    status: Optional[str] = Query(None, description="Filter by status"),
    skip: int = Query(0, ge=0, description="Number of jobs to skip"),
    limit: int = Query(50, ge=1, le=200, description="Number of jobs to return"),
) -> Dict[str, Any]:
    """List jobs with optional filters."""
    jobs: List[Dict[str, Any]] = []
    
    with db_session() as db:
        # If job_type specified, query only that type
        if job_type:
            config = _get_job_config(job_type)
            model = config["model"]
            status_field = config["status_field"]
            
            query = db.query(model)
            if status:
                query = query.filter(getattr(model, status_field) == status)
            
            results = query.order_by(getattr(model, "created_at", model.id).desc()).offset(skip).limit(limit).all()
            jobs.extend([_serialize_job(job, job_type) for job in results])
        
        else:
            # Query all job types
            for jt, config in JOB_TYPE_MAP.items():
                model = config["model"]
                status_field = config["status_field"]
                
                query = db.query(model)
                if status:
                    query = query.filter(getattr(model, status_field) == status)
                
                # Get a portion of each type (simple approach)
                type_limit = max(1, limit // len(JOB_TYPE_MAP))
                results = query.order_by(getattr(model, "created_at", model.id).desc()).offset(skip).limit(type_limit).all()
                jobs.extend([_serialize_job(job, jt) for job in results])
    
    # Sort by created_at if available
    jobs.sort(key=lambda x: x.get("created_at") or "", reverse=True)
    
    return {
        "jobs": jobs[:limit],
        "total": len(jobs),
        "skip": skip,
        "limit": limit,
    }


@router.get("/{job_type}/{job_id}")
def get_job(job_type: str, job_id: UUID) -> Dict[str, Any]:
    """Get job status by type and ID."""
    config = _get_job_config(job_type)
    model = config["model"]
    
    with db_session() as db:
        job = db.query(model).filter(getattr(model, "id") == job_id).first()
        if not job:
            raise HTTPException(
                status_code=404,
                detail=f"{job_type.title()} job {job_id} not found"
            )
        
        return _serialize_job(job, job_type)


@router.post("/{job_type}/{job_id}/cancel")
def cancel_job(job_type: str, job_id: UUID) -> Dict[str, Any]:
    """Cancel pending/running job."""
    config = _get_job_config(job_type)
    model = config["model"]
    status_field = config["status_field"]
    
    with db_session() as db:
        job = db.query(model).filter(getattr(model, "id") == job_id).first()
        if not job:
            raise HTTPException(
                status_code=404,
                detail=f"{job_type.title()} job {job_id} not found"
            )
        
        current_status = getattr(job, status_field)
        if current_status not in ["pending", "running"]:
            raise HTTPException(
                status_code=400,
                detail=f"Cannot cancel job with status: {current_status}"
            )
        
        # Update status to cancelled
        setattr(job, status_field, "cancelled")
        db.add(job)
        db.commit()
        
        # Try to revoke Celery task (best effort)
        try:
            from amprenta_rag.jobs.celery_app import celery_app
            # Note: This is simplified - in practice you'd need to track Celery task IDs
            # celery_app.control.revoke(task_id, terminate=True)
        except Exception:
            pass  # Celery revocation is best effort
        
        return {"message": f"Job {job_id} cancelled", "status": "cancelled"}


@router.post("/{job_type}/{job_id}/retry")
def retry_job(job_type: str, job_id: UUID) -> Dict[str, Any]:
    """Retry failed job."""
    config = _get_job_config(job_type)
    model = config["model"]
    status_field = config["status_field"]
    task_name = config["task"]
    
    with db_session() as db:
        job = db.query(model).filter(getattr(model, "id") == job_id).first()
        if not job:
            raise HTTPException(
                status_code=404,
                detail=f"{job_type.title()} job {job_id} not found"
            )
        
        current_status = getattr(job, status_field)
        if current_status != "failed":
            raise HTTPException(
                status_code=400,
                detail=f"Cannot retry job with status: {current_status}. Only failed jobs can be retried."
            )
        
        # Reset status to pending
        setattr(job, status_field, "pending")
        if hasattr(job, "error_message"):
            job.error_message = None
        if hasattr(job, "error_log"):
            job.error_log = None
        if hasattr(job, "processing_log"):
            job.processing_log = None
        db.add(job)
        db.commit()
        
        # Resubmit to Celery
        try:
            from amprenta_rag.jobs.celery_app import celery_app
            # Get the task function and submit
            task_parts = task_name.split(".")
            module_path = ".".join(task_parts[:-1])
            task_func_name = task_parts[-1]
            
            module = __import__(module_path, fromlist=[task_func_name])
            task_func = getattr(module, task_func_name)
            task_func.delay(str(job_id))
            
        except Exception as e:
            # Rollback status if task submission failed
            setattr(job, status_field, "failed")
            if hasattr(job, "error_message"):
                job.error_message = f"Failed to resubmit task: {e}"
            db.add(job)
            db.commit()
            raise HTTPException(
                status_code=500,
                detail=f"Failed to resubmit job: {e}"
            )
        
        return {"message": f"Job {job_id} resubmitted", "status": "pending"}
