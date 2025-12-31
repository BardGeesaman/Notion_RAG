"""External source sync API endpoints."""

from __future__ import annotations

import asyncio
from datetime import datetime, timezone
from typing import Any, Dict, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field

from amprenta_rag.database.models import SyncConflict, SyncJob, SyncRecord
from amprenta_rag.database.session import db_session
from amprenta_rag.sync.manager import SyncManager
from amprenta_rag.sync.adapters.chembl import ChEMBLAdapter
from amprenta_rag.sync.adapters.pubchem import PubChemAdapter
from amprenta_rag.sync.adapters.geo import GEOSyncAdapter
from amprenta_rag.ingestion.repositories.geo import GEORepository


router = APIRouter(prefix="/sync", tags=["Sync"])


class RunSyncRequest(BaseModel):
    source: str = Field(..., min_length=1)
    sync_type: str = Field("incremental")


class ResolveConflictRequest(BaseModel):
    resolution: str = Field(..., pattern="^(auto_merged|manual_override|ignored)$")


def _make_manager() -> SyncManager:
    import os
    
    mgr = SyncManager(db_session)
    mgr.register_adapter(ChEMBLAdapter())
    mgr.register_adapter(PubChemAdapter(db_session))
    
    # GEO adapter with NCBI-compliant rate limiting
    geo_repo = GEORepository(
        api_key=os.getenv("GEO_API_KEY"),
        email=os.getenv("NCBI_EMAIL", "admin@example.com"),
    )
    mgr.register_adapter(GEOSyncAdapter(geo_repo))
    
    return mgr


async def _run_sync_in_background(job_id: UUID) -> None:
    mgr = _make_manager()
    await asyncio.to_thread(mgr.run_sync, job_id)


@router.post("/run")
async def run_sync(req: RunSyncRequest) -> Dict[str, Any]:
    mgr = _make_manager()
    try:
        job = mgr.create_sync_job(req.source, sync_type=req.sync_type)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    # Submit job via Celery or fallback to asyncio
    import os
    if os.environ.get("USE_CELERY", "true").lower() == "true":
        from amprenta_rag.jobs.tasks.sync import run_sync_job
        run_sync_job.delay(str(job.id))
    else:
        # Fallback to asyncio for gradual rollout
        asyncio.create_task(_run_sync_in_background(job.id))
    
    return {"job_id": str(job.id), "status": job.status}


@router.get("/jobs/{job_id}")
def get_job(job_id: UUID) -> Dict[str, Any]:
    mgr = _make_manager()
    try:
        return mgr.get_job_status(job_id)
    except ValueError:
        raise HTTPException(status_code=404, detail="Job not found")


@router.get("/jobs")
def list_jobs(
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=200),
    source: Optional[str] = Query(None),
) -> Dict[str, Any]:
    with db_session() as db:
        q = db.query(SyncJob).order_by(SyncJob.created_at.desc())
        if source:
            q = q.filter(SyncJob.source == source.strip().lower())
        total = q.count()
        jobs = q.offset(skip).limit(limit).all()

        def _job(j: SyncJob) -> Dict[str, Any]:
            return {
                "id": str(j.id),
                "source": j.source,
                "sync_type": j.sync_type,
                "status": j.status,
                "records_synced": int(j.records_synced or 0),
                "records_updated": int(j.records_updated or 0),
                "records_new": int(j.records_new or 0),
                "conflicts_detected": int(j.conflicts_detected or 0),
                "started_at": j.started_at.isoformat() if j.started_at else None,
                "completed_at": j.completed_at.isoformat() if j.completed_at else None,
                "created_at": j.created_at.isoformat() if j.created_at else None,
                "updated_at": j.updated_at.isoformat() if j.updated_at else None,
            }

        return {"total": total, "skip": skip, "limit": limit, "jobs": [_job(j) for j in jobs]}


@router.get("/conflicts")
def list_conflicts(
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=200),
    source: Optional[str] = Query(None),
    resolution_status: Optional[str] = Query("pending"),
) -> Dict[str, Any]:
    with db_session() as db:
        q = db.query(SyncConflict, SyncRecord).join(SyncRecord, SyncRecord.id == SyncConflict.record_id)
        if source:
            q = q.filter(SyncRecord.source == source.strip().lower())
        if resolution_status:
            q = q.filter(SyncConflict.resolution_status == resolution_status)

        total = q.count()
        rows = q.order_by(SyncConflict.id.desc()).offset(skip).limit(limit).all()

        def _row(c: SyncConflict, r: SyncRecord) -> Dict[str, Any]:
            return {
                "id": str(c.id),
                "record_id": str(c.record_id),
                "source": r.source,
                "external_id": r.external_id,
                "conflict_type": c.conflict_type,
                "resolution_status": c.resolution_status,
                "resolved_at": c.resolved_at.isoformat() if c.resolved_at else None,
                "local_value": c.local_value,
                "external_value": c.external_value,
            }

        return {"total": total, "skip": skip, "limit": limit, "conflicts": [_row(c, r) for c, r in rows]}


@router.post("/conflicts/{conflict_id}/resolve")
def resolve_conflict(conflict_id: UUID, req: ResolveConflictRequest) -> Dict[str, Any]:
    with db_session() as db:
        conflict = db.query(SyncConflict).filter(SyncConflict.id == conflict_id).first()
        if not conflict:
            raise HTTPException(status_code=404, detail="Conflict not found")

        conflict.resolution_status = req.resolution
        conflict.resolved_at = datetime.now(timezone.utc)
        db.commit()

        return {
            "id": str(conflict.id),
            "resolution_status": conflict.resolution_status,
            "resolved_at": conflict.resolved_at.isoformat() if conflict.resolved_at else None,
        }


