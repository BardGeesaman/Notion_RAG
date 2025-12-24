"""Batch extraction API endpoints."""

from __future__ import annotations

import asyncio
import tempfile
from pathlib import Path
from typing import Any, Dict, List
from uuid import UUID

from fastapi import APIRouter, File, HTTPException, Query, UploadFile

from amprenta_rag.database.models import ExtractedDocument, ExtractionJob
from amprenta_rag.database.session import db_session
from amprenta_rag.extraction.batch_service import ExtractionBatchService


router = APIRouter(prefix="/extraction", tags=["Extraction"])


async def _run_job_in_background(job_id: UUID) -> None:
    svc = ExtractionBatchService(db_session)
    await asyncio.to_thread(svc.process_job, job_id)


@router.post("/upload-batch")
async def upload_batch(files: List[UploadFile] = File(...)) -> Dict[str, Any]:
    if not files:
        raise HTTPException(status_code=400, detail="files[] is required")

    tmpdir = Path(tempfile.mkdtemp(prefix="amprenta_extract_"))
    saved: List[str] = []
    for f in files:
        try:
            content = await f.read()
        except Exception as e:  # noqa: BLE001
            raise HTTPException(status_code=400, detail=f"Failed to read upload: {e}")

        name = (f.filename or "upload").replace("/", "_")
        path = tmpdir / name
        try:
            path.write_bytes(content)
        except Exception as e:  # noqa: BLE001
            raise HTTPException(status_code=500, detail=f"Failed to save upload: {e}")
        saved.append(str(path))

    svc = ExtractionBatchService(db_session)
    job = svc.create_job(saved)

    # Kick off background job processing
    asyncio.create_task(_run_job_in_background(job.id))

    return {"job_id": str(job.id), "status": job.status, "file_count": job.file_count}


@router.get("/jobs/{job_id}")
def get_job(job_id: UUID) -> Dict[str, Any]:
    svc = ExtractionBatchService(db_session)
    try:
        return svc.get_job_status(job_id)
    except ValueError:
        raise HTTPException(status_code=404, detail="Job not found")


@router.get("/jobs")
def list_jobs(
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=200),
) -> Dict[str, Any]:
    with db_session() as db:
        q = db.query(ExtractionJob).order_by(ExtractionJob.created_at.desc())
        total = q.count()
        jobs = q.offset(skip).limit(limit).all()

        ids = [j.id for j in jobs]
        docs_by_job: Dict[str, int] = {}
        if ids:
            # count documents per job for quick summary
            rows = (
                db.query(ExtractedDocument.job_id, db.func.count(ExtractedDocument.id))  # type: ignore[attr-defined]
                .filter(ExtractedDocument.job_id.in_(ids))
                .group_by(ExtractedDocument.job_id)
                .all()
            )
            docs_by_job = {str(jid): int(cnt) for jid, cnt in rows}

        def _job(j: ExtractionJob) -> Dict[str, Any]:
            file_count = int(j.file_count or 0) if j.file_count is not None else 0
            completed = int(j.completed_count or 0) if j.completed_count is not None else 0
            pct = 0.0
            if file_count > 0:
                pct = max(0.0, min(100.0, 100.0 * (completed / file_count)))
            return {
                "id": str(j.id),
                "batch_id": str(j.batch_id) if j.batch_id else None,
                "file_count": j.file_count,
                "completed_count": j.completed_count,
                "status": j.status,
                "progress_pct": pct,
                "document_count": docs_by_job.get(str(j.id), 0),
                "created_at": j.created_at.isoformat() if j.created_at else None,
                "updated_at": j.updated_at.isoformat() if j.updated_at else None,
            }

        return {"total": total, "skip": skip, "limit": limit, "jobs": [_job(j) for j in jobs]}


