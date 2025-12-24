"""Batch extraction service: create jobs, process documents, track progress."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List
from uuid import UUID

from amprenta_rag.database.models import ExtractedDocument, ExtractionJob
from amprenta_rag.extraction.parsers import get_parser
from amprenta_rag.extraction.structured_extractor import extract_from_parsed


class ExtractionBatchService:
    def __init__(self, db_session: Callable[[], Any]):
        # db_session is expected to be the context manager from amprenta_rag.database.session
        self._db_session = db_session

    def create_job(self, file_paths: List[str]) -> ExtractionJob:
        paths = [str(p) for p in file_paths if str(p).strip()]
        with self._db_session() as db:  # type: ignore[assignment]
            job = ExtractionJob(
                batch_id=None,
                file_count=len(paths),
                completed_count=0,
                status="pending",
                started_at=datetime.now(timezone.utc),
                completed_at=None,
            )
            db.add(job)
            db.flush()

            for p in paths:
                fp = Path(p)
                suffix = fp.suffix.lower().lstrip(".")
                doc = ExtractedDocument(
                    job_id=job.id,
                    file_path=str(fp),
                    original_filename=fp.name,
                    doc_type=suffix or "unknown",
                    extracted_entities={},  # required JSONB
                    extraction_config=None,
                    status="pending",
                    error_log=None,
                )
                db.add(doc)

            db.commit()
            db.refresh(job)
            return job

    def process_job(self, job_id: UUID) -> ExtractionJob:
        with self._db_session() as db:  # type: ignore[assignment]
            job = db.query(ExtractionJob).filter(ExtractionJob.id == job_id).first()
            if not job:
                raise ValueError("Job not found")

            job.status = "running"
            job.started_at = datetime.now(timezone.utc)
            job.completed_at = None
            job.completed_count = 0
            db.commit()
            db.refresh(job)

            docs = (
                db.query(ExtractedDocument)
                .filter(ExtractedDocument.job_id == job.id)
                .order_by(ExtractedDocument.created_at.asc())
                .all()
            )

            any_failed = False
            for doc in docs:
                try:
                    parser = get_parser(doc.file_path)
                    parsed = parser(doc.file_path)
                    res = extract_from_parsed(parsed, doc_type=doc.doc_type or "generic")

                    doc.extracted_entities = res.model_dump()
                    doc.status = "completed"
                    doc.error_log = None
                except Exception as e:  # noqa: BLE001
                    any_failed = True
                    doc.status = "failed"
                    doc.error_log = repr(e)
                    # keep extracted_entities best-effort
                    if not doc.extracted_entities:
                        doc.extracted_entities = {}

                job.completed_count = int(job.completed_count or 0) + 1
                db.commit()

            job.completed_at = datetime.now(timezone.utc)
            job.status = "failed" if any_failed else "completed"
            db.commit()
            db.refresh(job)
            return job

    def get_job_status(self, job_id: UUID) -> dict:
        with self._db_session() as db:  # type: ignore[assignment]
            job = db.query(ExtractionJob).filter(ExtractionJob.id == job_id).first()
            if not job:
                raise ValueError("Job not found")

            docs = (
                db.query(ExtractedDocument)
                .filter(ExtractedDocument.job_id == job.id)
                .order_by(ExtractedDocument.created_at.asc())
                .all()
            )

            file_count = int(job.file_count or 0) if job.file_count is not None else 0
            completed = int(job.completed_count or 0) if job.completed_count is not None else 0
            progress_pct = 0.0
            if file_count > 0:
                progress_pct = max(0.0, min(100.0, 100.0 * (completed / file_count)))

            def _job_dict(j: ExtractionJob) -> Dict[str, Any]:
                return {
                    "id": str(j.id),
                    "batch_id": str(j.batch_id) if j.batch_id else None,
                    "file_count": j.file_count,
                    "completed_count": j.completed_count,
                    "status": j.status,
                    "started_at": j.started_at.isoformat() if j.started_at else None,
                    "completed_at": j.completed_at.isoformat() if j.completed_at else None,
                    "created_at": j.created_at.isoformat() if j.created_at else None,
                    "updated_at": j.updated_at.isoformat() if j.updated_at else None,
                }

            def _doc_dict(d: ExtractedDocument) -> Dict[str, Any]:
                return {
                    "id": str(d.id),
                    "job_id": str(d.job_id),
                    "file_path": d.file_path,
                    "original_filename": d.original_filename,
                    "doc_type": d.doc_type,
                    "extracted_entities": d.extracted_entities or {},
                    "extraction_config": d.extraction_config,
                    "status": d.status,
                    "error_log": d.error_log,
                    "created_at": d.created_at.isoformat() if d.created_at else None,
                    "updated_at": d.updated_at.isoformat() if d.updated_at else None,
                }

            return {"job": _job_dict(job), "documents": [_doc_dict(d) for d in docs], "progress_pct": progress_pct}


