"""Sync orchestration for external sources (MVP)."""

from __future__ import annotations

import asyncio
from datetime import datetime, timezone
from typing import Any, Callable, Dict, Optional
from uuid import UUID

from amprenta_rag.database.models import SyncConflict, SyncJob, SyncRecord
from amprenta_rag.sync.adapters.base import BaseSyncAdapter


class SyncManager:
    def __init__(self, db_session: Callable[[], Any]):
        # db_session is expected to be the context manager from amprenta_rag.database.session
        self._db_session = db_session
        self.adapters: Dict[str, BaseSyncAdapter] = {}

    def register_adapter(self, adapter: BaseSyncAdapter) -> None:
        src = (adapter.source or "").strip().lower()
        if not src:
            raise ValueError("Adapter must define a non-empty `source` string")
        self.adapters[src] = adapter

    def create_sync_job(self, source: str, sync_type: str = "incremental") -> SyncJob:
        src = (source or "").strip().lower()
        if not src:
            raise ValueError("source is required")

        with self._db_session() as db:  # type: ignore[assignment]
            job = SyncJob(
                source=src,
                sync_type=sync_type,
                status="pending",
                records_synced=0,
                records_updated=0,
                records_new=0,
                conflicts_detected=0,
                started_at=None,
                completed_at=None,
                error_log=None,
            )
            db.add(job)
            db.commit()
            db.refresh(job)
            return job

    def run_sync(self, job_id: UUID) -> SyncJob:
        """Run the sync in a synchronous context.

        NOTE: This wrapper uses `asyncio.run()` and will raise if called from a running event loop.
        In async contexts, call `await run_sync_async(...)`.
        """

        try:
            asyncio.get_running_loop()
        except RuntimeError:
            return asyncio.run(self.run_sync_async(job_id))
        raise RuntimeError("SyncManager.run_sync() cannot be called from an event loop; use await run_sync_async().")

    async def run_sync_async(self, job_id: UUID) -> SyncJob:
        with self._db_session() as db:  # type: ignore[assignment]
            job: Optional[SyncJob] = db.query(SyncJob).filter(SyncJob.id == job_id).first()
            if not job:
                raise ValueError("Job not found")

            adapter = self.adapters.get((job.source or "").lower())
            if adapter is None:
                raise ValueError(f"No adapter registered for source={job.source!r}")

            # determine "since" from most recent completed job
            last_completed: Optional[SyncJob] = (
                db.query(SyncJob)
                .filter(SyncJob.source == job.source, SyncJob.status == "completed", SyncJob.completed_at.isnot(None))
                .order_by(SyncJob.completed_at.desc())
                .first()
            )
            since = last_completed.completed_at if last_completed else None

            job.status = "running"
            job.started_at = datetime.now(timezone.utc)
            job.completed_at = None
            job.error_log = None
            job.records_synced = 0
            job.records_updated = 0
            job.records_new = 0
            job.conflicts_detected = 0
            db.commit()
            db.refresh(job)

            any_failed = False
            try:
                async for record in adapter.fetch_records(since=since):
                    job.records_synced = int(job.records_synced or 0) + 1

                    external_id = (record.get("external_id") or "").strip()
                    if not external_id:
                        any_failed = True
                        job.conflicts_detected = int(job.conflicts_detected or 0) + 1
                        # create a lightweight conflict without a SyncRecord
                        job.error_log = (job.error_log or "") + "\nMissing external_id in record"
                        db.commit()
                        continue

                    checksum = adapter.compute_checksum(record)

                    sync_rec: Optional[SyncRecord] = (
                        db.query(SyncRecord)
                        .filter(SyncRecord.source == job.source, SyncRecord.external_id == external_id)
                        .first()
                    )

                    if sync_rec is None:
                        sync_rec = SyncRecord(
                            job_id=job.id,
                            source=job.source,
                            external_id=external_id,
                            entity_type="unknown",
                            entity_id=None,
                            checksum=checksum,
                            synced_at=datetime.now(timezone.utc),
                            metadata_=record,
                        )
                        db.add(sync_rec)
                        job.records_new = int(job.records_new or 0) + 1
                        db.flush()
                    else:
                        # existing record
                        sync_rec.job_id = job.id
                        if (sync_rec.checksum or "") != checksum:
                            sync_rec.checksum = checksum
                            sync_rec.synced_at = datetime.now(timezone.utc)
                            sync_rec.metadata_ = record
                            job.records_updated = int(job.records_updated or 0) + 1
                        else:
                            # unchanged -> skip mapping work
                            db.commit()
                            continue

                    # map to internal entity (best-effort)
                    try:
                        entity_type, entity_id = adapter.map_to_entity(record, db)
                        sync_rec.entity_type = entity_type
                        sync_rec.entity_id = entity_id
                    except Exception as e:  # noqa: BLE001
                        any_failed = True
                        job.conflicts_detected = int(job.conflicts_detected or 0) + 1
                        db.add(
                            SyncConflict(
                                record_id=sync_rec.id,
                                conflict_type="schema_mismatch",
                                local_value={"error": repr(e)},
                                external_value=record,
                                resolution_status="pending",
                                resolved_at=None,
                            )
                        )

                    db.commit()

                job.status = "failed" if any_failed else "completed"
                job.completed_at = datetime.now(timezone.utc)
                db.commit()
                db.refresh(job)
                return job
            except Exception as e:  # noqa: BLE001
                job.status = "failed"
                job.completed_at = datetime.now(timezone.utc)
                job.error_log = repr(e)
                db.commit()
                db.refresh(job)
                return job

    def get_job_status(self, job_id: UUID) -> dict:
        with self._db_session() as db:  # type: ignore[assignment]
            job = db.query(SyncJob).filter(SyncJob.id == job_id).first()
            if not job:
                raise ValueError("Job not found")

            def _job_dict(j: SyncJob) -> Dict[str, Any]:
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
                    "error_log": j.error_log,
                    "created_at": j.created_at.isoformat() if j.created_at else None,
                    "updated_at": j.updated_at.isoformat() if j.updated_at else None,
                }

            return {"job": _job_dict(job)}


