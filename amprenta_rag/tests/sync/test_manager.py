from __future__ import annotations

import asyncio
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Any, AsyncIterator, Callable, List, Optional
from uuid import UUID, uuid4

from amprenta_rag.database.models import SyncJob, SyncRecord
from amprenta_rag.sync.adapters.base import BaseSyncAdapter
from amprenta_rag.sync.manager import SyncManager


class _QueryMock:
    def __init__(self, session: "_SessionMock", model: Any):
        self._session = session
        self._model = model

    def filter(self, *args: Any, **kwargs: Any) -> "_QueryMock":  # noqa: ARG002
        return self

    def order_by(self, *args: Any, **kwargs: Any) -> "_QueryMock":  # noqa: ARG002
        return self

    def distinct(self, *args: Any, **kwargs: Any) -> "_QueryMock":  # noqa: ARG002
        return self

    def first(self) -> Any:
        if self._model is SyncJob:
            return self._session._syncjob_first()
        if self._model is SyncRecord:
            return self._session.sync_record
        return None

    def all(self) -> list[Any]:
        return []


class _SessionMock:
    def __init__(self, *, job: SyncJob, sync_record: Optional[SyncRecord]):
        self.job = job
        self.sync_record = sync_record
        self.added: List[Any] = []
        self.commit_calls = 0
        self.flush_calls = 0
        self.refresh_calls = 0
        self._syncjob_first_iter = iter([job, None])  # job lookup, then last_completed lookup

    def _syncjob_first(self) -> Any:
        try:
            return next(self._syncjob_first_iter)
        except StopIteration:
            return None

    def query(self, model: Any) -> _QueryMock:
        return _QueryMock(self, model)

    def add(self, obj: Any) -> None:
        self.added.append(obj)

    def flush(self) -> None:
        self.flush_calls += 1
        for obj in self.added:
            if hasattr(obj, "id") and getattr(obj, "id") is None:
                setattr(obj, "id", uuid4())

    def commit(self) -> None:
        self.commit_calls += 1

    def refresh(self, obj: Any) -> None:  # noqa: ARG002
        self.refresh_calls += 1


class _MockAdapter(BaseSyncAdapter):
    source = "chembl"

    def __init__(self, *, record: dict, checksum: str):
        self._record = record
        self._checksum = checksum
        self.map_calls = 0

    async def fetch_records(self, since: datetime | None) -> AsyncIterator[dict]:  # noqa: ARG002
        yield self._record

    def compute_checksum(self, record: dict) -> str:  # noqa: ARG002
        return self._checksum

    def map_to_entity(self, record: dict, db_session) -> tuple[str, UUID | None]:  # noqa: ARG002
        self.map_calls += 1
        return ("activity", None)


def _ctx_for(session: _SessionMock) -> Callable[[], Any]:
    @contextmanager
    def _cm():
        yield session

    return _cm


def test_create_sync_job() -> None:
    # Dummy placeholder session.job (create_sync_job does not query it, but our SessionMock requires one).
    job = SyncJob(
        id=uuid4(),
        source="chembl",
        sync_type="incremental",
        status="pending",
        records_synced=0,
        records_updated=0,
        records_new=0,
        conflicts_detected=0,
        started_at=None,
        completed_at=None,
        error_log=None,
        created_at=datetime.now(timezone.utc),
        updated_at=datetime.now(timezone.utc),
    )

    session = _SessionMock(job=job, sync_record=None)
    mgr = SyncManager(_ctx_for(session))
    created = mgr.create_sync_job("ChEMBL", sync_type="incremental")

    assert created.source == "chembl"
    assert created.status == "pending"
    assert session.commit_calls >= 1


def test_sync_record_upsert_updates_on_checksum_change() -> None:
    job = SyncJob(
        id=uuid4(),
        source="chembl",
        sync_type="incremental",
        status="pending",
        records_synced=0,
        records_updated=0,
        records_new=0,
        conflicts_detected=0,
        started_at=None,
        completed_at=None,
        error_log=None,
        created_at=datetime.now(timezone.utc),
        updated_at=datetime.now(timezone.utc),
    )

    existing = SyncRecord(
        id=uuid4(),
        job_id=None,
        source="chembl",
        external_id="CHEMBL_ACT_1",
        entity_type="activity",
        entity_id=None,
        checksum="old",
        synced_at=datetime.now(timezone.utc),
        metadata_={"external_id": "CHEMBL_ACT_1", "molecule_chembl_id": "CHEMBL1"},
    )

    session = _SessionMock(job=job, sync_record=existing)
    mgr = SyncManager(_ctx_for(session))
    adapter = _MockAdapter(record={"external_id": "CHEMBL_ACT_1", "molecule_chembl_id": "CHEMBL1"}, checksum="new")
    mgr.register_adapter(adapter)

    out = asyncio.run(mgr.run_sync_async(job.id))

    assert out.records_synced == 1
    assert out.records_updated == 1
    assert out.records_new == 0
    assert existing.checksum == "new"
    assert adapter.map_calls == 1


def test_sync_record_upsert_skips_when_checksum_unchanged() -> None:
    job = SyncJob(
        id=uuid4(),
        source="chembl",
        sync_type="incremental",
        status="pending",
        records_synced=0,
        records_updated=0,
        records_new=0,
        conflicts_detected=0,
        started_at=None,
        completed_at=None,
        error_log=None,
        created_at=datetime.now(timezone.utc),
        updated_at=datetime.now(timezone.utc),
    )

    existing = SyncRecord(
        id=uuid4(),
        job_id=None,
        source="chembl",
        external_id="CHEMBL_ACT_1",
        entity_type="activity",
        entity_id=None,
        checksum="same",
        synced_at=datetime.now(timezone.utc),
        metadata_={"external_id": "CHEMBL_ACT_1", "molecule_chembl_id": "CHEMBL1"},
    )

    session = _SessionMock(job=job, sync_record=existing)
    mgr = SyncManager(_ctx_for(session))
    adapter = _MockAdapter(record={"external_id": "CHEMBL_ACT_1", "molecule_chembl_id": "CHEMBL1"}, checksum="same")
    mgr.register_adapter(adapter)

    out = asyncio.run(mgr.run_sync_async(job.id))

    assert out.records_synced == 1
    assert out.records_updated == 0
    assert out.records_new == 0
    assert adapter.map_calls == 0


