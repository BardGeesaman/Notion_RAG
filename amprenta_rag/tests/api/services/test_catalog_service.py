from __future__ import annotations

from datetime import datetime, timedelta, timezone
from typing import List
from uuid import uuid4

import pytest

from amprenta_rag.api.services import catalog_service


class AttrField:
    def __init__(self, name: str):
        self.name = name

    def ilike(self, pattern: str):
        return ("ilike", self.name, pattern)

    def desc(self):
        return self


class FakeDataset:
    id = "id"
    data_origin = AttrField("data_origin")
    features = AttrField("features")
    name = AttrField("name")
    description = AttrField("description")
    updated_at = AttrField("updated_at")
    created_at = AttrField("created_at")

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.data_origin = kwargs.get("data_origin", "")
        self.updated_at = kwargs.get("updated_at")
        self.created_at = kwargs.get("created_at")
        self.name = kwargs.get("name")
        self.description = kwargs.get("description")
        self.ingestion_status = kwargs.get("ingestion_status")
        self.external_ids = kwargs.get("external_ids", {})
        self.features = kwargs.get("features", [])


class FakeRow:
    def __init__(self, dataset, total):
        self._dataset = dataset
        self.total_count = total

    def __getitem__(self, idx):
        if idx == 0:
            return self._dataset
        raise IndexError


class FakeQuery:
    def __init__(self, data: List[FakeDataset]):
        self._data = list(data)

    def options(self, *_):
        return self

    def filter(self, *predicates):
        filtered = []
        for d in self._data:
            ok = True
            for pred in predicates:
                if isinstance(pred, tuple):
                    if len(pred) == 3 and pred[0] == "ilike":
                        field, pattern = pred[1], pred[2]
                        term = pattern.strip("%").lower()
                        if term not in (getattr(d, field, "") or "").lower():
                            ok = False
                            break
                    if len(pred) == 3 and pred[0] == "eq":
                        field, val = pred[1], pred[2]
                        if getattr(d, field, None) != val:
                            ok = False
                            break
                    if len(pred) == 2 and pred[0] == "search":
                        term = pred[1].strip("%").lower()
                        if term not in ((d.name or "") + (d.description or "")).lower():
                            ok = False
                            break
            if ok:
                filtered.append(d)
        return FakeQuery(filtered)

    def with_entities(self, *args, **kwargs):
        return self

    def order_by(self, *_):
        self._data.sort(key=lambda d: d.created_at, reverse=True)
        return self

    def offset(self, n: int):
        return FakeQuery(self._data[n:])

    def limit(self, n: int):
        return FakeQuery(self._data[:n])

    def all(self):
        total = len(self._data)
        return [FakeRow(d, total) for d in self._data]


class FakeGroupQuery:
    def __init__(self, data: List[FakeDataset]):
        self._data = data

    def group_by(self, *args, **kwargs):
        return self

    def all(self):
        # group by data_origin
        groups = {}
        for d in self._data:
            groups.setdefault(d.data_origin, []).append(d)
        rows = []
        for origin, items in groups.items():
            count = len(items)
            last_sync = max((x.updated_at for x in items), default=None)
            rows.append((origin, count, last_sync))
        return rows


class FakeDB:
    def __init__(self, datasets: List[FakeDataset]):
        self.datasets = datasets

    def query(self, *args, **kwargs):
        if args and args[0] is FakeDataset:
            return FakeQuery(self.datasets)
        # group by query call
        return FakeGroupQuery(self.datasets)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


@pytest.fixture(autouse=True)
def patch_db(monkeypatch):
    monkeypatch.setattr(catalog_service, "Dataset", FakeDataset)
    monkeypatch.setattr(catalog_service, "db_session", lambda: FakeDB([]))
    monkeypatch.setattr(catalog_service, "or_", lambda a, b: ("search", a[2]))
    monkeypatch.setattr(catalog_service, "selectinload", lambda *args, **kwargs: None)
    yield


def test_normalize_dt_and_health_status():
    now = datetime.now(timezone.utc)
    stale = now - timedelta(days=40)
    warn = now - timedelta(days=10)
    assert catalog_service._health_status(stale) == "stale"
    assert catalog_service._health_status(warn) == "warning"
    assert catalog_service._health_status(now) == "healthy"
    assert catalog_service._health_status(None) == "stale"


def test_get_repository_summary_groups_and_health(monkeypatch):
    now = datetime.now(timezone.utc)
    ds1 = FakeDataset(data_origin="SRC", updated_at=now)
    ds2 = FakeDataset(data_origin="SRC", updated_at=now - timedelta(days=8))
    ds3 = FakeDataset(data_origin=None, updated_at=None)
    monkeypatch.setattr(catalog_service, "db_session", lambda: FakeDB([ds1, ds2, ds3]))

    rows = catalog_service.get_repository_summary()
    names = {r["name"] for r in rows}
    assert "SRC" in names and "Unknown" in names
    src_row = next(r for r in rows if r["name"] == "SRC")
    assert src_row["dataset_count"] == 2
    assert src_row["health_status"] in {"healthy", "warning"}


def test_get_catalog_datasets_filters_and_paginates(monkeypatch):
    now = datetime.now(timezone.utc)
    d1 = FakeDataset(name="Alpha", data_origin="SRC", description="first", created_at=now, features=[1, 2], external_ids={"accession": "A1"})
    d2 = FakeDataset(name="Beta", data_origin="OTHER", description="second", created_at=now - timedelta(days=1), features=[], external_ids={})
    monkeypatch.setattr(catalog_service, "db_session", lambda: FakeDB([d1, d2]))

    rows, total = catalog_service.get_catalog_datasets(source="SRC", search="alp", limit=1, offset=0)
    assert total == 1
    assert len(rows) == 1
    assert rows[0]["accession"] == "A1"
    assert rows[0]["feature_count"] == 2

