from __future__ import annotations

from datetime import datetime
from typing import Any, List
from uuid import uuid4

import pytest

from amprenta_rag.utils import activity


class _Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):  # type: ignore[override]
        return (self.name, other)

    def desc(self):
        return self


class FakeExperiment:
    created_at = _Field("created_at")
    id = _Field("id")

    def __init__(self, name: str, created_by=None, created_at: datetime | None = None):
        self.id = uuid4()
        self.name = name
        self.created_by = created_by
        self.created_at = created_at


class FakeDiscoveredStudy:
    discovered_at = _Field("discovered_at")

    def __init__(self, study_id: str, repository: str, title: str, status: str, discovered_at: datetime | None):
        self.id = uuid4()
        self.study_id = study_id
        self.repository = repository
        self.title = title
        self.status = status
        self.discovered_at = discovered_at


class FakeCompound:
    created_at = _Field("created_at")

    def __init__(self, compound_id: str, smiles: str, created_at: datetime | None):
        self.id = uuid4()
        self.compound_id = compound_id
        self.smiles = smiles
        self.created_at = created_at


class FakeDataset:
    def __init__(self):
        self.id = uuid4()


class FakeUser:
    id = _Field("id")
    is_active = _Field("is_active")

    def __init__(self, username: str, is_active: bool = True):
        self._id_val = uuid4()
        self.username = username
        self._is_active_val = is_active

    def __getattribute__(self, name: str):
        if name == "id":
            return object.__getattribute__(self, "_id_val")
        if name == "is_active":
            return object.__getattribute__(self, "_is_active_val")
        return object.__getattribute__(self, name)


class FakeQuery:
    def __init__(self, data: List[Any]):
        self._data = list(data)
        self._limit: int | None = None

    def filter(self, *predicates: Any) -> "FakeQuery":
        filtered = list(self._data)
        for pred in predicates:
            if isinstance(pred, tuple):
                field, val = pred
                filtered = [obj for obj in filtered if getattr(obj, field, None) == val]
            elif isinstance(pred, _Field):
                filtered = [obj for obj in filtered if getattr(obj, pred.name, None)]
            elif isinstance(pred, bool):
                if not pred:
                    filtered = []
        return FakeQuery(filtered)

    def order_by(self, *args: Any) -> "FakeQuery":
        return self

    def limit(self, n: int) -> "FakeQuery":
        self._limit = n
        return self

    def all(self) -> List[Any]:
        if self._limit is None:
            return list(self._data)
        return list(self._data)[: self._limit]

    def count(self) -> int:
        return len(self.all())


class FakeSession:
    def __init__(self):
        self.experiments: List[FakeExperiment] = []
        self.discoveries: List[FakeDiscoveredStudy] = []
        self.compounds: List[FakeCompound] = []
        self.datasets: List[FakeDataset] = []
        self.users: List[FakeUser] = []

    def query(self, model: Any) -> FakeQuery:
        if model is activity.Experiment:
            return FakeQuery(self.experiments)
        if model is activity.DiscoveredStudy:
            return FakeQuery(self.discoveries)
        if model is activity.Compound:
            return FakeQuery(self.compounds)
        if model is activity.Dataset:
            return FakeQuery(self.datasets)
        if model is activity.User:
            return FakeQuery(self.users)
        return FakeQuery([])


@pytest.fixture(autouse=True)
def patch_models(monkeypatch):
    monkeypatch.setattr(activity, "Experiment", FakeExperiment)
    monkeypatch.setattr(activity, "DiscoveredStudy", FakeDiscoveredStudy)
    monkeypatch.setattr(activity, "Compound", FakeCompound)
    monkeypatch.setattr(activity, "Dataset", FakeDataset)
    monkeypatch.setattr(activity, "User", FakeUser)
    yield


def test_get_recent_experiments():
    db = FakeSession()
    user = FakeUser(username="alice")
    db.users.append(user)
    db.experiments.extend(
        [
            FakeExperiment("exp1", created_by=user, created_at=datetime(2024, 1, 1)),
            FakeExperiment("exp2", created_by=None, created_at=None),
        ]
    )

    results = activity.get_recent_experiments(db, limit=2)
    assert len(results) == 2
    assert results[0]["name"] == "exp1"
    assert results[0]["created_by"] == "alice"
    assert results[1]["created_by"] is None


def test_get_recent_discoveries():
    db = FakeSession()
    db.discoveries.append(
        FakeDiscoveredStudy("ST001", "mw", "title", "new", datetime(2024, 2, 2))
    )

    results = activity.get_recent_discoveries(db)
    assert results[0]["study_id"] == "ST001"
    assert results[0]["status"] == "new"


def test_get_recent_compounds_truncates_smiles():
    db = FakeSession()
    long_smiles = "C" * 60
    db.compounds.append(FakeCompound("cmpd1", long_smiles, datetime(2024, 3, 3)))

    results = activity.get_recent_compounds(db)
    assert results[0]["compound_id"] == "cmpd1"
    assert results[0]["smiles"].endswith("...")
    assert len(results[0]["smiles"]) < len(long_smiles)


def test_get_activity_stats_counts_active_users_only():
    db = FakeSession()
    db.experiments.extend([FakeExperiment("e1"), FakeExperiment("e2")])
    db.compounds.append(FakeCompound("c1", "CC", None))
    db.datasets.append(FakeDataset())
    db.users.extend([FakeUser("active", True), FakeUser("inactive", False)])

    stats = activity.get_activity_stats(db)
    assert stats == {
        "total_experiments": 2,
        "total_compounds": 1,
        "total_datasets": 1,
        "total_users": 1,
    }

