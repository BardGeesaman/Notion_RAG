from __future__ import annotations

from datetime import datetime, timedelta
from uuid import uuid4

import pytest

from amprenta_rag.utils import widgets


class _Field:
    def __init__(self, name: str):
        self.name = name

    def __ge__(self, other):
        return (self.name, other, "ge")


class _Query:
    def __init__(self, data):
        self._data = list(data)

    def count(self):
        return len(self._data)

    def filter(self, predicate):
        if isinstance(predicate, tuple) and len(predicate) == 3 and predicate[2] == "ge":
            name, cutoff, _ = predicate
            filtered = [obj for obj in self._data if getattr(obj, name) >= cutoff]
            return _Query(filtered)
        return _Query(self._data)


class FakeSession:
    def __init__(self, mapping):
        self._mapping = mapping

    def query(self, model):
        return _Query(self._mapping.get(model, []))


class FakeExperiment:
    pass


class FakeCompound:
    pass


class FakeSample:
    pass


class FakeDiscovery:
    discovered_at = _Field("discovered_at")

    def __init__(self, discovered_at):
        self.discovered_at = discovered_at
        self.id = uuid4()


@pytest.fixture(autouse=True)
def patch_models(monkeypatch):
    monkeypatch.setattr(widgets, "Experiment", FakeExperiment)
    monkeypatch.setattr(widgets, "Compound", FakeCompound)
    monkeypatch.setattr(widgets, "Sample", FakeSample)
    monkeypatch.setattr(widgets, "DiscoveredStudy", FakeDiscovery)
    yield


def test_counts_for_experiment_compound_sample():
    db = FakeSession(
        {
            FakeExperiment: [1, 2, 3],
            FakeCompound: [1, 2],
            FakeSample: [1],
        }
    )
    assert widgets.get_experiment_count(db) == 3
    assert widgets.get_compound_count(db) == 2
    assert widgets.get_sample_count(db) == 1


def test_discovery_count_within_days():
    now = datetime.utcnow()
    recent = FakeDiscovery(now - timedelta(days=1))
    old = FakeDiscovery(now - timedelta(days=10))
    db = FakeSession({FakeDiscovery: [recent, old]})

    assert widgets.get_discovery_count(days=7, db=db) == 1
    assert widgets.get_discovery_count(days=30, db=db) == 2


def test_discovery_count_db_none_returns_zero():
    assert widgets.get_discovery_count(db=None) == 0

