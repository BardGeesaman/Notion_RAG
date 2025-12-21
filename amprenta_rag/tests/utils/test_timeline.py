from __future__ import annotations

from datetime import datetime, timedelta
from uuid import uuid4

import pytest

from amprenta_rag.utils import timeline


class _Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return (self.name, other)

    def desc(self):
        return self


class FakeUser:
    def __init__(self, username: str):
        self.username = username


class FakeExperiment:
    created_at = _Field("created_at")
    created_by_id = _Field("created_by_id")

    def __init__(self, id_, name, created_at, created_by, created_by_id):
        self.id = id_
        self.name = name
        self.created_at = created_at
        self.created_by = created_by
        self.created_by_id = created_by_id


class FakeCompound:
    created_at = _Field("created_at")

    def __init__(self, id_, compound_id, created_at):
        self.id = id_
        self.compound_id = compound_id
        self.created_at = created_at


class FakeSample:
    created_at = _Field("created_at")
    created_by_id = _Field("created_by_id")

    def __init__(self, id_, name, created_at, created_by, created_by_id):
        self.id = id_
        self.name = name
        self.created_at = created_at
        self.created_by = created_by
        self.created_by_id = created_by_id


class FakeDiscoveredStudy:
    discovered_at = _Field("discovered_at")

    def __init__(self, id_, study_id, discovered_at):
        self.id = id_
        self.study_id = study_id
        self.discovered_at = discovered_at


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def order_by(self, *_):
        return self

    def filter(self, *predicates):
        filtered = []
        for obj in self._data:
            ok = True
            for pred in predicates:
                if isinstance(pred, tuple):
                    field, val = pred
                    actual = getattr(obj, field, None)
                    if actual != val:
                        ok = False
                        break
            if ok:
                filtered.append(obj)
        return FakeQuery(filtered)

    def limit(self, n):
        return FakeQuery(self._data[:n])

    def all(self):
        return list(self._data)


class FakeSession:
    def __init__(self, experiments=None, compounds=None, samples=None, discoveries=None):
        self.store = {
            FakeExperiment: experiments or [],
            FakeCompound: compounds or [],
            FakeSample: samples or [],
            FakeDiscoveredStudy: discoveries or [],
        }

    def query(self, model):
        return FakeQuery(self.store.get(model, []))


@pytest.fixture(autouse=True)
def patch_models(monkeypatch):
    monkeypatch.setattr(timeline, "Experiment", FakeExperiment)
    monkeypatch.setattr(timeline, "Compound", FakeCompound)
    monkeypatch.setattr(timeline, "Sample", FakeSample)
    monkeypatch.setattr(timeline, "DiscoveredStudy", FakeDiscoveredStudy)
    yield


def test_get_timeline_combines_and_sorts():
    now = datetime.utcnow()
    user1 = uuid4()
    exp = FakeExperiment(uuid4(), "Exp1", now - timedelta(minutes=1), FakeUser("alice"), user1)
    comp = FakeCompound(uuid4(), "C001", now)
    sample = FakeSample(uuid4(), "S1", now - timedelta(minutes=2), FakeUser("bob"), uuid4())
    disc = FakeDiscoveredStudy(uuid4(), "DISC", now - timedelta(minutes=3))
    db = FakeSession(experiments=[exp], compounds=[comp], samples=[sample], discoveries=[disc])

    items = timeline.get_timeline(limit=10, db=db)

    assert [i["type"] for i in items] == ["compound", "experiment", "sample", "discovery"]
    assert items[0]["name"] == "C001"
    assert items[1]["user"] == "alice"
    assert items[2]["user"] == "bob"


def test_user_filter_applies_to_experiments_and_samples():
    now = datetime.utcnow()
    user1 = uuid4()
    user2 = uuid4()
    exp1 = FakeExperiment(uuid4(), "E1", now, FakeUser("alice"), user1)
    exp2 = FakeExperiment(uuid4(), "E2", now - timedelta(minutes=1), FakeUser("carol"), user2)
    sample1 = FakeSample(uuid4(), "S1", now - timedelta(minutes=2), FakeUser("bob"), user1)
    sample2 = FakeSample(uuid4(), "S2", now - timedelta(minutes=3), FakeUser("dave"), user2)
    comp = FakeCompound(uuid4(), "C1", now - timedelta(minutes=4))
    db = FakeSession(experiments=[exp1, exp2], compounds=[comp], samples=[sample1, sample2], discoveries=[])

    items = timeline.get_timeline(user_id=str(user1), db=db)

    types = [i["type"] for i in items]
    assert "experiment" in types and "sample" in types
    assert all(i["user"] in {"alice", "bob", None} for i in items)
    assert all(i["type"] != "experiment" or i["user"] == "alice" for i in items if i["type"] == "experiment")
    assert all(i["type"] != "sample" or i["user"] == "bob" for i in items if i["type"] == "sample")


def test_limit_respected_after_sort():
    now = datetime.utcnow()
    comp_new = FakeCompound(uuid4(), "C-new", now)
    comp_old = FakeCompound(uuid4(), "C-old", now - timedelta(days=1))
    db = FakeSession(experiments=[], compounds=[comp_new, comp_old], samples=[], discoveries=[])

    items = timeline.get_timeline(limit=1, db=db)
    assert len(items) == 1
    assert items[0]["name"] == "C-new"


def test_handles_none_timestamp():
    now = datetime.utcnow()
    comp_new = FakeCompound(uuid4(), "C-new", now)
    sample_none = FakeSample(uuid4(), "S-none", None, FakeUser("x"), uuid4())
    db = FakeSession(experiments=[], compounds=[comp_new], samples=[sample_none], discoveries=[])

    items = timeline.get_timeline(limit=5, db=db)
    assert len(items) == 2
    # None timestamp item should not break sort
    assert items[0]["name"] == "C-new"

