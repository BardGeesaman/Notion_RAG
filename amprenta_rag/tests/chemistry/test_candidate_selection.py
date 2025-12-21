from __future__ import annotations

from uuid import uuid4

import pytest

from amprenta_rag.chemistry import candidate_selection as cs


def test_get_traffic_light_ranges():
    assert cs.get_traffic_light(5, 0, 10) == "green"
    assert cs.get_traffic_light(11, 0, 10) == "yellow"  # within 20%
    assert cs.get_traffic_light(20, 0, 10) == "red"
    assert cs.get_traffic_light(5, 5, 5) == "green"


class _FakeQuery:
    def __init__(self, obj):
        self.obj = obj

    def filter(self, *args, **kwargs):
        return self

    def first(self):
        return self.obj


class _FakeDB:
    def __init__(self, compound=None, tpp=None):
        self.compound = compound
        self.tpp = tpp
        self.added = []
        self.committed = False
        self.refreshed = []

    def query(self, model):
        if model is cs.Compound:
            return _FakeQuery(self.compound)
        if model is cs.TargetProductProfile:
            return _FakeQuery(self.tpp)
        return _FakeQuery(None)

    def add(self, obj):
        self.added.append(obj)

    def commit(self):
        self.committed = True

    def refresh(self, obj):
        self.refreshed.append(obj)


class _FakeCompound:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


class _FakeTPP:
    def __init__(self, criteria):
        self.criteria = criteria


def test_score_compound_against_tpp_happy(monkeypatch):
    cmpd = _FakeCompound(molecular_weight=300.0)
    tpp = _FakeTPP(criteria=[{"property": "molecular_weight", "min": 250, "max": 350, "weight": 2.0, "unit": "Da"}])
    db = _FakeDB(compound=cmpd, tpp=tpp)

    res = cs.score_compound_against_tpp(uuid4(), uuid4(), db)
    assert res["overall_score"] == 100.0
    assert res["criteria_scores"][0]["traffic_light"] == "green"


def test_score_compound_missing_property_returns_zero(monkeypatch):
    cmpd = _FakeCompound()
    tpp = _FakeTPP(criteria=[{"property": "unknown", "min": 0, "max": 1}])
    db = _FakeDB(compound=cmpd, tpp=tpp)

    res = cs.score_compound_against_tpp(uuid4(), uuid4(), db)
    assert res["overall_score"] == 0.0
    assert res["criteria_scores"] == []


def test_nominate_candidate_creates_record(monkeypatch):
    fake_score = {
        "criteria_scores": [{"property": "mw", "value": 1, "score": 100, "traffic_light": "green"}],
        "overall_score": 90.0,
    }
    monkeypatch.setattr(cs, "score_compound_against_tpp", lambda compound_id, tpp_id, db: fake_score)
    monkeypatch.setattr(cs, "ensure_uuid", lambda x: x)

    created = {}

    class FakeNomination:
        def __init__(self, **kwargs):
            created.update(kwargs)

    monkeypatch.setattr(cs, "CandidateNomination", FakeNomination)

    db = _FakeDB()
    nomination = cs.nominate_candidate(uuid4(), uuid4(), uuid4(), "note", db)

    assert isinstance(nomination, FakeNomination)
    assert created["scores"]["mw"]["score"] == 100
    assert db.committed is True
    assert nomination is db.added[0]
from __future__ import annotations

from uuid import uuid4

import pytest

from amprenta_rag.chemistry import candidate_selection as cs


class FakeColumn:
    def __init__(self, attr: str):
        self.attr = attr

    def __eq__(self, other):
        return lambda obj: getattr(obj, self.attr) == other


class FakeCompound:
    id = FakeColumn("id")

    def __init__(self, id_, molecular_weight=300.0, logp=2.0, hbd_count=1, hba_count=2, rotatable_bonds=3, tpsa=10.0):
        self.id = id_
        self.molecular_weight = molecular_weight
        self.logp = logp
        self.hbd_count = hbd_count
        self.hba_count = hba_count
        self.rotatable_bonds = rotatable_bonds
        self.tpsa = tpsa


class FakeTPP:
    id = FakeColumn("id")

    def __init__(self, id_, criteria):
        self.id = id_
        self.criteria = criteria


class FakeNomination:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def filter(self, cond):
        self._data = [obj for obj in self._data if cond(obj)]
        return self

    def first(self):
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self, compounds=None, tpps=None):
        self.compounds = compounds or []
        self.tpps = tpps or []
        self.added = []
        self.commits = 0
        self.refreshed = []

    def query(self, model):
        if model is cs.Compound:
            return FakeQuery(self.compounds)
        if model is cs.TargetProductProfile:
            return FakeQuery(self.tpps)
        return FakeQuery([])

    def add(self, obj):
        self.added.append(obj)

    def commit(self):
        self.commits += 1

    def refresh(self, obj):
        self.refreshed.append(obj)


def test_get_traffic_light_boundaries():
    assert cs.get_traffic_light(5, 0, 10) == "green"
    assert cs.get_traffic_light(11, 0, 10) == "yellow"  # within 20%
    assert cs.get_traffic_light(20, 0, 10) == "red"
    assert cs.get_traffic_light(5, 5, 5) == "green"
    assert cs.get_traffic_light(6, 5, 5) == "red"


def test_score_compound_against_tpp_calculates_scores():
    cid = uuid4()
    tid = uuid4()
    compound = FakeCompound(cid, molecular_weight=300, logp=2.0)
    tpp = FakeTPP(
        tid,
        criteria=[
            {"property": "molecular_weight", "min": 250, "max": 350, "weight": 1},
            {"property": "logp", "min": 1, "max": 3, "weight": 1},
        ],
    )
    session = FakeSession(compounds=[compound], tpps=[tpp])

    result = cs.score_compound_against_tpp(cid, tid, session)
    assert result["compound_id"] == str(cid)
    assert result["overall_score"] == 100.0  # both green
    assert len(result["criteria_scores"]) == 2


def test_score_compound_missing_raises():
    session = FakeSession(compounds=[], tpps=[])
    with pytest.raises(ValueError):
        cs.score_compound_against_tpp(uuid4(), uuid4(), session)


def test_nominate_candidate_builds_scores(monkeypatch):
    session = FakeSession()
    cid = uuid4()
    tid = uuid4()
    uid = uuid4()

    fake_scores = {
        "compound_id": str(cid),
        "tpp_id": str(tid),
        "overall_score": 80.0,
        "criteria_scores": [
            {"property": "mw", "value": 300, "score": 100, "traffic_light": "green"},
            {"property": "logp", "value": 5, "score": 0, "traffic_light": "red"},
        ],
    }

    monkeypatch.setattr(cs, "score_compound_against_tpp", lambda compound_id, tpp_id, db: fake_scores)

    class FakeCandidateNomination:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    monkeypatch.setattr(cs, "CandidateNomination", FakeCandidateNomination)

    nomination = cs.nominate_candidate(cid, tid, uid, "note", session)
    assert nomination.overall_score == 80.0
    assert "mw" in nomination.scores
    assert nomination.nominated_by_id is not None
    assert session.commits == 1

