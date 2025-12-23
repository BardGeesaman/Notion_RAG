from __future__ import annotations

from uuid import uuid4


from amprenta_rag.chemistry import candidate_selection as cs


def test_get_traffic_light_ranges():
    assert cs.get_traffic_light(5, 0, 10) == "green"
    assert cs.get_traffic_light(11, 0, 10) == "yellow"  # within tolerance
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
        return obj


class _FakeCompound:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


class _FakeTPP:
    def __init__(self, criteria):
        self.criteria = criteria


def test_score_compound_against_tpp_happy():
    cmpd = _FakeCompound(molecular_weight=300.0)
    tpp = _FakeTPP(criteria=[{"property": "molecular_weight", "min": 250, "max": 350, "weight": 2.0, "unit": "Da"}])
    db = _FakeDB(compound=cmpd, tpp=tpp)

    res = cs.score_compound_against_tpp(uuid4(), uuid4(), db)
    assert res["overall_score"] == 100.0
    assert res["criteria_scores"][0]["traffic_light"] == "green"


def test_score_compound_missing_property_returns_zero():
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

