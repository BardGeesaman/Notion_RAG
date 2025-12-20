 # mypy: ignore-errors
from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4
from unittest.mock import Mock

import pytest
from fastapi.testclient import TestClient


@pytest.mark.unit
def test_diff_protocols_detects_added_removed_changed_steps_materials_and_params():
    from amprenta_rag.analysis import protocol_diff as pd

    id1 = uuid4()
    id2 = uuid4()

    v1 = SimpleNamespace(
        id=id1,
        steps=[{"order": 1, "name": "A"}, {"order": 2, "name": "B"}],
        materials=[{"name": "NaCl", "amount": "1g"}],
        parameters={"temp": 25, "time": 10},
    )
    v2 = SimpleNamespace(
        id=id2,
        steps=[{"order": 1, "name": "A"}, {"order": 2, "name": "B2"}, {"order": 3, "name": "C"}],
        materials=[{"name": "NaCl", "amount": "2g"}, {"name": "H2O"}],
        parameters={"temp": 30, "time": 10, "ph": 7.4},
    )

    diff = pd.diff_protocols(v1, v2)

    assert diff.protocol_id == id1
    assert diff.other_id == id2

    # Steps: added order 3, changed order 2, none removed
    assert diff.added_steps == [{"order": 3, "name": "C"}]
    assert diff.removed_steps == []
    assert diff.changed_steps == [{"from": {"order": 2, "name": "B"}, "to": {"order": 2, "name": "B2"}}]

    # Materials: keyed by order/id/name; only added/removed keys are tracked (no "changed materials" concept yet)
    # So NaCl amount change is not reported; only H2O is added.
    assert diff.materials_added == [{"name": "H2O"}]
    assert diff.materials_removed == []

    # Parameters: temp changed, ph added
    assert {"name": "temp", "from": 25, "to": 30} in diff.parameters_changed
    assert {"name": "ph", "from": None, "to": 7.4} in diff.parameters_changed
    assert all(item["name"] != "time" for item in diff.parameters_changed)


@pytest.mark.unit
def test_get_protocol_history_returns_version_chain_oldest_to_newest(monkeypatch):
    from amprenta_rag.analysis import protocol_diff as pd

    root_id = uuid4()
    mid_id = uuid4()
    leaf_id = uuid4()

    root = SimpleNamespace(id=root_id, version=1, parent_id=None, name="root")
    mid = SimpleNamespace(id=mid_id, version=2, parent_id=root_id, name="mid")
    leaf = SimpleNamespace(id=leaf_id, version=3, parent_id=mid_id, name="leaf")

    by_id = {root_id: root, mid_id: mid, leaf_id: leaf}

    class FakeQuery:
        def __init__(self):
            self._id = None

        def filter(self, expr):
            # SQLAlchemy expr comparison isn't evaluated; extract via string fallback isn't stable.
            # We rely on the test to set _id via a closure around protocol_id.
            return self

        def first(self):
            return by_id.get(self._id)

    class FakeDB:
        def __init__(self):
            self._query = FakeQuery()

        def query(self, _model):
            return self._query

    db = FakeDB()

    # Patch get_protocol_history's internal db.query(...).filter(...).first() to return
    # a chain by setting FakeQuery._id before each lookup.
    orig_query = db.query

    def query_with_id(_model):
        q = orig_query(_model)
        return q

    db.query = query_with_id  # type: ignore[assignment]

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(pd, "db_session", fake_db_session)

    def first_for_id(id_):
        db._query._id = id_
        return db._query.first()

    # Monkeypatch inside function by patching Protocol lookup via db.query().filter().first()
    # We can't intercept filter arg, so we override FakeQuery.first dynamically each call by setting _id.
    # Implement by wrapping pd.db_session consumers: set leaf first, then in loop set parent.
    # Achieve by patching pd.db_session to a generator that updates db._query._id from a stack.
    ids = [leaf_id, mid_id, root_id]

    @contextmanager
    def fake_db_session_seq():
        # first call should return leaf_id, subsequent parent lookups follow ids list
        original_first = db._query.first
        idx = {"i": 0}

        def seq_first():
            i = idx["i"]
            if i < len(ids):
                db._query._id = ids[i]
                idx["i"] = i + 1
            return original_first()

        db._query.first = seq_first  # type: ignore[assignment]
        try:
            yield db
        finally:
            db._query.first = original_first  # type: ignore[assignment]

    monkeypatch.setattr(pd, "db_session", fake_db_session_seq)

    history = pd.get_protocol_history(leaf_id)
    assert [h.protocol_id for h in history] == [root_id, mid_id, leaf_id]
    assert [h.version for h in history] == [1, 2, 3]
    assert [h.parent_id for h in history] == [None, root_id, mid_id]


@pytest.mark.unit
def test_audit_deviations_returns_deviation_reports(monkeypatch):
    from amprenta_rag.analysis import protocol_diff as pd

    exp_id = uuid4()
    p_id = uuid4()

    proto = SimpleNamespace(id=p_id, name="Proto", version=2)
    link1 = SimpleNamespace(experiment_id=exp_id, protocol_id=p_id, protocol=proto, deviations=[{"step": 1}])
    link2 = SimpleNamespace(experiment_id=exp_id, protocol_id=p_id, protocol=None, deviations={"note": "x"})

    db = SimpleNamespace(
        query=lambda model: SimpleNamespace(
            filter=lambda *args, **kwargs: SimpleNamespace(all=lambda: [link1, link2])
        )
    )

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(pd, "db_session", fake_db_session)

    reports = pd.audit_deviations(exp_id)
    assert len(reports) == 2

    r1 = reports[0]
    assert r1.experiment_id == exp_id
    assert r1.protocol_id == p_id
    assert r1.protocol_name == "Proto"
    assert r1.protocol_version == 2
    assert r1.deviations == [{"step": 1}]

    r2 = reports[1]
    assert r2.protocol_id == p_id
    assert r2.protocol_name == str(p_id)
    assert r2.protocol_version == 1
    assert r2.deviations == [{"note": "x"}]


@pytest.fixture()
def api_client():
    from amprenta_rag.api.main import app

    return TestClient(app)


@pytest.mark.api
def test_get_protocol_history_endpoint(api_client, monkeypatch):
    from amprenta_rag.api.routers import protocols as proto_router
    from amprenta_rag.analysis.protocol_diff import ProtocolHistoryItem

    pid = uuid4()
    items = [
        ProtocolHistoryItem(protocol_id=uuid4(), version=1, parent_id=None, name="root"),
        ProtocolHistoryItem(protocol_id=pid, version=2, parent_id=uuid4(), name="leaf"),
    ]
    monkeypatch.setattr(proto_router, "get_protocol_history", Mock(return_value=items))

    resp = api_client.get(f"/api/v1/protocols/{pid}/history")
    assert resp.status_code == 200
    body = resp.json()
    assert body[0]["name"] == "root"
    assert body[1]["protocol_id"] == str(pid)


@pytest.mark.api
def test_get_protocol_diff_endpoint_returns_schema(api_client, monkeypatch):
    from amprenta_rag.api.routers import protocols as proto_router
    from amprenta_rag.analysis.protocol_diff import ProtocolDiff

    p1 = uuid4()
    p2 = uuid4()

    # Patch db_session + protocol lookup to avoid DB
    class FakeQ:
        def __init__(self, obj):
            self._obj = obj

        def filter(self, *args, **kwargs):
            return self

        def first(self):
            return self._obj

    class FakeDB:
        def __init__(self, mapping):
            self._mapping = mapping
            self._last = None

        def query(self, _model):
            # Return a query that will be replaced after filter is called; easier:
            return self

        def filter(self, _expr):
            return self

        def first(self):
            # This FakeDB is used only after we set _last via closure
            return self._last

    proto1 = SimpleNamespace(id=p1, steps=[], materials=[], parameters={})
    proto2 = SimpleNamespace(id=p2, steps=[], materials=[], parameters={})
    db = FakeDB({})

    # Implement query() returning per-call query objects keyed by expected order (p1 then p2)
    calls = {"n": 0}

    def query(_model):
        calls["n"] += 1
        return FakeQ(proto1 if calls["n"] == 1 else proto2)

    db.query = query  # type: ignore[assignment]

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(proto_router, "db_session", fake_db_session)

    diff_obj = ProtocolDiff(
        protocol_id=p1,
        other_id=p2,
        added_steps=[{"order": 1}],
        removed_steps=[],
        changed_steps=[],
        materials_added=[],
        materials_removed=[],
        parameters_changed=[],
    )
    monkeypatch.setattr(proto_router, "diff_protocols", Mock(return_value=diff_obj))

    resp = api_client.get(f"/api/v1/protocols/{p1}/diff/{p2}")
    assert resp.status_code == 200
    assert resp.json()["protocol_id"] == str(p1)
    assert resp.json()["other_id"] == str(p2)
    assert resp.json()["added_steps"] == [{"order": 1}]


@pytest.mark.api
def test_get_protocol_diff_endpoint_404_when_missing(api_client, monkeypatch):
    from amprenta_rag.api.routers import protocols as proto_router

    p1 = uuid4()
    p2 = uuid4()

    class FakeQ:
        def __init__(self, obj):
            self._obj = obj

        def filter(self, *args, **kwargs):
            return self

        def first(self):
            return self._obj

    calls = {"n": 0}

    def query(_model):
        calls["n"] += 1
        # First protocol exists, second missing
        return FakeQ(SimpleNamespace(id=p1) if calls["n"] == 1 else None)

    db = SimpleNamespace(query=query)

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(proto_router, "db_session", fake_db_session)

    resp = api_client.get(f"/api/v1/protocols/{p1}/diff/{p2}")
    assert resp.status_code == 404
    assert resp.json()["detail"] == "Protocol not found"


@pytest.mark.api
def test_get_experiment_deviations_endpoint(api_client, monkeypatch):
    from amprenta_rag.api.routers import protocols as proto_router
    from amprenta_rag.analysis.protocol_diff import DeviationReport

    exp_id = uuid4()
    p_id = uuid4()

    monkeypatch.setattr(
        proto_router,
        "audit_deviations",
        Mock(
            return_value=[
                DeviationReport(
                    experiment_id=exp_id,
                    protocol_id=p_id,
                    protocol_name="Proto",
                    protocol_version=1,
                    deviations=[{"step": 1}],
                )
            ]
        ),
    )

    resp = api_client.get(f"/api/v1/experiments/{exp_id}/deviations")
    assert resp.status_code == 200
    assert resp.json()[0]["protocol_id"] == str(p_id)
    assert resp.json()[0]["deviations"] == [{"step": 1}]


