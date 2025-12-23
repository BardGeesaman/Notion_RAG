from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Dict, Generator, List, Type
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


client = TestClient(app)


class _FakeQuery:
    def __init__(self, rows: List[Any]):
        self._rows = list(rows)

    def filter_by(self, **kwargs):
        def ok(r: Any) -> bool:
            return all(getattr(r, k) == v for k, v in kwargs.items())

        return _FakeQuery([r for r in self._rows if ok(r)])

    def filter(self, *args, **kwargs):
        # Keep it simple: ignore complex filters in unit tests.
        return self

    def order_by(self, *args, **kwargs):
        return self

    def limit(self, n: int):
        return _FakeQuery(self._rows[:n])

    def all(self):
        return list(self._rows)

    def first(self):
        return self._rows[0] if self._rows else None


class _FakeSession:
    def __init__(self, state: Dict[str, List[Any]]):
        self.state = state

    def query(self, model: Type[Any]):
        name = getattr(model, "__name__", str(model))
        return _FakeQuery(self.state.get(name, []))

    def add(self, obj: Any):
        name = obj.__class__.__name__
        self.state.setdefault(name, []).append(obj)

    def commit(self):
        return None

    def refresh(self, obj: Any):
        return None


@pytest.fixture
def docking_db_stub(monkeypatch):
    # Import models locally so we can use their names
    from amprenta_rag.database import models as m

    # Minimal objects
    structure = m.ProteinStructure(id=uuid4(), source="pdb")
    pocket = m.BindingSite(id=uuid4(), structure_id=structure.id, pocket_rank=1, center_x=1.0, center_y=2.0, center_z=3.0)
    cmp = m.Compound(id=uuid4(), compound_id="CMPD-001", smiles="CCO")

    state: Dict[str, List[Any]] = {
        "ProteinStructure": [structure],
        "BindingSite": [pocket],
        "Compound": [cmp],
        "DockingRun": [],
        "DockingPose": [],
    }

    @contextmanager
    def _fake_db_session() -> Generator[_FakeSession, None, None]:
        yield _FakeSession(state)

    monkeypatch.setattr("amprenta_rag.api.routers.docking.db_session", _fake_db_session)

    # Avoid starting background threads
    monkeypatch.setattr("amprenta_rag.api.routers.docking.DOCKING_SERVICE.start_run", lambda run_id: None)

    return {"structure": structure, "pocket": pocket, "compound": cmp, "state": state}


def test_create_run_ok(docking_db_stub):
    structure = docking_db_stub["structure"]
    pocket = docking_db_stub["pocket"]

    payload = {
        "structure_id": str(structure.id),
        "binding_site_id": str(pocket.id),
        "compound_ids": ["CMPD-001"],
    }
    resp = client.post("/api/docking/runs", json=payload)
    assert resp.status_code == 200
    data = resp.json()
    assert data["structure_id"] == str(structure.id)
    assert data["binding_site_id"] == str(pocket.id)
    assert data["status"] in ("pending", "running", "completed", "failed")
    assert data["total_compounds"] == 1


def test_list_runs_empty_then_one(docking_db_stub):
    resp = client.get("/api/docking/runs")
    assert resp.status_code == 200
    assert resp.json() == []

    structure = docking_db_stub["structure"]
    payload = {
        "structure_id": str(structure.id),
        "binding_site_id": None,
        "compound_ids": ["CMPD-001"],
    }
    _ = client.post("/api/docking/runs", json=payload)

    resp2 = client.get("/api/docking/runs")
    assert resp2.status_code == 200
    runs = resp2.json()
    assert isinstance(runs, list)
    assert len(runs) == 1


