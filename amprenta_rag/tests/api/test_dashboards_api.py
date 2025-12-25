from __future__ import annotations

import uuid

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import text

from amprenta_rag.api.main import app
from amprenta_rag.database.models import Program, User
from amprenta_rag.database.session import db_session
from amprenta_rag.notebooks.registry import load_registry


client = TestClient(app)


def _db_available() -> bool:
    try:
        with db_session() as db:
            db.execute(text("SELECT 1"))
        return True
    except Exception:
        return False


def _table_exists(table: str) -> bool:
    try:
        with db_session() as db:
            ok = db.execute(
                text(
                    """
SELECT 1
FROM information_schema.tables
WHERE table_schema='public' AND table_name=:t
"""
                ),
                {"t": table},
            ).fetchone()
        return bool(ok)
    except Exception:
        return False


pytestmark = pytest.mark.requires_postgres


def _any_notebook_path_or_skip() -> str:
    try:
        reg = load_registry()
    except Exception:
        pytest.skip("Notebook registry.json not available")
    for it in reg:
        if isinstance(it, dict) and it.get("notebook_path"):
            return str(it["notebook_path"])
    pytest.skip("No notebook_path entries in registry.json")


def test_pin_requires_auth():
    if not _db_available():
        pytest.skip("Postgres not available/configured for dashboards API test")
    if not _table_exists("pinned_dashboards"):
        pytest.skip("pinned_dashboards table not present (migration not applied)")

    nb_path = _any_notebook_path_or_skip()
    resp = client.post(
        "/api/dashboards/pin",
        json={"notebook_path": nb_path, "program_id": str(uuid.uuid4())},
    )
    assert resp.status_code == 401


def test_pin_dashboard_to_program_and_lists_and_unpin():
    if not _db_available():
        pytest.skip("Postgres not available/configured for dashboards API test")
    for t in ("users", "programs", "pinned_dashboards"):
        if not _table_exists(t):
            pytest.skip(f"{t} table not present (migrations not applied)")

    nb_path = _any_notebook_path_or_skip()

    with db_session() as db:
        user = User(
            username=f"user_{uuid.uuid4().hex[:6]}",
            email=f"user_{uuid.uuid4().hex[:6]}@x.com",
            password_hash="x",
            role="researcher",
        )
        prog = Program(name=f"Program {uuid.uuid4().hex[:6]}", description="x")
        db.add(user)
        db.add(prog)
        db.commit()
        db.refresh(user)
        db.refresh(prog)
        uid = str(user.id)
        pid = str(prog.id)

    headers = {"X-User-Id": uid}

    # Pin
    resp = client.post(
        "/api/dashboards/pin",
        json={"notebook_path": nb_path, "program_id": pid, "display_name": "Demo Dashboard", "config": {"k": "v"}},
        headers=headers,
    )
    assert resp.status_code == 201, resp.text
    pin = resp.json()
    assert pin["program_id"] == pid
    assert pin["notebook_path"] == nb_path
    pin_id = pin["id"]

    # List by program
    resp2 = client.get(f"/api/dashboards/program/{pid}", headers=headers)
    assert resp2.status_code == 200
    items = resp2.json()
    assert any(x["id"] == pin_id for x in items)

    # List mine
    resp3 = client.get("/api/dashboards/mine", headers=headers)
    assert resp3.status_code == 200
    mine = resp3.json()
    assert any(x["id"] == pin_id for x in mine)

    # Unpin
    resp4 = client.delete(f"/api/dashboards/pin/{pin_id}", headers=headers)
    assert resp4.status_code == 200
    assert resp4.json()["unpinned"] is True


