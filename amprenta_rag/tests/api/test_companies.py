from __future__ import annotations

import uuid

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import text

from amprenta_rag.api.main import app
from amprenta_rag.database.models import Company, User
from amprenta_rag.database.session import db_session


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


def test_companies_crud_and_user_invite_remove():
    if not _db_available():
        pytest.skip("Postgres not available/configured for companies API test")
    if not _table_exists("companies"):
        pytest.skip("companies table not present (company migrations not applied)")

    # Create superadmin user
    with db_session() as db:
        admin = User(
            username=f"admin_{uuid.uuid4().hex[:6]}",
            email=f"admin_{uuid.uuid4().hex[:6]}@x.com",
            password_hash="x",
            role="admin",
        )
        # Assign to a default company context (not required but enables set_company_context)
        default_company = Company(name="Default Org", subdomain=f"default{uuid.uuid4().hex[:6]}")
        db.add(default_company)
        db.commit()
        db.refresh(default_company)
        admin.company_id = default_company.id
        admin.company_role = "owner"
        db.add(admin)
        db.commit()
        db.refresh(admin)
        admin_id = str(admin.id)

    headers = {"X-User-Id": admin_id}

    # Create company
    sub = f"acme{uuid.uuid4().hex[:6]}"
    resp = client.post("/api/companies", json={"name": "Acme", "subdomain": sub}, headers=headers)
    assert resp.status_code == 200
    comp = resp.json()
    cid = comp["id"]

    # Get company
    resp2 = client.get(f"/api/companies/{cid}", headers=headers)
    assert resp2.status_code == 200
    assert resp2.json()["subdomain"] == sub

    # Patch company
    resp3 = client.patch(f"/api/companies/{cid}", json={"primary_color": "#2E86AB"}, headers=headers)
    assert resp3.status_code == 200
    assert resp3.json()["primary_color"] == "#2E86AB"

    # Invite user (creates if missing)
    email = f"user_{uuid.uuid4().hex[:6]}@x.com"
    resp4 = client.post(f"/api/companies/{cid}/users/invite", json={"email": email, "role": "member"}, headers=headers)
    assert resp4.status_code == 200
    data = resp4.json()
    assert "user" in data
    invited = data["user"]
    assert invited["email"] == email
    invited_id = invited["id"]

    # List company users includes invited
    resp5 = client.get(f"/api/companies/{cid}/users", headers=headers)
    assert resp5.status_code == 200
    emails = {u["email"] for u in resp5.json()}
    assert email in emails

    # Update role
    resp6 = client.patch(f"/api/companies/{cid}/users/{invited_id}", json={"role": "admin"}, headers=headers)
    assert resp6.status_code == 200
    assert resp6.json()["company_role"] == "admin"

    # Remove user
    resp7 = client.delete(f"/api/companies/{cid}/users/{invited_id}", headers=headers)
    assert resp7.status_code == 200
    assert resp7.json()["removed"] is True


