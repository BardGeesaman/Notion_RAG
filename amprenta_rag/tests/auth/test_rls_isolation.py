from __future__ import annotations

import uuid

import pytest
from sqlalchemy import text

from amprenta_rag.database.models import Company, User
from amprenta_rag.database.session import db_session


def _db_available() -> bool:
    try:
        with db_session() as db:
            db.execute(text("SELECT 1"))
        return True
    except Exception:
        return False


pytestmark = pytest.mark.requires_postgres


def test_rls_isolation_blocks_cross_company_reads():
    if not _db_available():
        pytest.skip("Postgres not available/configured for RLS test")

    with db_session() as db:
        # Preconditions: datasets table must have company_id and RLS enabled
        cols = db.execute(
            text(
                """
SELECT 1
FROM information_schema.columns
WHERE table_schema='public' AND table_name='datasets' AND column_name='company_id'
"""
            )
        ).fetchall()
        if not cols:
            pytest.skip("datasets.company_id column not present (tenant migrations not applied)")

        rls = db.execute(text("SELECT relrowsecurity FROM pg_class WHERE relname='datasets'")).scalar()
        if not rls:
            pytest.skip("RLS is not enabled on datasets table")

        # Create 2 companies
        c1 = Company(name="A", subdomain=f"a{uuid.uuid4().hex[:8]}")
        c2 = Company(name="B", subdomain=f"b{uuid.uuid4().hex[:8]}")
        db.add(c1)
        db.add(c2)
        db.commit()
        db.refresh(c1)
        db.refresh(c2)

        # Create 2 users, each in a different company
        u1 = User(username=f"u1_{uuid.uuid4().hex[:6]}", email=f"u1_{uuid.uuid4().hex[:6]}@x.com", password_hash="x", role="admin")
        u2 = User(username=f"u2_{uuid.uuid4().hex[:6]}", email=f"u2_{uuid.uuid4().hex[:6]}@x.com", password_hash="x", role="admin")
        u1.company_id = c1.id
        u2.company_id = c2.id
        db.add(u1)
        db.add(u2)
        db.commit()
        db.refresh(u1)
        db.refresh(u2)

        # Insert a dataset row owned by company A (raw SQL because Dataset ORM may not yet have company_id)
        ds_id = uuid.uuid4()
        db.execute(
            text(
                """
INSERT INTO datasets (id, name, omics_type, description, created_at, updated_at, created_by_id, ingestion_status, company_id)
VALUES (:id, :name, :omics_type, :description, now(), now(), :created_by_id, 'complete', :company_id)
"""
            ),
            {
                "id": str(ds_id),
                "name": "RLS Test Dataset",
                "omics_type": "transcriptomics",
                "description": "rls",
                "created_by_id": str(u1.id),
                "company_id": str(c1.id),
            },
        )
        db.commit()

        # As company A, should see it
        db.execute(text(f"SET LOCAL app.current_company_id = '{c1.id}'"))
        n_a = db.execute(text("SELECT count(*) FROM datasets WHERE id = :id"), {"id": str(ds_id)}).scalar()
        assert int(n_a) == 1

        # As company B, should not see it
        db.execute(text(f"SET LOCAL app.current_company_id = '{c2.id}'"))
        n_b = db.execute(text("SELECT count(*) FROM datasets WHERE id = :id"), {"id": str(ds_id)}).scalar()
        assert int(n_b) == 0


