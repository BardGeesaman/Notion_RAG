"""Enable row-level security (RLS) tenant isolation policies.

Revision ID: g7h8i9j0k1l2
Revises: f6g7h8i9j0k1
Create Date: 2025-12-23
"""

from __future__ import annotations

from alembic import op


revision = "g7h8i9j0k1l2"
down_revision = "f6g7h8i9j0k1"
branch_labels = None
depends_on = None


TENANT_TABLES = [
    "datasets",
    "signatures",
    "experiments",
    "programs",
    "compounds",
    "hts_campaigns",
    "hts_plates",  # optional
    "variant_sets",
    "single_cell_datasets",
    "crispr_screens",
    "multi_omics_experiments",
    "ml_models",
    "docking_runs",
    "batch_corrections",
]


def _enable_rls(table: str) -> None:
    op.execute(f'ALTER TABLE IF EXISTS "{table}" ENABLE ROW LEVEL SECURITY;')
    op.execute(f'DROP POLICY IF EXISTS tenant_isolation ON "{table}";')
    op.execute(
        f"""
DO $$
BEGIN
  IF EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = '{table}') THEN
    CREATE POLICY tenant_isolation ON "{table}"
      USING (company_id = current_setting('app.current_company_id', true)::uuid);
  END IF;
END $$;
"""
    )


def upgrade() -> None:
    for t in TENANT_TABLES:
        _enable_rls(t)


def downgrade() -> None:
    for t in TENANT_TABLES:
        op.execute(f'DROP POLICY IF EXISTS tenant_isolation ON "{t}";')
        op.execute(f'ALTER TABLE IF EXISTS "{t}" DISABLE ROW LEVEL SECURITY;')


