"""Add company_id to tenant-scoped tables.

Revision ID: f6g7h8i9j0k1
Revises: e1f2g3h4i5j6
Create Date: 2025-12-23
"""

from __future__ import annotations

from alembic import op


revision = "f6g7h8i9j0k1"
down_revision = "e1f2g3h4i5j6"
branch_labels = None
depends_on = None


TENANT_TABLES = [
    # Core
    "datasets",
    "signatures",
    "experiments",
    "programs",
    # Chemistry / HTS
    "compounds",
    "hts_campaigns",
    "hts_plates",  # optional (may not exist in this repo yet)
    # Newer modules
    "variant_sets",
    "single_cell_datasets",
    "crispr_screens",
    "multi_omics_experiments",
    "ml_models",
    "docking_runs",
    "batch_corrections",
]


def _add_company_id(table: str) -> None:
    # Add column
    op.execute(f'ALTER TABLE IF EXISTS "{table}" ADD COLUMN IF NOT EXISTS company_id UUID NULL;')
    # Index
    op.execute(f'CREATE INDEX IF NOT EXISTS ix_{table}_company_id ON "{table}" (company_id);')
    # FK constraint (guarded)
    op.execute(
        f"""
DO $$
BEGIN
  IF EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = '{table}') THEN
    IF NOT EXISTS (SELECT 1 FROM pg_constraint WHERE conname = 'fk_{table}_company_id') THEN
      ALTER TABLE "{table}"
        ADD CONSTRAINT fk_{table}_company_id
        FOREIGN KEY (company_id) REFERENCES companies(id);
    END IF;
  END IF;
END $$;
"""
    )


def upgrade() -> None:
    for t in TENANT_TABLES:
        _add_company_id(t)


def downgrade() -> None:
    for t in TENANT_TABLES:
        op.execute(f'DROP INDEX IF EXISTS ix_{t}_company_id;')
        op.execute(
            f"""
DO $$
BEGIN
  IF EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = '{t}') THEN
    IF EXISTS (SELECT 1 FROM pg_constraint WHERE conname = 'fk_{t}_company_id') THEN
      ALTER TABLE "{t}" DROP CONSTRAINT fk_{t}_company_id;
    END IF;
  END IF;
END $$;
"""
        )
        op.execute(f'ALTER TABLE IF EXISTS "{t}" DROP COLUMN IF EXISTS company_id;')


