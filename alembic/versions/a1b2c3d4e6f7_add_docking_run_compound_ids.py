"""Add docking_runs.compound_ids column.

Revision ID: a1b2c3d4e6f7
Revises: f0a1b2c3d4e5
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op


revision = "a1b2c3d4e6f7"
down_revision = "f0a1b2c3d4e5"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("docking_runs", sa.Column("compound_ids", sa.JSON(), nullable=True))


def downgrade() -> None:
    op.drop_column("docking_runs", "compound_ids")


