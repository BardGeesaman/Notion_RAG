"""Add docking_runs.exhaustiveness column.

Revision ID: b2c3d4e6f7a8
Revises: a1b2c3d4e6f7
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op


revision = "b2c3d4e6f7a8"
down_revision = "a1b2c3d4e6f7"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "docking_runs",
        sa.Column("exhaustiveness", sa.Integer(), nullable=False, server_default="8"),
    )


def downgrade() -> None:
    op.drop_column("docking_runs", "exhaustiveness")


