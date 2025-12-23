"""Add DVC fields to datasets table.

Revision ID: d1e2f3a4b5c6
Revises: c2d3e4f5a6b7
Create Date: 2025-12-23
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "d1e2f3a4b5c6"
down_revision = "c2d3e4f5a6b7"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("datasets", sa.Column("dvc_version", sa.String(length=32), nullable=True))
    op.add_column("datasets", sa.Column("dvc_metadata", sa.JSON(), nullable=True))
    op.add_column(
        "datasets",
        sa.Column("dvc_pushed", sa.Boolean(), server_default=sa.text("false"), nullable=False),
    )


def downgrade() -> None:
    op.drop_column("datasets", "dvc_pushed")
    op.drop_column("datasets", "dvc_metadata")
    op.drop_column("datasets", "dvc_version")


