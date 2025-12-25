"""Add pinned dashboards table for Voila notebook dashboards.

Revision ID: p1n2d3a4s5h6
Revises: e3b4c5d6e7f8
Create Date: 2025-12-25
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB, UUID


revision = "p1n2d3a4s5h6"
down_revision = "e3b4c5d6e7f8"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "pinned_dashboards",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("notebook_path", sa.String(length=500), nullable=False),
        sa.Column("program_id", UUID(as_uuid=True), sa.ForeignKey("programs.id", ondelete="CASCADE"), nullable=False),
        sa.Column("pinned_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id", ondelete="SET NULL"), nullable=True),
        sa.Column("pinned_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("display_name", sa.String(length=500), nullable=True),
        sa.Column("config", JSONB(), nullable=True),
    )
    op.create_index("ix_pinned_dashboards_program_id", "pinned_dashboards", ["program_id"])
    op.create_index("ix_pinned_dashboards_pinned_by_id", "pinned_dashboards", ["pinned_by_id"])
    op.create_index("ix_pinned_dashboards_notebook_path", "pinned_dashboards", ["notebook_path"])


def downgrade() -> None:
    op.drop_index("ix_pinned_dashboards_notebook_path", table_name="pinned_dashboards")
    op.drop_index("ix_pinned_dashboards_pinned_by_id", table_name="pinned_dashboards")
    op.drop_index("ix_pinned_dashboards_program_id", table_name="pinned_dashboards")
    op.drop_table("pinned_dashboards")


