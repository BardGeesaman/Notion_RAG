"""Add notebook review workflow with signed review cards.

Revision ID: r3v1e2w3s4
Revises: d1g2e3s4t5
Create Date: 2025-12-25
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "r3v1e2w3s4"
down_revision = "d1g2e3s4t5"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "notebook_reviews",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("notebook_path", sa.String(length=500), nullable=False, index=True),
        sa.Column("version_hash", sa.String(length=64), nullable=False, index=True),
        sa.Column("reviewer_id", UUID(as_uuid=True), sa.ForeignKey("users.id", ondelete="SET NULL"), nullable=True),
        sa.Column("status", sa.String(length=32), nullable=False, server_default="pending", index=True),
        sa.Column("comments", sa.Text(), nullable=True),
        sa.Column("signature", sa.String(length=128), nullable=True),
        sa.Column("reviewed_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_notebook_reviews_notebook_status", "notebook_reviews", ["notebook_path", "status"])


def downgrade() -> None:
    op.drop_index("ix_notebook_reviews_notebook_status", table_name="notebook_reviews")
    op.drop_table("notebook_reviews")



