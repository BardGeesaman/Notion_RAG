"""Add batch_corrections table.

Revision ID: f1a2b3c4d5e6
Revises: e7f8g9h0i1j2
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "f1a2b3c4d5e6"
down_revision = "e7f8g9h0i1j2"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "batch_corrections",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("method", sa.String(length=50), nullable=False),
        sa.Column("batch_map", sa.JSON(), nullable=True),
        sa.Column("corrected_dataset_id", UUID(as_uuid=True), sa.ForeignKey("datasets.id"), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )


def downgrade() -> None:
    op.drop_table("batch_corrections")


