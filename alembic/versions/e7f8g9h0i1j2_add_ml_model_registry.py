"""Add ML model registry table (ml_models).

Revision ID: e7f8g9h0i1j2
Revises: d1e2f3a4b5c6
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "e7f8g9h0i1j2"
down_revision = "d1e2f3a4b5c6"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "ml_models",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("version", sa.String(length=50), nullable=False),
        sa.Column("model_type", sa.String(length=100), nullable=False),
        sa.Column("framework", sa.String(length=50), nullable=False),
        sa.Column("artifact_path", sa.String(length=500), nullable=False),
        sa.Column("training_dataset_id", UUID(as_uuid=True), sa.ForeignKey("datasets.id"), nullable=True),
        sa.Column("features", sa.JSON(), nullable=True),
        sa.Column("hyperparameters", sa.JSON(), nullable=True),
        sa.Column("metrics", sa.JSON(), nullable=True),
        sa.Column("status", sa.String(length=50), server_default="active", nullable=True),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.UniqueConstraint("name", "version", name="uq_ml_model_name_version"),
    )


def downgrade() -> None:
    op.drop_table("ml_models")


