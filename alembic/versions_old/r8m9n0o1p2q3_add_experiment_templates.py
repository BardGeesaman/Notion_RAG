"""Add experiment_templates table for reusable experiment templates.

Revision ID: r8m9n0o1p2q3
Revises: q7l8m9n0o1p2
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, ARRAY, JSON

revision = "r8m9n0o1p2q3"
down_revision = "q7l8m9n0o1p2"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "experiment_templates",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(255), nullable=False),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("design_type", sa.String(50), nullable=True),
        sa.Column("organism", ARRAY(sa.String()), nullable=True),
        sa.Column("sample_groups", JSON, nullable=True),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now()),
    )


def downgrade() -> None:
    op.drop_table("experiment_templates")
