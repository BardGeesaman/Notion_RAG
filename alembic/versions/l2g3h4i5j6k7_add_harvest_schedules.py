"""Add harvest_schedules table for automated repository harvesting.

Revision ID: l2g3h4i5j6k7
Revises: k1f2g3h4i5j6
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = "l2g3h4i5j6k7"
down_revision = "k1f2g3h4i5j6"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "harvest_schedules",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(255), nullable=False),
        sa.Column("repository", sa.String(50), nullable=False),
        sa.Column("query", sa.String(500), nullable=False),
        sa.Column("interval_hours", sa.Integer(), default=24),
        sa.Column("is_active", sa.Boolean(), default=True),
        sa.Column("last_run", sa.DateTime(timezone=True), nullable=True),
        sa.Column("next_run", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index("ix_harvest_schedules_is_active", "harvest_schedules", ["is_active"])


def downgrade() -> None:
    op.drop_index("ix_harvest_schedules_is_active")
    op.drop_table("harvest_schedules")
