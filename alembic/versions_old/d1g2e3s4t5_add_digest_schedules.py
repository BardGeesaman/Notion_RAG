"""Add digest schedules for weekly executive notebook digests.

Revision ID: d1g2e3s4t5
Revises: p1n2d3a4s5h6
Create Date: 2025-12-25
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB, UUID


revision = "d1g2e3s4t5"
down_revision = "p1n2d3a4s5h6"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "digest_schedules",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("program_id", UUID(as_uuid=True), sa.ForeignKey("programs.id", ondelete="CASCADE"), nullable=False),
        sa.Column("notebook_path", sa.String(length=500), nullable=False),
        sa.Column("schedule_cron", sa.String(length=100), nullable=False),
        sa.Column("recipients", JSONB(), nullable=False),
        sa.Column("last_run_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("last_status", sa.String(length=20), nullable=True),
        sa.Column("enabled", sa.Boolean(), nullable=False, server_default=sa.text("true")),
    )
    op.create_index("ix_digest_schedules_program_id", "digest_schedules", ["program_id"])
    op.create_index("ix_digest_schedules_enabled", "digest_schedules", ["enabled"])
    op.create_index("ix_digest_schedules_program_enabled", "digest_schedules", ["program_id", "enabled"])


def downgrade() -> None:
    op.drop_index("ix_digest_schedules_program_enabled", table_name="digest_schedules")
    op.drop_index("ix_digest_schedules_enabled", table_name="digest_schedules")
    op.drop_index("ix_digest_schedules_program_id", table_name="digest_schedules")
    op.drop_table("digest_schedules")



