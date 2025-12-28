"""add report_schedules table

Revision ID: ab_report_schedules
Revises: aa_report_artifacts
Create Date: 2025-12-17
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID


revision = "ab_report_schedules"
down_revision = "aa_report_artifacts"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "report_schedules",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("entity_type", sa.String(length=50), nullable=False),
        sa.Column("entity_id", UUID(as_uuid=True), nullable=True),
        sa.Column("format", sa.String(length=10), nullable=False),
        sa.Column("cron_expression", sa.String(length=100), nullable=False),
        sa.Column("enabled", sa.Boolean(), nullable=False, server_default=sa.sql.expression.true()),
        sa.Column("last_run_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_report_schedules_entity_type", "report_schedules", ["entity_type"])
    op.create_index("ix_report_schedules_entity_id", "report_schedules", ["entity_id"])


def downgrade():
    op.drop_index("ix_report_schedules_entity_id", table_name="report_schedules")
    op.drop_index("ix_report_schedules_entity_type", table_name="report_schedules")
    op.drop_table("report_schedules")

