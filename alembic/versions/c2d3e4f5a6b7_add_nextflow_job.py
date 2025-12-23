"""Add nextflow_jobs table.

Revision ID: c2d3e4f5a6b7
Revises: b7c8d9e0f1g2
Create Date: 2025-12-22
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID


revision = "c2d3e4f5a6b7"
down_revision = "b7c8d9e0f1g2"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "nextflow_jobs",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("pipeline_name", sa.String(length=100), nullable=False),
        sa.Column("pipeline_version", sa.String(length=20), nullable=True),
        sa.Column("status", sa.String(length=20), server_default="pending", nullable=True),
        sa.Column("sample_sheet_path", sa.String(length=500), nullable=True),
        sa.Column("genome", sa.String(length=50), nullable=True),
        sa.Column("output_dir", sa.String(length=500), nullable=True),
        sa.Column("work_dir", sa.String(length=500), nullable=True),
        sa.Column("nextflow_log", sa.String(length=500), nullable=True),
        sa.Column("multiqc_report", sa.String(length=500), nullable=True),
        sa.Column("progress_percent", sa.Integer(), server_default="0", nullable=True),
        sa.Column("error_message", sa.Text(), nullable=True),
        sa.Column("started_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("completed_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_by", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("params", sa.JSON(), nullable=True),
    )
    op.create_index("ix_nextflow_jobs_created_by", "nextflow_jobs", ["created_by"])
    op.create_index("ix_nextflow_jobs_status", "nextflow_jobs", ["status"])
    op.create_index("ix_nextflow_jobs_created_at", "nextflow_jobs", ["created_at"])


def downgrade() -> None:
    op.drop_index("ix_nextflow_jobs_created_at", table_name="nextflow_jobs")
    op.drop_index("ix_nextflow_jobs_status", table_name="nextflow_jobs")
    op.drop_index("ix_nextflow_jobs_created_by", table_name="nextflow_jobs")
    op.drop_table("nextflow_jobs")


