"""Add genomics pipeline tables (genomics_indices, pipeline_jobs).

Revision ID: b7c8d9e0f1g2
Revises: a1b2c3
Create Date: 2025-12-22
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID


revision = "b7c8d9e0f1g2"
down_revision = "a1b2c3"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "genomics_indices",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("organism", sa.String(length=100), nullable=False),
        sa.Column("tool", sa.String(length=20), nullable=False),
        sa.Column("version", sa.String(length=50), nullable=False),
        sa.Column("file_path", sa.String(length=500), nullable=False),
        sa.Column("file_size_bytes", sa.BigInteger(), nullable=True),
        sa.Column("uploaded_by", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("uploaded_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=True),
        sa.Column("metadata", sa.JSON(), nullable=True),
        sa.UniqueConstraint("organism", "tool", "version", name="uix_index_org_tool_ver"),
    )
    op.create_index("ix_genomics_indices_tool", "genomics_indices", ["tool"])
    op.create_index("ix_genomics_indices_organism", "genomics_indices", ["organism"])

    op.create_table(
        "pipeline_jobs",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("status", sa.String(length=20), server_default="pending", nullable=True),
        sa.Column("tool", sa.String(length=20), nullable=False),
        sa.Column("input_fastq_path", sa.String(length=500), nullable=True),
        sa.Column("index_id", UUID(as_uuid=True), sa.ForeignKey("genomics_indices.id"), nullable=True),
        sa.Column("output_dir", sa.String(length=500), nullable=True),
        sa.Column("result_file", sa.String(length=500), nullable=True),
        sa.Column("progress_percent", sa.Integer(), server_default="0", nullable=True),
        sa.Column("error_message", sa.Text(), nullable=True),
        sa.Column("started_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("completed_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_by", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_pipeline_jobs_created_by", "pipeline_jobs", ["created_by"])
    op.create_index("ix_pipeline_jobs_status", "pipeline_jobs", ["status"])
    op.create_index("ix_pipeline_jobs_created_at", "pipeline_jobs", ["created_at"])


def downgrade() -> None:
    op.drop_index("ix_pipeline_jobs_created_at", table_name="pipeline_jobs")
    op.drop_index("ix_pipeline_jobs_status", table_name="pipeline_jobs")
    op.drop_index("ix_pipeline_jobs_created_by", table_name="pipeline_jobs")
    op.drop_table("pipeline_jobs")

    op.drop_index("ix_genomics_indices_organism", table_name="genomics_indices")
    op.drop_index("ix_genomics_indices_tool", table_name="genomics_indices")
    op.drop_table("genomics_indices")


