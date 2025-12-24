"""Add extraction_jobs and extracted_documents tables, and link rag_chunks to extraction jobs.

Revision ID: e2a3b4c5d6f7
Revises: d8e9f0a1b2c3
Create Date: 2025-12-24
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB, UUID


revision = "e2a3b4c5d6f7"
down_revision = "d8e9f0a1b2c3"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "extraction_jobs",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("batch_id", UUID(as_uuid=True), nullable=True),
        sa.Column("file_count", sa.Integer(), nullable=False),
        sa.Column("completed_count", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("status", sa.String(length=50), nullable=False, server_default="pending"),
        sa.Column("started_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("completed_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_extraction_jobs_batch_id", "extraction_jobs", ["batch_id"])
    op.create_index("ix_extraction_jobs_status", "extraction_jobs", ["status"])
    op.create_index("ix_extraction_jobs_batch_status", "extraction_jobs", ["batch_id", "status"])

    op.create_table(
        "extracted_documents",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column(
            "job_id",
            UUID(as_uuid=True),
            sa.ForeignKey("extraction_jobs.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("file_path", sa.String(length=500), nullable=False),
        sa.Column("original_filename", sa.String(length=500), nullable=False),
        sa.Column("doc_type", sa.String(length=50), nullable=False),
        sa.Column("extracted_entities", JSONB(), nullable=False),
        sa.Column("extraction_config", JSONB(), nullable=True),
        sa.Column("status", sa.String(length=50), nullable=False, server_default="pending"),
        sa.Column("error_log", sa.Text(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_extracted_documents_job_id", "extracted_documents", ["job_id"])
    op.create_index("ix_extracted_documents_doc_type", "extracted_documents", ["doc_type"])
    op.create_index("ix_extracted_documents_status", "extracted_documents", ["status"])
    op.create_index(
        "ix_extracted_documents_job_status",
        "extracted_documents",
        ["job_id", "status"],
    )

    op.add_column(
        "rag_chunks",
        sa.Column("extraction_job_id", UUID(as_uuid=True), sa.ForeignKey("extraction_jobs.id"), nullable=True),
    )
    op.create_index("ix_rag_chunks_extraction_job_id", "rag_chunks", ["extraction_job_id"])


def downgrade() -> None:
    op.drop_index("ix_rag_chunks_extraction_job_id", table_name="rag_chunks")
    op.drop_column("rag_chunks", "extraction_job_id")

    op.drop_index("ix_extracted_documents_job_status", table_name="extracted_documents")
    op.drop_index("ix_extracted_documents_status", table_name="extracted_documents")
    op.drop_index("ix_extracted_documents_doc_type", table_name="extracted_documents")
    op.drop_index("ix_extracted_documents_job_id", table_name="extracted_documents")
    op.drop_table("extracted_documents")

    op.drop_index("ix_extraction_jobs_batch_status", table_name="extraction_jobs")
    op.drop_index("ix_extraction_jobs_status", table_name="extraction_jobs")
    op.drop_index("ix_extraction_jobs_batch_id", table_name="extraction_jobs")
    op.drop_table("extraction_jobs")


