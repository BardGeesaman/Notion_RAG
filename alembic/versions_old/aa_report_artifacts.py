"""add report_artifacts table

Revision ID: aa_report_artifacts
Revises: z6u7v8w9x0y1
Create Date: 2025-12-17
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID


revision = "aa_report_artifacts"
down_revision = "z6u7v8w9x0y1"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "report_artifacts",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("entity_type", sa.String(length=50), nullable=False),
        sa.Column("entity_id", UUID(as_uuid=True), nullable=False),
        sa.Column("format", sa.String(length=10), nullable=False),
        sa.Column("file_path", sa.String(), nullable=False),
        sa.Column("params_hash", sa.String(length=128), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
    )
    op.create_index("ix_report_artifacts_entity_type", "report_artifacts", ["entity_type"])
    op.create_index("ix_report_artifacts_entity_id", "report_artifacts", ["entity_id"])
    op.create_index("ix_report_artifacts_params_hash", "report_artifacts", ["params_hash"])


def downgrade():
    op.drop_index("ix_report_artifacts_params_hash", table_name="report_artifacts")
    op.drop_index("ix_report_artifacts_entity_id", table_name="report_artifacts")
    op.drop_index("ix_report_artifacts_entity_type", table_name="report_artifacts")
    op.drop_table("report_artifacts")

