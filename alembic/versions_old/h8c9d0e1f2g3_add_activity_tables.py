"""Add activity tables for biochemical assays and results.

Revision ID: h8c9d0e1f2g3
Revises: g7b8c9d0e1f2
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = "h8c9d0e1f2g3"
down_revision = "g7b8c9d0e1f2"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Biochemical assays table
    op.create_table(
        "biochemical_assays",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(200), nullable=False),
        sa.Column("target", sa.String(200), nullable=True),
        sa.Column("assay_type", sa.String(50), nullable=True),
        sa.Column("unit", sa.String(20), nullable=True),
        sa.Column("description", sa.Text, nullable=True),
        sa.Column("created_at", sa.DateTime, nullable=False),
    )

    # Activity results table
    op.create_table(
        "activity_results",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("compound_id", UUID(as_uuid=True), sa.ForeignKey("compounds.id"), nullable=False),
        sa.Column("assay_id", UUID(as_uuid=True), sa.ForeignKey("biochemical_assays.id"), nullable=False),
        sa.Column("value", sa.Float, nullable=False),
        sa.Column("qualifier", sa.String(5), nullable=True),
        sa.Column("created_at", sa.DateTime, nullable=False),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
    )

    # Indexes
    op.create_index("ix_activity_results_compound_id", "activity_results", ["compound_id"])
    op.create_index("ix_activity_results_assay_id", "activity_results", ["assay_id"])


def downgrade() -> None:
    op.drop_index("ix_activity_results_assay_id")
    op.drop_index("ix_activity_results_compound_id")
    op.drop_table("activity_results")
    op.drop_table("biochemical_assays")
