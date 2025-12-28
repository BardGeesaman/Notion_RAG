"""Add notes table for user notes on entities.

Revision ID: s9n0o1p2q3r4
Revises: r8m9n0o1p2q3
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = "s9n0o1p2q3r4"
down_revision = "r8m9n0o1p2q3"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "notes",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("entity_type", sa.String(50), nullable=False),
        sa.Column("entity_id", UUID(as_uuid=True), nullable=False),
        sa.Column("content", sa.Text(), nullable=False),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index("ix_notes_entity", "notes", ["entity_type", "entity_id"])


def downgrade() -> None:
    op.drop_index("ix_notes_entity")
    op.drop_table("notes")
