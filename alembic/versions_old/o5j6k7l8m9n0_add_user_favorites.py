"""Add user_favorites table for user bookmarks.

Revision ID: o5j6k7l8m9n0
Revises: n4i5j6k7l8m9
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = "o5j6k7l8m9n0"
down_revision = "n4i5j6k7l8m9"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "user_favorites",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("user_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=False),
        sa.Column("page_name", sa.String(100), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index("ix_user_favorites_user_id", "user_favorites", ["user_id"])


def downgrade() -> None:
    op.drop_index("ix_user_favorites_user_id")
    op.drop_table("user_favorites")
