"""add repository_subscription table

Revision ID: c1a2b3c4d5e6
Revises: 0225cbc008f1
Create Date: 2025-12-18
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON


revision = "c1a2b3c4d5e6"
down_revision = "0225cbc008f1"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "repository_subscriptions",
        sa.Column("id", UUID(as_uuid=True), primary_key=True, server_default=sa.text("gen_random_uuid()")),
        sa.Column("user_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("repository_source", sa.String(length=50), nullable=False),
        sa.Column("query_params", JSON, nullable=True),
        sa.Column("notify_email", sa.Boolean(), nullable=False, server_default=sa.false()),
        sa.Column("notify_in_app", sa.Boolean(), nullable=False, server_default=sa.true()),
        sa.Column("is_active", sa.Boolean(), nullable=False, server_default=sa.true()),
        sa.Column("last_checked", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column(
            "updated_at",
            sa.DateTime(timezone=True),
            server_default=sa.func.now(),
            onupdate=sa.func.now(),
            nullable=False,
        ),
    )
    op.create_index("ix_repository_subscriptions_user_id", "repository_subscriptions", ["user_id"])
    op.create_index("ix_repository_subscriptions_source", "repository_subscriptions", ["repository_source"])


def downgrade():
    op.drop_index("ix_repository_subscriptions_source", table_name="repository_subscriptions")
    op.drop_index("ix_repository_subscriptions_user_id", table_name="repository_subscriptions")
    op.drop_table("repository_subscriptions")

