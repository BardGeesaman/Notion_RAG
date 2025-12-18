"""add alerts table

Revision ID: alerts_2024_001
Revises: c1a2b3c4d5e6
Create Date: 2025-12-18
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID


revision = "alerts_2024_001"
down_revision = "c1a2b3c4d5e6"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "alerts",
        sa.Column("id", UUID(as_uuid=True), primary_key=True, server_default=sa.text("gen_random_uuid()")),
        sa.Column("subscription_id", UUID(as_uuid=True), sa.ForeignKey("repository_subscriptions.id"), nullable=False),
        sa.Column("dataset_id", UUID(as_uuid=True), sa.ForeignKey("datasets.id"), nullable=False),
        sa.Column("is_read", sa.Boolean(), nullable=False, server_default=sa.false()),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
    )
    op.create_index("ix_alerts_subscription_id", "alerts", ["subscription_id"])
    op.create_index("ix_alerts_dataset_id", "alerts", ["dataset_id"])


def downgrade():
    op.drop_index("ix_alerts_dataset_id", table_name="alerts")
    op.drop_index("ix_alerts_subscription_id", table_name="alerts")
    op.drop_table("alerts")

