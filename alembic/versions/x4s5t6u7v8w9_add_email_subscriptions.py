"""Add email_subscriptions table

Revision ID: x4s5t6u7v8w9
Revises: w3r4s5t6u7v8
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'x4s5t6u7v8w9'
down_revision = 'w3r4s5t6u7v8'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'email_subscriptions',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('user_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=False),
        sa.Column('subscription_type', sa.String(50), nullable=False),
        sa.Column('frequency', sa.String(50), nullable=False),
        sa.Column('is_active', sa.Boolean(), default=True),
        sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index('ix_email_subscriptions_user_id', 'email_subscriptions', ['user_id'])


def downgrade():
    op.drop_index('ix_email_subscriptions_user_id', 'email_subscriptions')
    op.drop_table('email_subscriptions')

