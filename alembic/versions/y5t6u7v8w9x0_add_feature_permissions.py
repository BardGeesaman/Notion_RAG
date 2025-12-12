"""Add feature_permissions table

Revision ID: y5t6u7v8w9x0
Revises: x4s5t6u7v8w9
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'y5t6u7v8w9x0'
down_revision = 'x4s5t6u7v8w9'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'feature_permissions',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('role', sa.String(50), nullable=False),
        sa.Column('page_name', sa.String(100), nullable=False),
        sa.Column('is_visible', sa.Boolean(), default=True),
        sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
        sa.UniqueConstraint('role', 'page_name', name='uq_feature_permissions_role_page'),
    )
    op.create_index('ix_feature_permissions_role', 'feature_permissions', ['role'])


def downgrade():
    op.drop_index('ix_feature_permissions_role', 'feature_permissions')
    op.drop_table('feature_permissions')

