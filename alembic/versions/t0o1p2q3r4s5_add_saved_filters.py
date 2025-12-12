"""Add saved_filters table

Revision ID: t0o1p2q3r4s5
Revises: s9n0o1p2q3r4
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON

revision = 't0o1p2q3r4s5'
down_revision = 's9n0o1p2q3r4'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'saved_filters',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('name', sa.String(255), nullable=False),
        sa.Column('entity_type', sa.String(50), nullable=False),
        sa.Column('filters', JSON, nullable=False),
        sa.Column('user_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=False),
        sa.Column('created_at', sa.DateTime(), default=sa.func.now()),
    )
    op.create_index('ix_saved_filters_user_id', 'saved_filters', ['user_id'])
    op.create_index('ix_saved_filters_entity_type', 'saved_filters', ['entity_type'])


def downgrade():
    op.drop_index('ix_saved_filters_entity_type')
    op.drop_index('ix_saved_filters_user_id')
    op.drop_table('saved_filters')
