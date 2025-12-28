"""Add created_by_id to core entities

Revision ID: d4e5f6a7b8c9
Revises: c3d4e5f6a7b8
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'd4e5f6a7b8c9'
down_revision = 'c3d4e5f6a7b8'
branch_labels = None
depends_on = None


def upgrade():
    op.add_column('programs', sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True))
    op.add_column('experiments', sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True))
    op.add_column('datasets', sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True))
    op.add_column('signatures', sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True))


def downgrade():
    op.drop_column('signatures', 'created_by_id')
    op.drop_column('datasets', 'created_by_id')
    op.drop_column('experiments', 'created_by_id')
    op.drop_column('programs', 'created_by_id')
