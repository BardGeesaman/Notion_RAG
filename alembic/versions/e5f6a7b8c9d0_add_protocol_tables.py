"""Add protocol management tables

Revision ID: e5f6a7b8c9d0
Revises: d4e5f6a7b8c9
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON

revision = 'e5f6a7b8c9d0'
down_revision = 'd4e5f6a7b8c9'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'protocols',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('name', sa.String(255), nullable=False),
        sa.Column('version', sa.Integer(), default=1),
        sa.Column('category', sa.String(100), nullable=True),
        sa.Column('description', sa.Text(), nullable=True),
        sa.Column('steps', JSON, nullable=True),
        sa.Column('materials', JSON, nullable=True),
        sa.Column('parameters', JSON, nullable=True),
        sa.Column('is_template', sa.Boolean(), default=True),
        sa.Column('is_active', sa.Boolean(), default=True),
        sa.Column('parent_id', UUID(as_uuid=True), sa.ForeignKey('protocols.id'), nullable=True),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(), default=sa.func.now()),
        sa.Column('updated_at', sa.DateTime(), default=sa.func.now()),
    )
    
    op.create_table(
        'experiment_protocols',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('experiment_id', UUID(as_uuid=True), sa.ForeignKey('experiments.id'), nullable=False),
        sa.Column('protocol_id', UUID(as_uuid=True), sa.ForeignKey('protocols.id'), nullable=False),
        sa.Column('executed_at', sa.DateTime(), nullable=True),
        sa.Column('notes', sa.Text(), nullable=True),
        sa.Column('deviations', JSON, nullable=True),
    )


def downgrade():
    op.drop_table('experiment_protocols')
    op.drop_table('protocols')
