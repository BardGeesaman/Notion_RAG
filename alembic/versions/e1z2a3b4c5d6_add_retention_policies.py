"""Add retention_policies table and archive fields

Revision ID: e1z2a3b4c5d6
Revises: d0y1z2a3b4c5
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'e1z2a3b4c5d6'
down_revision = 'd0y1z2a3b4c5'
branch_labels = None
depends_on = None


def upgrade():
    # Create retention_policies table
    op.create_table(
        'retention_policies',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('name', sa.String(200), nullable=False),
        sa.Column('entity_type', sa.String(50), nullable=False),
        sa.Column('retention_days', sa.Integer(), nullable=False),
        sa.Column('action', sa.String(50), nullable=False),
        sa.Column('is_active', sa.Boolean(), nullable=False, server_default='true'),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.func.now()),
    )
    
    # Create index on retention_policies
    op.create_index('ix_retention_policies_entity_type', 'retention_policies', ['entity_type'])
    op.create_index('ix_retention_policies_is_active', 'retention_policies', ['is_active'])
    
    # Add archive fields to experiments
    op.add_column('experiments', sa.Column('is_archived', sa.Boolean(), nullable=False, server_default='false'))
    op.add_column('experiments', sa.Column('archived_at', sa.DateTime(), nullable=True))
    
    # Add archive fields to datasets
    op.add_column('datasets', sa.Column('is_archived', sa.Boolean(), nullable=False, server_default='false'))
    op.add_column('datasets', sa.Column('archived_at', sa.DateTime(), nullable=True))


def downgrade():
    # Remove archive fields from datasets
    op.drop_column('datasets', 'archived_at')
    op.drop_column('datasets', 'is_archived')
    
    # Remove archive fields from experiments
    op.drop_column('experiments', 'archived_at')
    op.drop_column('experiments', 'is_archived')
    
    # Drop retention_policies table
    op.drop_index('ix_retention_policies_is_active', 'retention_policies')
    op.drop_index('ix_retention_policies_entity_type', 'retention_policies')
    op.drop_table('retention_policies')

