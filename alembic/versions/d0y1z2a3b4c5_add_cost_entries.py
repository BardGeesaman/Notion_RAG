"""Add cost_entries table

Revision ID: d0y1z2a3b4c5
Revises: c9x0y1z2a3b4
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'd0y1z2a3b4c5'
down_revision = 'c9x0y1z2a3b4'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'cost_entries',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('project_id', UUID(as_uuid=True), sa.ForeignKey('projects.id'), nullable=True),
        sa.Column('experiment_id', UUID(as_uuid=True), sa.ForeignKey('experiments.id'), nullable=True),
        sa.Column('category', sa.String(100), nullable=False),
        sa.Column('description', sa.Text(), nullable=False),
        sa.Column('amount', sa.Float(), nullable=False),
        sa.Column('currency', sa.String(10), nullable=False, server_default='USD'),
        sa.Column('entry_date', sa.DateTime(), nullable=False),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.func.now()),
    )

    # Create indexes
    op.create_index('ix_cost_entries_project_id', 'cost_entries', ['project_id'])
    op.create_index('ix_cost_entries_experiment_id', 'cost_entries', ['experiment_id'])
    op.create_index('ix_cost_entries_category', 'cost_entries', ['category'])
    op.create_index('ix_cost_entries_entry_date', 'cost_entries', ['entry_date'])


def downgrade():
    op.drop_index('ix_cost_entries_entry_date', 'cost_entries')
    op.drop_index('ix_cost_entries_category', 'cost_entries')
    op.drop_index('ix_cost_entries_experiment_id', 'cost_entries')
    op.drop_index('ix_cost_entries_project_id', 'cost_entries')
    op.drop_table('cost_entries')

