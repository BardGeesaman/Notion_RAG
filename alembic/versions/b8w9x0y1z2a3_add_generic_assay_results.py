"""Add generic_assay_results table

Revision ID: b8w9x0y1z2a3
Revises: a7v8w9x0y1z2
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON

revision = 'b8w9x0y1z2a3'
down_revision = 'a7v8w9x0y1z2'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'generic_assay_results',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('assay_name', sa.String(200), nullable=False),
        sa.Column('assay_type', sa.String(100), nullable=False),
        sa.Column('experiment_id', UUID(as_uuid=True), sa.ForeignKey('experiments.id'), nullable=True),
        sa.Column('compound_id', UUID(as_uuid=True), sa.ForeignKey('compounds.id'), nullable=True),
        sa.Column('result_data', JSON(), nullable=False),
        sa.Column('metadata', JSON(), nullable=True),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.func.now()),
    )
    
    # Create indexes
    op.create_index('ix_generic_assay_results_assay_name', 'generic_assay_results', ['assay_name'])
    op.create_index('ix_generic_assay_results_assay_type', 'generic_assay_results', ['assay_type'])
    op.create_index('ix_generic_assay_results_experiment_id', 'generic_assay_results', ['experiment_id'])
    op.create_index('ix_generic_assay_results_compound_id', 'generic_assay_results', ['compound_id'])


def downgrade():
    op.drop_index('ix_generic_assay_results_compound_id', 'generic_assay_results')
    op.drop_index('ix_generic_assay_results_experiment_id', 'generic_assay_results')
    op.drop_index('ix_generic_assay_results_assay_type', 'generic_assay_results')
    op.drop_index('ix_generic_assay_results_assay_name', 'generic_assay_results')
    op.drop_table('generic_assay_results')

