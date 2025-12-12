"""Add genetic_variants table

Revision ID: g3b4c5d6e7f8
Revises: f2a3b4c5d6e7
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'g3b4c5d6e7f8'
down_revision = 'f2a3b4c5d6e7'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'genetic_variants',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('gene', sa.String(), nullable=False),
        sa.Column('variant', sa.String(), nullable=False),
        sa.Column('zygosity', sa.String(), nullable=True),
        sa.Column('organism', sa.String(), nullable=False),
        sa.Column('experiment_id', UUID(as_uuid=True), sa.ForeignKey('experiments.id'), nullable=True),
        sa.Column('notes', sa.Text(), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.func.now()),
    )
    
    # Create index on gene
    op.create_index('ix_genetic_variants_gene', 'genetic_variants', ['gene'])


def downgrade():
    op.drop_index('ix_genetic_variants_gene', 'genetic_variants')
    op.drop_table('genetic_variants')

