"""Add ontology_terms table

Revision ID: f2a3b4c5d6e7
Revises: e1z2a3b4c5d6
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'f2a3b4c5d6e7'
down_revision = 'e1z2a3b4c5d6'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'ontology_terms',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('vocabulary', sa.String(100), nullable=False),
        sa.Column('term', sa.String(500), nullable=False),
        sa.Column('description', sa.Text(), nullable=True),
        sa.Column('parent_id', UUID(as_uuid=True), sa.ForeignKey('ontology_terms.id'), nullable=True),
        sa.Column('is_active', sa.Boolean(), nullable=False, server_default='true'),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.func.now()),
    )
    
    # Create indexes
    op.create_index('ix_ontology_terms_vocabulary', 'ontology_terms', ['vocabulary'])
    op.create_index('ix_ontology_terms_term', 'ontology_terms', ['term'])
    op.create_index('ix_ontology_terms_parent_id', 'ontology_terms', ['parent_id'])
    
    # Create unique constraint
    op.create_unique_constraint('uq_ontology_terms_vocab_term', 'ontology_terms', ['vocabulary', 'term'])


def downgrade():
    op.drop_constraint('uq_ontology_terms_vocab_term', 'ontology_terms', type_='unique')
    op.drop_index('ix_ontology_terms_parent_id', 'ontology_terms')
    op.drop_index('ix_ontology_terms_term', 'ontology_terms')
    op.drop_index('ix_ontology_terms_vocabulary', 'ontology_terms')
    op.drop_table('ontology_terms')

