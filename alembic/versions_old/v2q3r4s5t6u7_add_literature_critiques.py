"""Add literature_critiques table

Revision ID: v2q3r4s5t6u7
Revises: u1p2q3r4s5t6
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON

revision = 'v2q3r4s5t6u7'
down_revision = 'u1p2q3r4s5t6'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'literature_critiques',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('source_type', sa.String(50), nullable=False),
        sa.Column('source_id', UUID(as_uuid=True), nullable=True),
        sa.Column('source_text', sa.Text(), nullable=False),
        sa.Column('critique_type', sa.String(50), nullable=False),
        sa.Column('content', JSON, nullable=False),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index('ix_literature_critiques_critique_type', 'literature_critiques', ['critique_type'])
    op.create_index('ix_literature_critiques_source_type', 'literature_critiques', ['source_type'])


def downgrade():
    op.drop_index('ix_literature_critiques_source_type')
    op.drop_index('ix_literature_critiques_critique_type')
    op.drop_table('literature_critiques')
