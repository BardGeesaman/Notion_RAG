"""Add candidate selection tables

Revision ID: w3r4s5t6u7v8
Revises: v2q3r4s5t6u7
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON

revision = 'w3r4s5t6u7v8'
down_revision = 'v2q3r4s5t6u7'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'target_product_profiles',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('name', sa.String(255), nullable=False),
        sa.Column('description', sa.Text(), nullable=True),
        sa.Column('criteria', JSON, nullable=False),
        sa.Column('is_active', sa.Boolean(), default=True),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index('ix_target_product_profiles_is_active', 'target_product_profiles', ['is_active'])

    op.create_table(
        'candidate_nominations',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('compound_id', UUID(as_uuid=True), sa.ForeignKey('compounds.id'), nullable=False),
        sa.Column('tpp_id', UUID(as_uuid=True), sa.ForeignKey('target_product_profiles.id'), nullable=False),
        sa.Column('scores', JSON, nullable=True),
        sa.Column('overall_score', sa.Float(), nullable=True),
        sa.Column('status', sa.String(50), default='nominated'),
        sa.Column('notes', sa.Text(), nullable=True),
        sa.Column('nominated_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index('ix_candidate_nominations_compound_id', 'candidate_nominations', ['compound_id'])
    op.create_index('ix_candidate_nominations_tpp_id', 'candidate_nominations', ['tpp_id'])
    op.create_index('ix_candidate_nominations_status', 'candidate_nominations', ['status'])


def downgrade():
    op.drop_index('ix_candidate_nominations_status')
    op.drop_index('ix_candidate_nominations_tpp_id')
    op.drop_index('ix_candidate_nominations_compound_id')
    op.drop_table('candidate_nominations')
    op.drop_index('ix_target_product_profiles_is_active')
    op.drop_table('target_product_profiles')
