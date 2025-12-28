"""Add discovery job tables

Revision ID: f6a7b8c9d0e1
Revises: e5f6a7b8c9d0
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'f6a7b8c9d0e1'
down_revision = 'e5f6a7b8c9d0'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'discovery_jobs',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('repository', sa.String(50), nullable=False),
        sa.Column('query', sa.String(500), nullable=False),
        sa.Column('status', sa.String(50), default='pending'),
        sa.Column('studies_found', sa.Integer(), default=0),
        sa.Column('studies_imported', sa.Integer(), default=0),
        sa.Column('error_message', sa.Text(), nullable=True),
        sa.Column('started_at', sa.DateTime(), nullable=True),
        sa.Column('completed_at', sa.DateTime(), nullable=True),
        sa.Column('created_at', sa.DateTime(), default=sa.func.now()),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
    )

    op.create_table(
        'discovered_studies',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('job_id', UUID(as_uuid=True), sa.ForeignKey('discovery_jobs.id'), nullable=False),
        sa.Column('study_id', sa.String(100), nullable=False),
        sa.Column('repository', sa.String(50), nullable=False),
        sa.Column('title', sa.String(500), nullable=True),
        sa.Column('description', sa.Text(), nullable=True),
        sa.Column('organism', sa.String(200), nullable=True),
        sa.Column('omics_type', sa.String(100), nullable=True),
        sa.Column('status', sa.String(50), default='new'),
        sa.Column('imported_experiment_id', UUID(as_uuid=True), sa.ForeignKey('experiments.id'), nullable=True),
        sa.Column('discovered_at', sa.DateTime(), default=sa.func.now()),
    )

    op.create_index('ix_discovered_studies_job_id', 'discovered_studies', ['job_id'])
    op.create_index('ix_discovered_studies_status', 'discovered_studies', ['status'])


def downgrade():
    op.drop_index('ix_discovered_studies_status')
    op.drop_index('ix_discovered_studies_job_id')
    op.drop_table('discovered_studies')
    op.drop_table('discovery_jobs')
