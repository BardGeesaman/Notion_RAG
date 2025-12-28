"""Add scheduled_events table

Revision ID: c9x0y1z2a3b4
Revises: b8w9x0y1z2a3
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'c9x0y1z2a3b4'
down_revision = 'b8w9x0y1z2a3'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'scheduled_events',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('title', sa.String(500), nullable=False),
        sa.Column('event_type', sa.String(100), nullable=False),
        sa.Column('resource_name', sa.String(200), nullable=False),
        sa.Column('start_time', sa.DateTime(), nullable=False),
        sa.Column('end_time', sa.DateTime(), nullable=False),
        sa.Column('experiment_id', UUID(as_uuid=True), sa.ForeignKey('experiments.id'), nullable=True),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('notes', sa.Text(), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.func.now()),
    )

    # Create indexes
    op.create_index('ix_scheduled_events_event_type', 'scheduled_events', ['event_type'])
    op.create_index('ix_scheduled_events_resource_name', 'scheduled_events', ['resource_name'])
    op.create_index('ix_scheduled_events_start_time', 'scheduled_events', ['start_time'])
    op.create_index('ix_scheduled_events_end_time', 'scheduled_events', ['end_time'])
    op.create_index('ix_scheduled_events_experiment_id', 'scheduled_events', ['experiment_id'])


def downgrade():
    op.drop_index('ix_scheduled_events_experiment_id', 'scheduled_events')
    op.drop_index('ix_scheduled_events_end_time', 'scheduled_events')
    op.drop_index('ix_scheduled_events_start_time', 'scheduled_events')
    op.drop_index('ix_scheduled_events_resource_name', 'scheduled_events')
    op.drop_index('ix_scheduled_events_event_type', 'scheduled_events')
    op.drop_table('scheduled_events')

