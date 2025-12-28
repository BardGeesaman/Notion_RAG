"""Add sample inventory tables

Revision ID: g7b8c9d0e1f2
Revises: f6a7b8c9d0e1
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'g7b8c9d0e1f2'
down_revision = 'f6a7b8c9d0e1'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'storage_locations',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('name', sa.String(255), nullable=False),
        sa.Column('location_type', sa.String(50), nullable=True),
        sa.Column('parent_id', UUID(as_uuid=True), sa.ForeignKey('storage_locations.id'), nullable=True),
        sa.Column('temperature', sa.String(50), nullable=True),
        sa.Column('capacity', sa.Integer(), nullable=True),
        sa.Column('description', sa.Text(), nullable=True),
        sa.Column('created_at', sa.DateTime(), default=sa.func.now()),
    )

    op.create_table(
        'samples',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('name', sa.String(255), nullable=False),
        sa.Column('sample_type', sa.String(100), nullable=True),
        sa.Column('barcode', sa.String(255), unique=True, nullable=True),
        sa.Column('storage_location_id', UUID(as_uuid=True), sa.ForeignKey('storage_locations.id'), nullable=True),
        sa.Column('position', sa.String(100), nullable=True),
        sa.Column('parent_sample_id', UUID(as_uuid=True), sa.ForeignKey('samples.id'), nullable=True),
        sa.Column('experiment_id', UUID(as_uuid=True), sa.ForeignKey('experiments.id'), nullable=True),
        sa.Column('quantity', sa.Float(), nullable=True),
        sa.Column('unit', sa.String(50), nullable=True),
        sa.Column('status', sa.String(50), default='available'),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(), default=sa.func.now()),
        sa.Column('notes', sa.Text(), nullable=True),
    )
    op.create_index('ix_samples_barcode', 'samples', ['barcode'], unique=True)
    op.create_index('ix_samples_storage_location_id', 'samples', ['storage_location_id'])

    op.create_table(
        'sample_transfers',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('sample_id', UUID(as_uuid=True), sa.ForeignKey('samples.id'), nullable=False),
        sa.Column('from_location_id', UUID(as_uuid=True), sa.ForeignKey('storage_locations.id'), nullable=True),
        sa.Column('to_location_id', UUID(as_uuid=True), sa.ForeignKey('storage_locations.id'), nullable=True),
        sa.Column('transferred_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('transferred_at', sa.DateTime(), default=sa.func.now()),
        sa.Column('notes', sa.Text(), nullable=True),
    )


def downgrade():
    op.drop_table('sample_transfers')
    op.drop_index('ix_samples_storage_location_id')
    op.drop_index('ix_samples_barcode')
    op.drop_table('samples')
    op.drop_table('storage_locations')
