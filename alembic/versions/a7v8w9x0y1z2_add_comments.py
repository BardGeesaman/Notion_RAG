"""Add comments table for contextual commenting

Revision ID: a7v8w9x0y1z2
Revises: z6u7v8w9x0y1
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'a7v8w9x0y1z2'
down_revision = 'z6u7v8w9x0y1'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'comments',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('entity_type', sa.String(50), nullable=False),
        sa.Column('entity_id', UUID(as_uuid=True), nullable=False),
        sa.Column('parent_id', UUID(as_uuid=True), sa.ForeignKey('comments.id'), nullable=True),
        sa.Column('content', sa.Text(), nullable=False),
        sa.Column('created_by_id', UUID(as_uuid=True), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.func.now()),
        sa.Column('updated_at', sa.DateTime(), nullable=True),
    )

    # Create indexes
    op.create_index('ix_comments_entity', 'comments', ['entity_type', 'entity_id'])
    op.create_index('ix_comments_parent_id', 'comments', ['parent_id'])


def downgrade():
    op.drop_index('ix_comments_parent_id', 'comments')
    op.drop_index('ix_comments_entity', 'comments')
    op.drop_table('comments')

