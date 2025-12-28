"""Add version columns for concurrent editing safety

Revision ID: z6u7v8w9x0y1
Revises: y5t6u7v8w9x0
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa

revision = 'z6u7v8w9x0y1'
down_revision = 'y5t6u7v8w9x0'
branch_labels = None
depends_on = None


def upgrade():
    # Add version column to experiments
    op.add_column('experiments', sa.Column('version', sa.Integer(), nullable=False, server_default='1'))

    # Add version column to datasets
    op.add_column('datasets', sa.Column('version', sa.Integer(), nullable=False, server_default='1'))

    # Add version column to compounds
    op.add_column('compounds', sa.Column('version', sa.Integer(), nullable=False, server_default='1'))

    # Update version column in protocols (ensure nullable=False)
    op.alter_column('protocols', 'version',
                    existing_type=sa.Integer(),
                    nullable=False,
                    server_default='1')


def downgrade():
    # Remove version column from experiments
    op.drop_column('experiments', 'version')

    # Remove version column from datasets
    op.drop_column('datasets', 'version')

    # Remove version column from compounds
    op.drop_column('compounds', 'version')

    # Revert version column in protocols (make nullable again)
    op.alter_column('protocols', 'version',
                    existing_type=sa.Integer(),
                    nullable=True,
                    server_default=None)

