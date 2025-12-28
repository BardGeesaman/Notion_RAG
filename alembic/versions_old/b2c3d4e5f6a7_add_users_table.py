"""Add users table for authentication

Revision ID: b2c3d4e5f6a7
Revises: f1b8a8e4a9e2
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID

revision = 'b2c3d4e5f6a7'
down_revision = 'f1b8a8e4a9e2'  # Set to actual current head
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'users',
        sa.Column('id', UUID(as_uuid=True), primary_key=True),
        sa.Column('username', sa.String(100), unique=True, nullable=False, index=True),
        sa.Column('email', sa.String(255), unique=True, nullable=False),
        sa.Column('password_hash', sa.String(255), nullable=False),
        sa.Column('role', sa.String(50), default='researcher'),
        sa.Column('is_active', sa.Boolean(), default=True),
        sa.Column('created_at', sa.DateTime(), default=sa.func.now()),
        sa.Column('last_login', sa.DateTime(), nullable=True),
    )


def downgrade():
    op.drop_table('users')
