"""add_auth_lockout_table

Revision ID: aea7b4fd12f7
Revises: 6327e1a98d8a
Create Date: 2026-01-02 09:00:45.549887

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'aea7b4fd12f7'
down_revision = '6327e1a98d8a'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        'auth_lockout_attempts',
        sa.Column('id', sa.Integer(), primary_key=True),
        sa.Column('identifier', sa.String(255), nullable=False),
        sa.Column('attempt_time', sa.DateTime(timezone=True), nullable=False),
        sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
    )
    op.create_index('ix_lockout_identifier', 'auth_lockout_attempts', ['identifier'])
    op.create_index('ix_lockout_identifier_time', 'auth_lockout_attempts', ['identifier', 'attempt_time'])


def downgrade() -> None:
    op.drop_index('ix_lockout_identifier_time')
    op.drop_index('ix_lockout_identifier')
    op.drop_table('auth_lockout_attempts')

