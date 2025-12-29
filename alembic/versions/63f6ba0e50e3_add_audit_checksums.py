"""add audit checksums

Revision ID: 63f6ba0e50e3
Revises: 8f4f671e4f2f
Create Date: 2025-12-29 11:01:46.408197

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '63f6ba0e50e3'
down_revision = '8f4f671e4f2f'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column('audit_logs', sa.Column('old_checksum', sa.String(length=64), nullable=True))
    op.add_column('audit_logs', sa.Column('new_checksum', sa.String(length=64), nullable=True))
    op.add_column('audit_logs', sa.Column('changes', sa.JSON(), nullable=True))


def downgrade() -> None:
    op.drop_column('audit_logs', 'changes')
    op.drop_column('audit_logs', 'new_checksum')
    op.drop_column('audit_logs', 'old_checksum')

