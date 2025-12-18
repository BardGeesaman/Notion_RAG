"""merge_all_heads

Revision ID: 59dd5c3fc6a8
Revises: alerts_2024_001, b18cbc8ad60c
Create Date: 2025-12-18 09:36:39.004351

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '59dd5c3fc6a8'
down_revision = ('alerts_2024_001', 'b18cbc8ad60c')
branch_labels = None
depends_on = None


def upgrade() -> None:
    pass


def downgrade() -> None:
    pass

