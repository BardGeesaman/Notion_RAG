"""Add annotation_type to notes table.

Revision ID: i0j1k2l3m4n5
Revises: h9i0j1k2l3m4
Create Date: 2025-12-14
"""

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "i0j1k2l3m4n5"
down_revision = "h9i0j1k2l3m4"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("notes", sa.Column("annotation_type", sa.String(length=50), nullable=True))


def downgrade() -> None:
    op.drop_column("notes", "annotation_type")


