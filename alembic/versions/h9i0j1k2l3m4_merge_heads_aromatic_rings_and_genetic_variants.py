"""Merge heads: aromatic_rings and genetic_variants

Revision ID: h9i0j1k2l3m4
Revises: aa1b2c3d4e5f, g3b4c5d6e7f8
Create Date: 2025-12-13
"""


# revision identifiers, used by Alembic.
revision = "h9i0j1k2l3m4"
down_revision = ("aa1b2c3d4e5f", "g3b4c5d6e7f8")
branch_labels = None
depends_on = None


def upgrade() -> None:
    # This is a merge migration; no schema changes.
    pass


def downgrade() -> None:
    # This is a merge migration; no schema changes.
    pass


