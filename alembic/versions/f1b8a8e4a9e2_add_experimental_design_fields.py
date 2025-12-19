"""Add experimental design fields to experiments and datasets

Revision ID: f1b8a8e4a9e2
Revises: 8283ae19e3a9
Create Date: 2025-12-10

"""

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "f1b8a8e4a9e2"
down_revision = "8283ae19e3a9"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Add design fields to experiments table
    op.add_column("experiments", sa.Column("design_type", sa.String(length=50), nullable=True))
    op.add_column("experiments", sa.Column("design_metadata", sa.JSON(), nullable=True))
    op.add_column("experiments", sa.Column("sample_groups", sa.JSON(), nullable=True))

    # Add sample-level design fields to datasets table
    op.add_column("datasets", sa.Column("sample_group", sa.String(length=100), nullable=True))
    op.add_column("datasets", sa.Column("timepoint", sa.String(length=50), nullable=True))
    op.add_column("datasets", sa.Column("intervention", sa.String(length=200), nullable=True))
    op.add_column("datasets", sa.Column("dose", sa.String(length=50), nullable=True))
    op.add_column("datasets", sa.Column("replicate_id", sa.String(length=50), nullable=True))


def downgrade() -> None:
    # Remove from datasets
    op.drop_column("datasets", "replicate_id")
    op.drop_column("datasets", "dose")
    op.drop_column("datasets", "intervention")
    op.drop_column("datasets", "timepoint")
    op.drop_column("datasets", "sample_group")

    # Remove from experiments
    op.drop_column("experiments", "sample_groups")
    op.drop_column("experiments", "design_metadata")
    op.drop_column("experiments", "design_type")

