"""add mwtab_json to datasets

Revision ID: ac_dataset_mwtab_json
Revises: ab_report_schedules
Create Date: 2025-12-17
"""
from alembic import op
import sqlalchemy as sa


revision = "ac_dataset_mwtab_json"
down_revision = "ab_report_schedules"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("datasets", sa.Column("mwtab_json", sa.JSON(), nullable=True))


def downgrade():
    op.drop_column("datasets", "mwtab_json")

