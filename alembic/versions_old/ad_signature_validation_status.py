"""add validation_status to signatures

Revision ID: ad_signature_validation_status
Revises: 0225cbc008f1
Create Date: 2025-12-17
"""
from alembic import op
import sqlalchemy as sa


revision = "ad_signature_validation_status"
down_revision = "0225cbc008f1"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("signatures", sa.Column("validation_status", sa.String(length=20), nullable=True))
    op.create_index("ix_signatures_validation_status", "signatures", ["validation_status"])


def downgrade():
    op.drop_index("ix_signatures_validation_status", table_name="signatures")
    op.drop_column("signatures", "validation_status")

