"""Add binding_sites table.

Revision ID: e9f0a1b2c3d4
Revises: d8e9f0a1b2c3
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "e9f0a1b2c3d4"
down_revision = "d8e9f0a1b2c3"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "binding_sites",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("structure_id", UUID(as_uuid=True), sa.ForeignKey("protein_structures.id"), nullable=False),
        sa.Column("pocket_rank", sa.Integer(), nullable=False),
        sa.Column("score", sa.Float(), nullable=True),
        sa.Column("volume", sa.Float(), nullable=True),
        sa.Column("center_x", sa.Float(), nullable=True),
        sa.Column("center_y", sa.Float(), nullable=True),
        sa.Column("center_z", sa.Float(), nullable=True),
        sa.Column("residues", sa.JSON(), nullable=True),
        sa.Column("pocket_pdb_path", sa.String(length=500), nullable=True),
        sa.Column("detection_method", sa.String(length=50), nullable=False, server_default="fpocket"),
        sa.Column("detected_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.UniqueConstraint("structure_id", "detection_method", "pocket_rank", name="uq_binding_site_rank"),
    )
    op.create_index("ix_binding_sites_structure_id", "binding_sites", ["structure_id"])


def downgrade() -> None:
    op.drop_index("ix_binding_sites_structure_id", table_name="binding_sites")
    op.drop_table("binding_sites")


