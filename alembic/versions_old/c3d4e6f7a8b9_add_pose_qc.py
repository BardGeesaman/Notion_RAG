"""Add pose_qualities and pose_interactions tables.

Revision ID: c3d4e6f7a8b9
Revises: b2c3d4e6f7a8
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "c3d4e6f7a8b9"
down_revision = "b2c3d4e6f7a8"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "pose_qualities",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("pose_id", UUID(as_uuid=True), sa.ForeignKey("docking_poses.id"), nullable=False),
        sa.Column("num_hbonds", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("num_hydrophobic", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("num_salt_bridges", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("num_pi_stacking", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("num_pi_cation", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("num_halogen", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("num_metal", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("total_interactions", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("has_clashes", sa.Boolean(), nullable=False, server_default=sa.text("false")),
        sa.Column("ligand_efficiency", sa.Float(), nullable=True),
        sa.Column("analyzed_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
        sa.UniqueConstraint("pose_id", name="uq_pose_quality_pose_id"),
    )
    op.create_index("ix_pose_qualities_pose_id", "pose_qualities", ["pose_id"])

    op.create_table(
        "pose_interactions",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("pose_id", UUID(as_uuid=True), sa.ForeignKey("docking_poses.id"), nullable=False),
        sa.Column("interaction_type", sa.String(length=50), nullable=False),
        sa.Column("ligand_atom", sa.String(length=200), nullable=True),
        sa.Column("protein_residue", sa.String(length=200), nullable=True),
        sa.Column("distance", sa.Float(), nullable=True),
        sa.Column("angle", sa.Float(), nullable=True),
    )
    op.create_index("ix_pose_interactions_pose_id", "pose_interactions", ["pose_id"])
    op.create_index(
        "ix_pose_interactions_pose_type",
        "pose_interactions",
        ["pose_id", "interaction_type"],
    )


def downgrade() -> None:
    op.drop_index("ix_pose_interactions_pose_type", table_name="pose_interactions")
    op.drop_index("ix_pose_interactions_pose_id", table_name="pose_interactions")
    op.drop_table("pose_interactions")

    op.drop_index("ix_pose_qualities_pose_id", table_name="pose_qualities")
    op.drop_table("pose_qualities")


