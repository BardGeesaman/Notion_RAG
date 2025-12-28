"""Add docking_runs and docking_poses tables.

Revision ID: f0a1b2c3d4e5
Revises: e9f0a1b2c3d4
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "f0a1b2c3d4e5"
down_revision = "e9f0a1b2c3d4"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "docking_runs",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("structure_id", UUID(as_uuid=True), sa.ForeignKey("protein_structures.id"), nullable=False),
        sa.Column("binding_site_id", UUID(as_uuid=True), sa.ForeignKey("binding_sites.id"), nullable=True),
        sa.Column("center_x", sa.Float(), nullable=True),
        sa.Column("center_y", sa.Float(), nullable=True),
        sa.Column("center_z", sa.Float(), nullable=True),
        sa.Column("size_x", sa.Float(), nullable=True),
        sa.Column("size_y", sa.Float(), nullable=True),
        sa.Column("size_z", sa.Float(), nullable=True),
        sa.Column("status", sa.String(length=50), nullable=False, server_default="pending"),
        sa.Column("total_compounds", sa.Integer(), nullable=True),
        sa.Column("completed_compounds", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("started_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("completed_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("error_log", sa.Text(), nullable=True),
    )
    op.create_index("ix_docking_runs_structure_id", "docking_runs", ["structure_id"])
    op.create_index("ix_docking_runs_binding_site_id", "docking_runs", ["binding_site_id"])

    op.create_table(
        "docking_poses",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("docking_run_id", UUID(as_uuid=True), sa.ForeignKey("docking_runs.id"), nullable=False),
        sa.Column("compound_id", UUID(as_uuid=True), sa.ForeignKey("compounds.id"), nullable=False),
        sa.Column("binding_affinity", sa.Float(), nullable=True),
        sa.Column("pose_rank", sa.Integer(), nullable=True),
        sa.Column("rmsd_lb", sa.Float(), nullable=True),
        sa.Column("rmsd_ub", sa.Float(), nullable=True),
        sa.Column("pose_pdbqt_path", sa.String(length=500), nullable=True),
        sa.Column("docked_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_docking_poses_docking_run_id", "docking_poses", ["docking_run_id"])
    op.create_index("ix_docking_poses_compound_id", "docking_poses", ["compound_id"])
    op.create_index("ix_docking_poses_run_rank", "docking_poses", ["docking_run_id", "pose_rank"])


def downgrade() -> None:
    op.drop_index("ix_docking_poses_run_rank", table_name="docking_poses")
    op.drop_index("ix_docking_poses_compound_id", table_name="docking_poses")
    op.drop_index("ix_docking_poses_docking_run_id", table_name="docking_poses")
    op.drop_table("docking_poses")

    op.drop_index("ix_docking_runs_binding_site_id", table_name="docking_runs")
    op.drop_index("ix_docking_runs_structure_id", table_name="docking_runs")
    op.drop_table("docking_runs")


