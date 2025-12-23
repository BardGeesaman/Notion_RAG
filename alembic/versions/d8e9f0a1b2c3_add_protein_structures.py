"""Add protein_structures and structure_files tables.

Revision ID: d8e9f0a1b2c3
Revises: c3d4e5f6a7b9
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import ARRAY, UUID


revision = "d8e9f0a1b2c3"
down_revision = "c3d4e5f6a7b9"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "protein_structures",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=True),
        sa.Column("pdb_id", sa.String(length=10), nullable=True),
        sa.Column("alphafold_uniprot_id", sa.String(length=20), nullable=True),
        sa.Column("source", sa.String(length=50), nullable=False, server_default="unknown"),
        sa.Column("resolution", sa.Float(), nullable=True),
        sa.Column("method", sa.String(length=100), nullable=True),
        sa.Column("chain_ids", ARRAY(sa.String()), nullable=True),
        sa.Column("prep_status", sa.String(length=50), nullable=True, server_default="raw"),
        sa.Column("prep_log", sa.Text(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_protein_structures_feature_id", "protein_structures", ["feature_id"])
    op.create_index("ix_protein_structures_pdb_id", "protein_structures", ["pdb_id"])
    op.create_index(
        "ix_protein_structures_alphafold_uniprot_id",
        "protein_structures",
        ["alphafold_uniprot_id"],
    )

    op.create_table(
        "structure_files",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("structure_id", UUID(as_uuid=True), sa.ForeignKey("protein_structures.id"), nullable=False),
        sa.Column("file_type", sa.String(length=50), nullable=False),
        sa.Column("file_path", sa.String(length=500), nullable=False),
        sa.Column("file_size_bytes", sa.BigInteger(), nullable=True),
        sa.Column("md5_hash", sa.String(length=32), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_structure_files_structure_id", "structure_files", ["structure_id"])


def downgrade() -> None:
    op.drop_index("ix_structure_files_structure_id", table_name="structure_files")
    op.drop_table("structure_files")

    op.drop_index("ix_protein_structures_alphafold_uniprot_id", table_name="protein_structures")
    op.drop_index("ix_protein_structures_pdb_id", table_name="protein_structures")
    op.drop_index("ix_protein_structures_feature_id", table_name="protein_structures")
    op.drop_table("protein_structures")


