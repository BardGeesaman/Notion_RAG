"""Add lipidomics spectral matching tables.

Revision ID: f7a8b9c0d1e2
Revises: e6f7a8b9c0d1
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "f7a8b9c0d1e2"
down_revision = "e6f7a8b9c0d1"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "spectral_libraries",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(length=200), nullable=False),
        sa.Column("version", sa.String(length=100), nullable=True),
        sa.Column("source_url", sa.String(length=500), nullable=True),
        sa.Column("file_path", sa.String(length=500), nullable=True),
        sa.Column("n_spectra", sa.Integer(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
        sa.UniqueConstraint("name", "version", name="uq_spectral_library_name_version"),
    )
    op.create_index("ix_spectral_libraries_name", "spectral_libraries", ["name"])
    op.create_index("ix_spectral_libraries_version", "spectral_libraries", ["version"])

    op.create_table(
        "spectral_references",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("library_id", UUID(as_uuid=True), sa.ForeignKey("spectral_libraries.id"), nullable=False),
        sa.Column("lipid_name", sa.String(length=500), nullable=False),
        sa.Column("lipid_class", sa.String(length=200), nullable=True),
        sa.Column("smiles", sa.Text(), nullable=True),
        sa.Column("inchi_key", sa.String(length=50), nullable=True),
        sa.Column("precursor_mz", sa.Float(), nullable=False),
        sa.Column("precursor_type", sa.String(length=50), nullable=True),
        sa.Column("collision_energy", sa.Float(), nullable=True),
        sa.Column("spectrum", sa.JSON(), nullable=False),
    )
    op.create_index("ix_spectral_references_library_id", "spectral_references", ["library_id"])
    op.create_index("ix_spectral_references_lipid_name", "spectral_references", ["lipid_name"])
    op.create_index("ix_spectral_references_lipid_class", "spectral_references", ["lipid_class"])
    op.create_index("ix_spectral_references_inchi_key", "spectral_references", ["inchi_key"])
    op.create_index("ix_spectral_references_precursor_mz", "spectral_references", ["precursor_mz"])
    op.create_index(
        "ix_spectral_ref_library_precursor",
        "spectral_references",
        ["library_id", "precursor_mz"],
    )

    op.create_table(
        "lipid_annotations",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=False),
        sa.Column("spectral_reference_id", UUID(as_uuid=True), sa.ForeignKey("spectral_references.id"), nullable=False),
        sa.Column("spectral_score", sa.Float(), nullable=False),
        sa.Column("mz_error_ppm", sa.Float(), nullable=True),
        sa.Column("matched_peaks", sa.Integer(), nullable=True),
        sa.Column("total_peaks", sa.Integer(), nullable=True),
        sa.Column("is_confident", sa.Boolean(), nullable=False, server_default=sa.text("false")),
        sa.Column("is_ambiguous", sa.Boolean(), nullable=False, server_default=sa.text("false")),
        sa.Column("rank", sa.Integer(), nullable=True),
        sa.Column("manually_reviewed", sa.Boolean(), nullable=False, server_default=sa.text("false")),
        sa.Column("review_status", sa.String(length=50), nullable=True),
        sa.Column("reviewed_at", sa.DateTime(timezone=True), nullable=True),
    )
    op.create_index("ix_lipid_annotations_feature_id", "lipid_annotations", ["feature_id"])
    op.create_index("ix_lipid_annotations_spectral_reference_id", "lipid_annotations", ["spectral_reference_id"])
    op.create_index("ix_lipid_annotations_feature_rank", "lipid_annotations", ["feature_id", "rank"])


def downgrade() -> None:
    op.drop_index("ix_lipid_annotations_feature_rank", table_name="lipid_annotations")
    op.drop_index("ix_lipid_annotations_spectral_reference_id", table_name="lipid_annotations")
    op.drop_index("ix_lipid_annotations_feature_id", table_name="lipid_annotations")
    op.drop_table("lipid_annotations")

    op.drop_index("ix_spectral_ref_library_precursor", table_name="spectral_references")
    op.drop_index("ix_spectral_references_precursor_mz", table_name="spectral_references")
    op.drop_index("ix_spectral_references_inchi_key", table_name="spectral_references")
    op.drop_index("ix_spectral_references_lipid_class", table_name="spectral_references")
    op.drop_index("ix_spectral_references_lipid_name", table_name="spectral_references")
    op.drop_index("ix_spectral_references_library_id", table_name="spectral_references")
    op.drop_table("spectral_references")

    op.drop_index("ix_spectral_libraries_version", table_name="spectral_libraries")
    op.drop_index("ix_spectral_libraries_name", table_name="spectral_libraries")
    op.drop_table("spectral_libraries")


