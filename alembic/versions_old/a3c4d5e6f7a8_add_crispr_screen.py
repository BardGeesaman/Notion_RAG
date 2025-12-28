"""Add CRISPR screen tables.

Revision ID: a3c4d5e6f7a8
Revises: f7a8b9c0d1e2
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "a3c4d5e6f7a8"
down_revision = "f7a8b9c0d1e2"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "crispr_screens",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("dataset_id", UUID(as_uuid=True), sa.ForeignKey("datasets.id"), nullable=False),
        sa.Column("name", sa.String(length=500), nullable=False),
        sa.Column("library_type", sa.String(length=100), nullable=True),
        sa.Column("cell_line", sa.String(length=200), nullable=True),
        sa.Column("treatment", sa.String(length=200), nullable=True),
        sa.Column("control_label", sa.String(length=200), nullable=True),
        sa.Column("treatment_label", sa.String(length=200), nullable=True),
        sa.Column("status", sa.String(length=50), nullable=True),
        sa.Column("created_at", sa.DateTime(), nullable=False, server_default=sa.func.now()),
    )
    op.create_index("ix_crispr_screens_dataset_id", "crispr_screens", ["dataset_id"])
    op.create_index("ix_crispr_screens_name", "crispr_screens", ["name"])
    op.create_index("ix_crispr_screens_library_type", "crispr_screens", ["library_type"])
    op.create_index("ix_crispr_screens_cell_line", "crispr_screens", ["cell_line"])
    op.create_index("ix_crispr_screens_treatment", "crispr_screens", ["treatment"])
    op.create_index("ix_crispr_screens_status", "crispr_screens", ["status"])

    op.create_table(
        "crispr_guides",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("screen_id", UUID(as_uuid=True), sa.ForeignKey("crispr_screens.id"), nullable=False),
        sa.Column("guide_seq", sa.String(length=200), nullable=False),
        sa.Column("gene_symbol", sa.String(length=100), nullable=True),
        sa.Column("gene_id", sa.String(length=100), nullable=True),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=True),
        sa.UniqueConstraint("screen_id", "guide_seq", name="uq_crispr_guide_screen_seq"),
    )
    op.create_index("ix_crispr_guides_screen_id", "crispr_guides", ["screen_id"])
    op.create_index("ix_crispr_guides_guide_seq", "crispr_guides", ["guide_seq"])
    op.create_index("ix_crispr_guides_gene_symbol", "crispr_guides", ["gene_symbol"])
    op.create_index("ix_crispr_guides_gene_id", "crispr_guides", ["gene_id"])
    op.create_index("ix_crispr_guides_feature_id", "crispr_guides", ["feature_id"])
    op.create_index("ix_crispr_guides_screen_gene", "crispr_guides", ["screen_id", "gene_symbol"])

    op.create_table(
        "crispr_results",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("screen_id", UUID(as_uuid=True), sa.ForeignKey("crispr_screens.id"), nullable=False),
        sa.Column("gene_symbol", sa.String(length=100), nullable=True),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=True),
        sa.Column("beta_score", sa.Float(), nullable=True),
        sa.Column("p_value", sa.Float(), nullable=True),
        sa.Column("fdr", sa.Float(), nullable=True),
        sa.Column("neg_lfc", sa.Float(), nullable=True),
        sa.Column("pos_lfc", sa.Float(), nullable=True),
        sa.Column("rank", sa.Integer(), nullable=True),
        sa.Column("is_hit", sa.Boolean(), nullable=False, server_default=sa.text("false")),
    )
    op.create_index("ix_crispr_results_screen_id", "crispr_results", ["screen_id"])
    op.create_index("ix_crispr_results_gene_symbol", "crispr_results", ["gene_symbol"])
    op.create_index("ix_crispr_results_feature_id", "crispr_results", ["feature_id"])
    op.create_index("ix_crispr_results_rank", "crispr_results", ["rank"])
    op.create_index("ix_crispr_results_is_hit", "crispr_results", ["is_hit"])
    op.create_index("ix_crispr_results_screen_gene", "crispr_results", ["screen_id", "gene_symbol"])


def downgrade() -> None:
    op.drop_index("ix_crispr_results_screen_gene", table_name="crispr_results")
    op.drop_index("ix_crispr_results_is_hit", table_name="crispr_results")
    op.drop_index("ix_crispr_results_rank", table_name="crispr_results")
    op.drop_index("ix_crispr_results_feature_id", table_name="crispr_results")
    op.drop_index("ix_crispr_results_gene_symbol", table_name="crispr_results")
    op.drop_index("ix_crispr_results_screen_id", table_name="crispr_results")
    op.drop_table("crispr_results")

    op.drop_index("ix_crispr_guides_screen_gene", table_name="crispr_guides")
    op.drop_index("ix_crispr_guides_feature_id", table_name="crispr_guides")
    op.drop_index("ix_crispr_guides_gene_id", table_name="crispr_guides")
    op.drop_index("ix_crispr_guides_gene_symbol", table_name="crispr_guides")
    op.drop_index("ix_crispr_guides_guide_seq", table_name="crispr_guides")
    op.drop_index("ix_crispr_guides_screen_id", table_name="crispr_guides")
    op.drop_table("crispr_guides")

    op.drop_index("ix_crispr_screens_status", table_name="crispr_screens")
    op.drop_index("ix_crispr_screens_treatment", table_name="crispr_screens")
    op.drop_index("ix_crispr_screens_cell_line", table_name="crispr_screens")
    op.drop_index("ix_crispr_screens_library_type", table_name="crispr_screens")
    op.drop_index("ix_crispr_screens_name", table_name="crispr_screens")
    op.drop_index("ix_crispr_screens_dataset_id", table_name="crispr_screens")
    op.drop_table("crispr_screens")


