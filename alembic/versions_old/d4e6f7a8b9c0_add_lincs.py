"""Add LINCS tables for connectivity mapping.

Revision ID: d4e6f7a8b9c0
Revises: c3d4e6f7a8b9
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "d4e6f7a8b9c0"
down_revision = "c3d4e6f7a8b9"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "lincs_genes",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("entrez_id", sa.Integer(), nullable=False),
        sa.Column("gene_symbol", sa.String(length=50), nullable=True),
        sa.Column("gene_title", sa.String(length=500), nullable=True),
        sa.Column("is_landmark", sa.Boolean(), nullable=False, server_default=sa.text("false")),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=True),
        sa.UniqueConstraint("entrez_id", name="uq_lincs_genes_entrez_id"),
    )
    op.create_index("ix_lincs_genes_entrez_id", "lincs_genes", ["entrez_id"])
    op.create_index("ix_lincs_genes_gene_symbol", "lincs_genes", ["gene_symbol"])
    op.create_index("ix_lincs_genes_feature_id", "lincs_genes", ["feature_id"])

    op.create_table(
        "lincs_signatures",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("sig_id", sa.String(length=100), nullable=False),
        sa.Column("pert_iname", sa.String(length=200), nullable=True),
        sa.Column("pert_id", sa.String(length=100), nullable=True),
        sa.Column("pert_type", sa.String(length=100), nullable=True),
        sa.Column("cell_id", sa.String(length=100), nullable=True),
        sa.Column("pert_time", sa.Float(), nullable=True),
        sa.Column("pert_dose", sa.Float(), nullable=True),
        sa.Column("gene_expression", sa.JSON(), nullable=True),
        sa.Column("tas", sa.Float(), nullable=True),
        sa.Column("distil_id", sa.String(length=500), nullable=True),
        sa.Column("ingested_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
        sa.UniqueConstraint("sig_id", name="uq_lincs_signatures_sig_id"),
    )
    op.create_index("ix_lincs_signatures_sig_id", "lincs_signatures", ["sig_id"])
    op.create_index("ix_lincs_signatures_pert_iname", "lincs_signatures", ["pert_iname"])
    op.create_index("ix_lincs_signatures_pert_id", "lincs_signatures", ["pert_id"])
    op.create_index("ix_lincs_signatures_cell_id", "lincs_signatures", ["cell_id"])

    op.create_table(
        "connectivity_scores",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("query_signature_id", UUID(as_uuid=True), sa.ForeignKey("signatures.id"), nullable=False),
        sa.Column("lincs_signature_id", UUID(as_uuid=True), sa.ForeignKey("lincs_signatures.id"), nullable=False),
        sa.Column("score", sa.Float(), nullable=False),
        sa.Column("score_type", sa.String(length=50), nullable=False, server_default="connectivity"),
        sa.Column("p_value", sa.Float(), nullable=True),
        sa.Column("computed_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
        sa.UniqueConstraint(
            "query_signature_id",
            "lincs_signature_id",
            "score_type",
            name="uq_connectivity_score",
        ),
    )
    op.create_index("ix_connectivity_scores_query_signature_id", "connectivity_scores", ["query_signature_id"])
    op.create_index("ix_connectivity_scores_lincs_signature_id", "connectivity_scores", ["lincs_signature_id"])
    op.create_index(
        "ix_connectivity_scores_query_score_type",
        "connectivity_scores",
        ["query_signature_id", "score_type"],
    )


def downgrade() -> None:
    op.drop_index("ix_connectivity_scores_query_score_type", table_name="connectivity_scores")
    op.drop_index("ix_connectivity_scores_lincs_signature_id", table_name="connectivity_scores")
    op.drop_index("ix_connectivity_scores_query_signature_id", table_name="connectivity_scores")
    op.drop_table("connectivity_scores")

    op.drop_index("ix_lincs_signatures_cell_id", table_name="lincs_signatures")
    op.drop_index("ix_lincs_signatures_pert_id", table_name="lincs_signatures")
    op.drop_index("ix_lincs_signatures_pert_iname", table_name="lincs_signatures")
    op.drop_index("ix_lincs_signatures_sig_id", table_name="lincs_signatures")
    op.drop_table("lincs_signatures")

    op.drop_index("ix_lincs_genes_feature_id", table_name="lincs_genes")
    op.drop_index("ix_lincs_genes_gene_symbol", table_name="lincs_genes")
    op.drop_index("ix_lincs_genes_entrez_id", table_name="lincs_genes")
    op.drop_table("lincs_genes")


