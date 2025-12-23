"""Add single-cell omics tables (Phase 1).

Revision ID: e6f7a8b9c0d1
Revises: d4e6f7a8b9c0
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "e6f7a8b9c0d1"
down_revision = "d4e6f7a8b9c0"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "single_cell_datasets",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("dataset_id", UUID(as_uuid=True), sa.ForeignKey("datasets.id"), nullable=False),
        sa.Column("h5ad_path", sa.String(length=500), nullable=False),
        sa.Column("file_size_bytes", sa.BigInteger(), nullable=True),
        sa.Column("n_cells", sa.Integer(), nullable=True),
        sa.Column("n_genes", sa.Integer(), nullable=True),
        sa.Column("normalization_method", sa.String(length=100), nullable=True),
        sa.Column("hvg_count", sa.Integer(), nullable=True),
        sa.Column("clustering_method", sa.String(length=100), nullable=True),
        sa.Column("clustering_resolution", sa.Float(), nullable=True),
        sa.Column("has_pca", sa.Boolean(), nullable=False, server_default=sa.text("false")),
        sa.Column("has_umap", sa.Boolean(), nullable=False, server_default=sa.text("false")),
        sa.Column("processing_status", sa.String(length=50), nullable=False, server_default="pending"),
        sa.Column("processing_log", sa.Text(), nullable=True),
        sa.Column("ingested_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
        sa.Column("processed_at", sa.DateTime(timezone=True), nullable=True),
    )
    op.create_index("ix_single_cell_datasets_dataset_id", "single_cell_datasets", ["dataset_id"])

    op.create_table(
        "cell_annotations",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column(
            "single_cell_dataset_id",
            UUID(as_uuid=True),
            sa.ForeignKey("single_cell_datasets.id"),
            nullable=False,
        ),
        sa.Column("barcode", sa.String(length=200), nullable=False),
        sa.Column("cluster_id", sa.Integer(), nullable=True),
        sa.Column("cell_type", sa.String(length=200), nullable=True),
        sa.Column("umap_1", sa.Float(), nullable=True),
        sa.Column("umap_2", sa.Float(), nullable=True),
        sa.Column("n_genes_detected", sa.Integer(), nullable=True),
        sa.Column("total_counts", sa.Float(), nullable=True),
        sa.Column("pct_mitochondrial", sa.Float(), nullable=True),
        sa.UniqueConstraint("single_cell_dataset_id", "barcode", name="uq_cell_annotation_barcode"),
    )
    op.create_index("ix_cell_annotations_single_cell_dataset_id", "cell_annotations", ["single_cell_dataset_id"])
    op.create_index("ix_cell_annotations_dataset_cluster", "cell_annotations", ["single_cell_dataset_id", "cluster_id"])

    op.create_table(
        "cell_clusters",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column(
            "single_cell_dataset_id",
            UUID(as_uuid=True),
            sa.ForeignKey("single_cell_datasets.id"),
            nullable=False,
        ),
        sa.Column("cluster_id", sa.Integer(), nullable=False),
        sa.Column("cell_type", sa.String(length=200), nullable=True),
        sa.Column("n_cells", sa.Integer(), nullable=True),
        sa.Column("marker_feature_ids", sa.JSON(), nullable=True),
        sa.Column("description", sa.Text(), nullable=True),
        sa.UniqueConstraint("single_cell_dataset_id", "cluster_id", name="uq_cell_cluster"),
    )
    op.create_index("ix_cell_clusters_single_cell_dataset_id", "cell_clusters", ["single_cell_dataset_id"])

    op.create_table(
        "cell_type_markers",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column(
            "single_cell_dataset_id",
            UUID(as_uuid=True),
            sa.ForeignKey("single_cell_datasets.id"),
            nullable=False,
        ),
        sa.Column("cluster_id", sa.Integer(), nullable=False),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=True),
        sa.Column("gene_symbol", sa.String(length=50), nullable=True),
        sa.Column("log2_fold_change", sa.Float(), nullable=True),
        sa.Column("pval", sa.Float(), nullable=True),
        sa.Column("pval_adj", sa.Float(), nullable=True),
        sa.Column("pct_in_cluster", sa.Float(), nullable=True),
        sa.Column("pct_out_cluster", sa.Float(), nullable=True),
    )
    op.create_index("ix_cell_type_markers_single_cell_dataset_id", "cell_type_markers", ["single_cell_dataset_id"])
    op.create_index("ix_cell_type_markers_feature_id", "cell_type_markers", ["feature_id"])
    op.create_index("ix_cell_type_markers_gene_symbol", "cell_type_markers", ["gene_symbol"])
    op.create_index(
        "ix_cell_type_markers_dataset_cluster",
        "cell_type_markers",
        ["single_cell_dataset_id", "cluster_id"],
    )


def downgrade() -> None:
    op.drop_index("ix_cell_type_markers_dataset_cluster", table_name="cell_type_markers")
    op.drop_index("ix_cell_type_markers_gene_symbol", table_name="cell_type_markers")
    op.drop_index("ix_cell_type_markers_feature_id", table_name="cell_type_markers")
    op.drop_index("ix_cell_type_markers_single_cell_dataset_id", table_name="cell_type_markers")
    op.drop_table("cell_type_markers")

    op.drop_index("ix_cell_clusters_single_cell_dataset_id", table_name="cell_clusters")
    op.drop_table("cell_clusters")

    op.drop_index("ix_cell_annotations_dataset_cluster", table_name="cell_annotations")
    op.drop_index("ix_cell_annotations_single_cell_dataset_id", table_name="cell_annotations")
    op.drop_table("cell_annotations")

    op.drop_index("ix_single_cell_datasets_dataset_id", table_name="single_cell_datasets")
    op.drop_table("single_cell_datasets")


