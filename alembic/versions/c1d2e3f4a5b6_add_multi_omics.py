"""Add multi-omics latent factor tables.

Revision ID: c1d2e3f4a5b6
Revises: a3c4d5e6f7a8
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "c1d2e3f4a5b6"
down_revision = "a3c4d5e6f7a8"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "multi_omics_experiments",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(length=500), nullable=False),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("dataset_ids", sa.JSON(), nullable=True),
        sa.Column("sample_mapping", sa.JSON(), nullable=True),
        sa.Column("n_factors", sa.Integer(), nullable=True),
        sa.Column("convergence_mode", sa.String(length=100), nullable=True),
        sa.Column("status", sa.String(length=50), nullable=True),
        sa.Column("processing_log", sa.Text(), nullable=True),
        sa.Column("created_at", sa.DateTime(), nullable=False, server_default=sa.func.now()),
        sa.Column("processed_at", sa.DateTime(), nullable=True),
    )
    op.create_index("ix_multi_omics_experiments_name", "multi_omics_experiments", ["name"])
    op.create_index("ix_multi_omics_experiments_status", "multi_omics_experiments", ["status"])

    op.create_table(
        "latent_factors",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("experiment_id", UUID(as_uuid=True), sa.ForeignKey("multi_omics_experiments.id"), nullable=False),
        sa.Column("factor_index", sa.Integer(), nullable=False),
        sa.Column("variance_explained", sa.JSON(), nullable=True),
        sa.Column("description", sa.Text(), nullable=True),
        sa.UniqueConstraint("experiment_id", "factor_index", name="uq_latent_factor_experiment_index"),
    )
    op.create_index("ix_latent_factors_experiment_id", "latent_factors", ["experiment_id"])
    op.create_index("ix_latent_factors_factor_index", "latent_factors", ["factor_index"])

    op.create_table(
        "factor_loadings",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("factor_id", UUID(as_uuid=True), sa.ForeignKey("latent_factors.id"), nullable=False),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=False),
        sa.Column("loading", sa.Float(), nullable=False),
        sa.Column("abs_loading", sa.Float(), nullable=False),
        sa.Column("omics_type", sa.String(length=50), nullable=True),
        sa.UniqueConstraint("factor_id", "feature_id", name="uq_factor_loading_factor_feature"),
    )
    op.create_index("ix_factor_loadings_factor_id", "factor_loadings", ["factor_id"])
    op.create_index("ix_factor_loadings_feature_id", "factor_loadings", ["feature_id"])
    op.create_index("ix_factor_loadings_abs_loading", "factor_loadings", ["abs_loading"])
    op.create_index("ix_factor_loadings_omics_type", "factor_loadings", ["omics_type"])
    op.create_index("ix_factor_loadings_factor_abs", "factor_loadings", ["factor_id", "abs_loading"])

    op.create_table(
        "factor_scores",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("factor_id", UUID(as_uuid=True), sa.ForeignKey("latent_factors.id"), nullable=False),
        sa.Column("sample_id", sa.String(length=200), nullable=False),
        sa.Column("score", sa.Float(), nullable=False),
        sa.UniqueConstraint("factor_id", "sample_id", name="uq_factor_score_factor_sample"),
    )
    op.create_index("ix_factor_scores_factor_id", "factor_scores", ["factor_id"])
    op.create_index("ix_factor_scores_sample_id", "factor_scores", ["sample_id"])
    op.create_index("ix_factor_scores_factor_sample", "factor_scores", ["factor_id", "sample_id"])


def downgrade() -> None:
    op.drop_index("ix_factor_scores_factor_sample", table_name="factor_scores")
    op.drop_index("ix_factor_scores_sample_id", table_name="factor_scores")
    op.drop_index("ix_factor_scores_factor_id", table_name="factor_scores")
    op.drop_table("factor_scores")

    op.drop_index("ix_factor_loadings_factor_abs", table_name="factor_loadings")
    op.drop_index("ix_factor_loadings_omics_type", table_name="factor_loadings")
    op.drop_index("ix_factor_loadings_abs_loading", table_name="factor_loadings")
    op.drop_index("ix_factor_loadings_feature_id", table_name="factor_loadings")
    op.drop_index("ix_factor_loadings_factor_id", table_name="factor_loadings")
    op.drop_table("factor_loadings")

    op.drop_index("ix_latent_factors_factor_index", table_name="latent_factors")
    op.drop_index("ix_latent_factors_experiment_id", table_name="latent_factors")
    op.drop_table("latent_factors")

    op.drop_index("ix_multi_omics_experiments_status", table_name="multi_omics_experiments")
    op.drop_index("ix_multi_omics_experiments_name", table_name="multi_omics_experiments")
    op.drop_table("multi_omics_experiments")


