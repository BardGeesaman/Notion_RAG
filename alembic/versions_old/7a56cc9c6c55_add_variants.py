"""Add variant interpretation tables.

Revision ID: 7a56cc9c6c55
Revises: c1d2e3f4a5b6
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "7a56cc9c6c55"
down_revision = "c1d2e3f4a5b6"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "variant_sets",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("name", sa.String(length=500), nullable=False),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("source_file", sa.String(length=500), nullable=True),
        sa.Column("source_type", sa.String(length=100), nullable=True),
        sa.Column("n_variants", sa.Integer(), nullable=True),
        sa.Column("n_genes", sa.Integer(), nullable=True),
        sa.Column("status", sa.String(length=50), nullable=True),
        sa.Column("created_at", sa.DateTime(), nullable=False, server_default=sa.func.now()),
    )
    op.create_index("ix_variant_sets_name", "variant_sets", ["name"])
    op.create_index("ix_variant_sets_source_type", "variant_sets", ["source_type"])
    op.create_index("ix_variant_sets_status", "variant_sets", ["status"])

    op.create_table(
        "variants",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("variant_set_id", UUID(as_uuid=True), sa.ForeignKey("variant_sets.id"), nullable=False),
        sa.Column("chromosome", sa.String(length=20), nullable=False),
        sa.Column("position", sa.Integer(), nullable=False),
        sa.Column("ref_allele", sa.String(length=500), nullable=False),
        sa.Column("alt_allele", sa.String(length=500), nullable=False),
        sa.Column("rs_id", sa.String(length=50), nullable=True),
        sa.Column("hgvs_genomic", sa.String(length=500), nullable=True),
        sa.Column("hgvs_coding", sa.String(length=500), nullable=True),
        sa.Column("hgvs_protein", sa.String(length=500), nullable=True),
        sa.Column("gene_symbol", sa.String(length=100), nullable=True),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=True),
        sa.Column("consequence", sa.String(length=200), nullable=True),
        sa.Column("impact", sa.String(length=50), nullable=True),
        sa.Column("gnomad_af", sa.Float(), nullable=True),
        sa.UniqueConstraint(
            "variant_set_id",
            "chromosome",
            "position",
            "ref_allele",
            "alt_allele",
            name="uq_variant_set_locus",
        ),
    )
    op.create_index("ix_variants_variant_set_id", "variants", ["variant_set_id"])
    op.create_index("ix_variants_chromosome", "variants", ["chromosome"])
    op.create_index("ix_variants_position", "variants", ["position"])
    op.create_index("ix_variants_rs_id", "variants", ["rs_id"])
    op.create_index("ix_variants_hgvs_genomic", "variants", ["hgvs_genomic"])
    op.create_index("ix_variants_hgvs_coding", "variants", ["hgvs_coding"])
    op.create_index("ix_variants_hgvs_protein", "variants", ["hgvs_protein"])
    op.create_index("ix_variants_gene_symbol", "variants", ["gene_symbol"])
    op.create_index("ix_variants_feature_id", "variants", ["feature_id"])
    op.create_index("ix_variants_consequence", "variants", ["consequence"])
    op.create_index("ix_variants_impact", "variants", ["impact"])
    op.create_index("ix_variants_set_gene", "variants", ["variant_set_id", "gene_symbol"])

    op.create_table(
        "variant_annotations",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("variant_id", UUID(as_uuid=True), sa.ForeignKey("variants.id"), nullable=False),
        sa.Column("annotation_source", sa.String(length=100), nullable=False),
        sa.Column("clinvar_id", sa.String(length=100), nullable=True),
        sa.Column("clinical_significance", sa.String(length=200), nullable=True),
        sa.Column("review_status", sa.String(length=200), nullable=True),
        sa.Column("condition", sa.String(length=500), nullable=True),
        sa.Column("cadd_phred", sa.Float(), nullable=True),
        sa.Column("revel_score", sa.Float(), nullable=True),
        sa.Column("sift_prediction", sa.String(length=50), nullable=True),
        sa.Column("polyphen_prediction", sa.String(length=50), nullable=True),
    )
    op.create_index("ix_variant_annotations_variant_id", "variant_annotations", ["variant_id"])
    op.create_index("ix_variant_annotations_annotation_source", "variant_annotations", ["annotation_source"])
    op.create_index("ix_variant_annotations_clinvar_id", "variant_annotations", ["clinvar_id"])
    op.create_index("ix_variant_annotations_clinical_significance", "variant_annotations", ["clinical_significance"])
    op.create_index("ix_variant_annotations_variant_source", "variant_annotations", ["variant_id", "annotation_source"])

    op.create_table(
        "gene_burdens",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("variant_set_id", UUID(as_uuid=True), sa.ForeignKey("variant_sets.id"), nullable=False),
        sa.Column("gene_symbol", sa.String(length=100), nullable=True),
        sa.Column("feature_id", UUID(as_uuid=True), sa.ForeignKey("features.id"), nullable=True),
        sa.Column("n_variants", sa.Integer(), nullable=True),
        sa.Column("n_pathogenic", sa.Integer(), nullable=True),
        sa.Column("n_vus", sa.Integer(), nullable=True),
        sa.Column("n_benign", sa.Integer(), nullable=True),
        sa.Column("burden_score", sa.Float(), nullable=True),
        sa.UniqueConstraint("variant_set_id", "gene_symbol", name="uq_gene_burden_set_gene"),
    )
    op.create_index("ix_gene_burdens_variant_set_id", "gene_burdens", ["variant_set_id"])
    op.create_index("ix_gene_burdens_gene_symbol", "gene_burdens", ["gene_symbol"])
    op.create_index("ix_gene_burdens_feature_id", "gene_burdens", ["feature_id"])
    op.create_index("ix_gene_burdens_set_gene", "gene_burdens", ["variant_set_id", "gene_symbol"])


def downgrade() -> None:
    op.drop_index("ix_gene_burdens_set_gene", table_name="gene_burdens")
    op.drop_index("ix_gene_burdens_feature_id", table_name="gene_burdens")
    op.drop_index("ix_gene_burdens_gene_symbol", table_name="gene_burdens")
    op.drop_index("ix_gene_burdens_variant_set_id", table_name="gene_burdens")
    op.drop_table("gene_burdens")

    op.drop_index("ix_variant_annotations_variant_source", table_name="variant_annotations")
    op.drop_index("ix_variant_annotations_clinical_significance", table_name="variant_annotations")
    op.drop_index("ix_variant_annotations_clinvar_id", table_name="variant_annotations")
    op.drop_index("ix_variant_annotations_annotation_source", table_name="variant_annotations")
    op.drop_index("ix_variant_annotations_variant_id", table_name="variant_annotations")
    op.drop_table("variant_annotations")

    op.drop_index("ix_variants_set_gene", table_name="variants")
    op.drop_index("ix_variants_impact", table_name="variants")
    op.drop_index("ix_variants_consequence", table_name="variants")
    op.drop_index("ix_variants_feature_id", table_name="variants")
    op.drop_index("ix_variants_gene_symbol", table_name="variants")
    op.drop_index("ix_variants_hgvs_protein", table_name="variants")
    op.drop_index("ix_variants_hgvs_coding", table_name="variants")
    op.drop_index("ix_variants_hgvs_genomic", table_name="variants")
    op.drop_index("ix_variants_rs_id", table_name="variants")
    op.drop_index("ix_variants_position", table_name="variants")
    op.drop_index("ix_variants_chromosome", table_name="variants")
    op.drop_index("ix_variants_variant_set_id", table_name="variants")
    op.drop_table("variants")

    op.drop_index("ix_variant_sets_status", table_name="variant_sets")
    op.drop_index("ix_variant_sets_source_type", table_name="variant_sets")
    op.drop_index("ix_variant_sets_name", table_name="variant_sets")
    op.drop_table("variant_sets")


