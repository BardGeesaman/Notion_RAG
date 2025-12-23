"""Add phenotype_gene_associations table (HPO phenotype -> gene mapping).

Revision ID: b2c3d4e5f6a7
Revises: f1a2b3c4d5e6
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "b2c3d4e5f6a7"
down_revision = "f1a2b3c4d5e6"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "phenotype_gene_associations",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("hpo_id", sa.String(length=20), nullable=False),
        sa.Column("hpo_name", sa.String(length=500), nullable=True),
        sa.Column("gene_symbol", sa.String(length=50), nullable=False),
        sa.Column("ncbi_gene_id", sa.String(length=30), nullable=True),
        sa.Column("disease_id", sa.String(length=50), nullable=True),
        sa.Column("frequency", sa.String(length=50), nullable=True),
        sa.Column("source", sa.String(length=100), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.UniqueConstraint("hpo_id", "gene_symbol", "disease_id", name="uq_hpo_gene_disease"),
    )
    op.create_index("ix_phenotype_gene_hpo_id", "phenotype_gene_associations", ["hpo_id"])
    op.create_index("ix_phenotype_gene_symbol", "phenotype_gene_associations", ["gene_symbol"])


def downgrade() -> None:
    op.drop_index("ix_phenotype_gene_symbol", table_name="phenotype_gene_associations")
    op.drop_index("ix_phenotype_gene_hpo_id", table_name="phenotype_gene_associations")
    op.drop_table("phenotype_gene_associations")


