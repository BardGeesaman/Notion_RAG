"""Add ADME, PK, and Toxicology tables.

Revision ID: i9d0e1f2g3h4
Revises: h8c9d0e1f2g3
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID, JSON

revision = "i9d0e1f2g3h4"
down_revision = "h8c9d0e1f2g3"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ADME results
    op.create_table(
        "adme_results",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("compound_id", UUID(as_uuid=True), sa.ForeignKey("compounds.id"), nullable=False),
        sa.Column("assay_type", sa.String(50), nullable=False),
        sa.Column("value", sa.Float, nullable=False),
        sa.Column("unit", sa.String(20), nullable=True),
        sa.Column("conditions", JSON, nullable=True),
        sa.Column("created_at", sa.DateTime, nullable=False),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
    )
    op.create_index("ix_adme_results_compound_id", "adme_results", ["compound_id"])
    
    # PK studies
    op.create_table(
        "pk_studies",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("compound_id", UUID(as_uuid=True), sa.ForeignKey("compounds.id"), nullable=False),
        sa.Column("species", sa.String(50), nullable=False),
        sa.Column("route", sa.String(50), nullable=True),
        sa.Column("dose", sa.Float, nullable=True),
        sa.Column("dose_unit", sa.String(20), nullable=True),
        sa.Column("auc", sa.Float, nullable=True),
        sa.Column("c_max", sa.Float, nullable=True),
        sa.Column("t_max", sa.Float, nullable=True),
        sa.Column("half_life", sa.Float, nullable=True),
        sa.Column("bioavailability", sa.Float, nullable=True),
        sa.Column("clearance", sa.Float, nullable=True),
        sa.Column("vd", sa.Float, nullable=True),
        sa.Column("created_at", sa.DateTime, nullable=False),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
    )
    op.create_index("ix_pk_studies_compound_id", "pk_studies", ["compound_id"])
    
    # Toxicology results
    op.create_table(
        "toxicology_results",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("compound_id", UUID(as_uuid=True), sa.ForeignKey("compounds.id"), nullable=False),
        sa.Column("assay_type", sa.String(50), nullable=False),
        sa.Column("value", sa.Float, nullable=True),
        sa.Column("unit", sa.String(20), nullable=True),
        sa.Column("result", sa.String(50), nullable=True),
        sa.Column("conditions", JSON, nullable=True),
        sa.Column("created_at", sa.DateTime, nullable=False),
        sa.Column("created_by_id", UUID(as_uuid=True), sa.ForeignKey("users.id"), nullable=True),
    )
    op.create_index("ix_toxicology_results_compound_id", "toxicology_results", ["compound_id"])


def downgrade() -> None:
    op.drop_index("ix_toxicology_results_compound_id")
    op.drop_table("toxicology_results")
    op.drop_index("ix_pk_studies_compound_id")
    op.drop_table("pk_studies")
    op.drop_index("ix_adme_results_compound_id")
    op.drop_table("adme_results")
