"""Add external sync tables (sync_jobs/records/conflicts) and activity stores (ChEMBL/PubChem).

Revision ID: e3b4c5d6e7f8
Revises: e2a3b4c5d6f7
Create Date: 2025-12-24
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB, UUID


revision = "e3b4c5d6e7f8"
down_revision = "e2a3b4c5d6f7"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "sync_jobs",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("source", sa.String(length=50), nullable=False),
        sa.Column("sync_type", sa.String(length=50), nullable=False),
        sa.Column("status", sa.String(length=50), nullable=False, server_default="pending"),
        sa.Column("records_synced", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("records_updated", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("records_new", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("conflicts_detected", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("started_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("completed_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("error_log", sa.Text(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_sync_jobs_source", "sync_jobs", ["source"])
    op.create_index("ix_sync_jobs_status", "sync_jobs", ["status"])
    op.create_index("ix_sync_jobs_source_status", "sync_jobs", ["source", "status"])

    op.create_table(
        "sync_records",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("job_id", UUID(as_uuid=True), sa.ForeignKey("sync_jobs.id", ondelete="SET NULL"), nullable=True),
        sa.Column("source", sa.String(length=50), nullable=False),
        sa.Column("external_id", sa.String(length=200), nullable=False),
        sa.Column("entity_type", sa.String(length=50), nullable=False),
        sa.Column("entity_id", UUID(as_uuid=True), nullable=True),
        sa.Column("checksum", sa.String(length=32), nullable=False),
        sa.Column("synced_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("metadata", JSONB(), nullable=False),
        sa.UniqueConstraint("source", "external_id", name="uq_sync_record_source_external_id"),
    )
    op.create_index("ix_sync_records_job_id", "sync_records", ["job_id"])
    op.create_index("ix_sync_records_source", "sync_records", ["source"])
    op.create_index("ix_sync_records_external_id", "sync_records", ["external_id"])
    op.create_index("ix_sync_records_source_external_id", "sync_records", ["source", "external_id"])

    op.create_table(
        "sync_conflicts",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("record_id", UUID(as_uuid=True), sa.ForeignKey("sync_records.id", ondelete="CASCADE"), nullable=False),
        sa.Column("conflict_type", sa.String(length=50), nullable=False),
        sa.Column("local_value", JSONB(), nullable=False),
        sa.Column("external_value", JSONB(), nullable=False),
        sa.Column("resolution_status", sa.String(length=50), nullable=False, server_default="pending"),
        sa.Column("resolved_at", sa.DateTime(timezone=True), nullable=True),
    )
    op.create_index("ix_sync_conflicts_record_id", "sync_conflicts", ["record_id"])
    op.create_index("ix_sync_conflicts_resolution_status", "sync_conflicts", ["resolution_status"])
    op.create_index(
        "ix_sync_conflicts_record_status", "sync_conflicts", ["record_id", "resolution_status"]
    )

    op.create_table(
        "chembl_activities",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("compound_id", UUID(as_uuid=True), sa.ForeignKey("compounds.id", ondelete="SET NULL"), nullable=True),
        sa.Column("chembl_molecule_id", sa.String(length=50), nullable=False),
        sa.Column("chembl_assay_id", sa.String(length=50), nullable=False),
        sa.Column("activity_type", sa.String(length=50), nullable=True),
        sa.Column("value", sa.Float(), nullable=True),
        sa.Column("units", sa.String(length=50), nullable=True),
        sa.Column("relation", sa.String(length=10), nullable=True),
        sa.Column("target_chembl_id", sa.String(length=50), nullable=True),
        sa.Column("target_name", sa.String(length=500), nullable=True),
        sa.Column("assay_type", sa.String(length=50), nullable=True),
        sa.Column("synced_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_chembl_activities_compound_id", "chembl_activities", ["compound_id"])
    op.create_index("ix_chembl_activities_chembl_molecule_id", "chembl_activities", ["chembl_molecule_id"])
    op.create_index("ix_chembl_activities_chembl_assay_id", "chembl_activities", ["chembl_assay_id"])
    op.create_index("ix_chembl_activities_target_chembl_id", "chembl_activities", ["target_chembl_id"])
    op.create_index(
        "ix_chembl_activities_molecule_assay", "chembl_activities", ["chembl_molecule_id", "chembl_assay_id"]
    )

    op.create_table(
        "pubchem_bioassays",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("compound_id", UUID(as_uuid=True), sa.ForeignKey("compounds.id", ondelete="SET NULL"), nullable=True),
        sa.Column("pubchem_cid", sa.Integer(), nullable=False),
        sa.Column("aid", sa.Integer(), nullable=False),
        sa.Column("activity_outcome", sa.String(length=50), nullable=True),
        sa.Column("activity_score", sa.Float(), nullable=True),
        sa.Column("assay_name", sa.String(length=500), nullable=True),
        sa.Column("synced_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_pubchem_bioassays_compound_id", "pubchem_bioassays", ["compound_id"])
    op.create_index("ix_pubchem_bioassays_pubchem_cid", "pubchem_bioassays", ["pubchem_cid"])
    op.create_index("ix_pubchem_bioassays_aid", "pubchem_bioassays", ["aid"])
    op.create_index("ix_pubchem_bioassays_cid_aid", "pubchem_bioassays", ["pubchem_cid", "aid"])

    op.add_column("compounds", sa.Column("chembl_id", sa.String(length=50), nullable=True))
    op.add_column("compounds", sa.Column("pubchem_cid", sa.Integer(), nullable=True))
    op.add_column("compounds", sa.Column("last_synced_at", sa.DateTime(timezone=True), nullable=True))
    op.create_index("ix_compounds_chembl_id", "compounds", ["chembl_id"])
    op.create_index("ix_compounds_pubchem_cid", "compounds", ["pubchem_cid"])
    op.create_index("ix_compounds_last_synced_at", "compounds", ["last_synced_at"])


def downgrade() -> None:
    op.drop_index("ix_compounds_last_synced_at", table_name="compounds")
    op.drop_index("ix_compounds_pubchem_cid", table_name="compounds")
    op.drop_index("ix_compounds_chembl_id", table_name="compounds")
    op.drop_column("compounds", "last_synced_at")
    op.drop_column("compounds", "pubchem_cid")
    op.drop_column("compounds", "chembl_id")

    op.drop_index("ix_pubchem_bioassays_cid_aid", table_name="pubchem_bioassays")
    op.drop_index("ix_pubchem_bioassays_aid", table_name="pubchem_bioassays")
    op.drop_index("ix_pubchem_bioassays_pubchem_cid", table_name="pubchem_bioassays")
    op.drop_index("ix_pubchem_bioassays_compound_id", table_name="pubchem_bioassays")
    op.drop_table("pubchem_bioassays")

    op.drop_index("ix_chembl_activities_molecule_assay", table_name="chembl_activities")
    op.drop_index("ix_chembl_activities_target_chembl_id", table_name="chembl_activities")
    op.drop_index("ix_chembl_activities_chembl_assay_id", table_name="chembl_activities")
    op.drop_index("ix_chembl_activities_chembl_molecule_id", table_name="chembl_activities")
    op.drop_index("ix_chembl_activities_compound_id", table_name="chembl_activities")
    op.drop_table("chembl_activities")

    op.drop_index("ix_sync_conflicts_record_status", table_name="sync_conflicts")
    op.drop_index("ix_sync_conflicts_resolution_status", table_name="sync_conflicts")
    op.drop_index("ix_sync_conflicts_record_id", table_name="sync_conflicts")
    op.drop_table("sync_conflicts")

    op.drop_index("ix_sync_records_source_external_id", table_name="sync_records")
    op.drop_index("ix_sync_records_external_id", table_name="sync_records")
    op.drop_index("ix_sync_records_source", table_name="sync_records")
    op.drop_index("ix_sync_records_job_id", table_name="sync_records")
    op.drop_table("sync_records")

    op.drop_index("ix_sync_jobs_source_status", table_name="sync_jobs")
    op.drop_index("ix_sync_jobs_status", table_name="sync_jobs")
    op.drop_index("ix_sync_jobs_source", table_name="sync_jobs")
    op.drop_table("sync_jobs")


