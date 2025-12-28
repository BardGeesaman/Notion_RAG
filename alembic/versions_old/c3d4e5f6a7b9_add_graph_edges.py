"""Add graph_edges table (generic evidence graph edges).

Revision ID: c3d4e5f6a7b9
Revises: f1a2b3c4d5e6
Create Date: 2025-12-23
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import UUID


revision = "c3d4e5f6a7b9"
down_revision = "f1a2b3c4d5e6"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "graph_edges",
        sa.Column("id", UUID(as_uuid=True), primary_key=True),
        sa.Column("source_entity_type", sa.String(length=50), nullable=False),
        sa.Column("source_entity_id", UUID(as_uuid=True), nullable=False),
        sa.Column("target_entity_type", sa.String(length=50), nullable=False),
        sa.Column("target_entity_id", UUID(as_uuid=True), nullable=False),
        sa.Column("relationship_type", sa.String(length=100), nullable=False),
        sa.Column("confidence", sa.Float(), nullable=True),
        sa.Column("evidence_source", sa.String(length=100), nullable=True),
        sa.Column("provenance", sa.JSON(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.UniqueConstraint(
            "source_entity_type",
            "source_entity_id",
            "target_entity_type",
            "target_entity_id",
            "relationship_type",
            "evidence_source",
            name="uq_graph_edge",
        ),
    )

    op.create_index("ix_graph_edges_source_entity_type", "graph_edges", ["source_entity_type"])
    op.create_index("ix_graph_edges_source_entity_id", "graph_edges", ["source_entity_id"])
    op.create_index("ix_graph_edges_target_entity_type", "graph_edges", ["target_entity_type"])
    op.create_index("ix_graph_edges_target_entity_id", "graph_edges", ["target_entity_id"])
    op.create_index("ix_graph_edges_relationship_type", "graph_edges", ["relationship_type"])
    op.create_index("ix_graph_edges_evidence_source", "graph_edges", ["evidence_source"])

    op.create_index(
        "ix_graph_edges_src_rel",
        "graph_edges",
        ["source_entity_type", "source_entity_id", "relationship_type"],
    )
    op.create_index(
        "ix_graph_edges_tgt_rel",
        "graph_edges",
        ["target_entity_type", "target_entity_id", "relationship_type"],
    )


def downgrade() -> None:
    op.drop_index("ix_graph_edges_tgt_rel", table_name="graph_edges")
    op.drop_index("ix_graph_edges_src_rel", table_name="graph_edges")

    op.drop_index("ix_graph_edges_evidence_source", table_name="graph_edges")
    op.drop_index("ix_graph_edges_relationship_type", table_name="graph_edges")
    op.drop_index("ix_graph_edges_target_entity_id", table_name="graph_edges")
    op.drop_index("ix_graph_edges_target_entity_type", table_name="graph_edges")
    op.drop_index("ix_graph_edges_source_entity_id", table_name="graph_edges")
    op.drop_index("ix_graph_edges_source_entity_type", table_name="graph_edges")

    op.drop_table("graph_edges")


