"""Add search_vector to rag_chunks for BM25 search.

Revision ID: k1f2g3h4i5j6
Revises: j0e1f2g3h4i5
Create Date: 2025-12-11
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import TSVECTOR

revision = "k1f2g3h4i5j6"
down_revision = "j0e1f2g3h4i5"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("rag_chunks", sa.Column("search_vector", TSVECTOR, nullable=True))
    op.create_index(
        "ix_rag_chunks_search_vector",
        "rag_chunks",
        ["search_vector"],
        postgresql_using="gin"
    )


def downgrade() -> None:
    op.drop_index("ix_rag_chunks_search_vector")
    op.drop_column("rag_chunks", "search_vector")
