"""Add pgvector embedding columns + index to rag_chunks.

Revision ID: a1b2c3
Revises: 59dd5c3fc6a8
Create Date: 2025-12-22
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op

try:
    from pgvector.sqlalchemy import Vector
except Exception:  # pragma: no cover
    Vector = None  # type: ignore


# revision identifiers, used by Alembic.
revision = "a1b2c3"
down_revision = "59dd5c3fc6a8"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Try to enable pgvector extension, but skip if not available
    try:
        op.execute("CREATE EXTENSION IF NOT EXISTS vector")
    except Exception as e:
        # pgvector not available - skip this migration
        import warnings
        warnings.warn(f"pgvector extension not available, skipping embedding migration: {e}")
        return

    if Vector is None:  # pragma: no cover
        import warnings
        warnings.warn("pgvector Python package not installed, skipping embedding migration")
        return

    op.add_column("rag_chunks", sa.Column("embedding", Vector(3072), nullable=True))
    op.add_column("rag_chunks", sa.Column("embedding_model", sa.String(length=100), nullable=True))

    # Create HNSW index concurrently (must run outside a transaction).
    try:
        with op.get_context().autocommit_block():
            op.execute(
                """
                CREATE INDEX CONCURRENTLY IF NOT EXISTS ix_rag_chunks_embedding_hnsw
                ON rag_chunks
                USING hnsw (embedding vector_cosine_ops)
                """
            )
    except Exception as e:
        import warnings
        warnings.warn(f"Failed to create HNSW index: {e}")


def downgrade() -> None:
    # Drop index concurrently (must run outside a transaction).
    with op.get_context().autocommit_block():
        op.execute("DROP INDEX CONCURRENTLY IF EXISTS ix_rag_chunks_embedding_hnsw")

    op.drop_column("rag_chunks", "embedding_model")
    op.drop_column("rag_chunks", "embedding")

    # Drop extension (requested for rollback completeness).
    op.execute("DROP EXTENSION IF EXISTS vector")


