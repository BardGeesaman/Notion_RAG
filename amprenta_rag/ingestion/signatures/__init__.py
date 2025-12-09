"""
Signature utilities.

This package provides utilities for signature handling.
Notion CRUD operations have been removed - Postgres is now the source of truth.
"""

from __future__ import annotations

from amprenta_rag.ingestion.signatures.short_id import generate_signature_short_id

__all__ = [
    "generate_signature_short_id",
]
