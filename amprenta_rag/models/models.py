"""Shared SQLAlchemy base and helpers for split model modules."""

from __future__ import annotations

import uuid

from amprenta_rag.database.base import Base


def generate_uuid() -> uuid.UUID:
    """Generate a new UUID4 value."""
    return uuid.uuid4()


__all__ = ["Base", "generate_uuid"]

