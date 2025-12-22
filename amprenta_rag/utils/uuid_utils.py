"""UUID type guard utilities for type-safe UUID handling."""
from __future__ import annotations

from typing import overload
from uuid import UUID


@overload
def ensure_uuid(value: UUID) -> UUID: ...


@overload
def ensure_uuid(value: str) -> UUID: ...


@overload
def ensure_uuid(value: str | UUID) -> UUID: ...


@overload
def ensure_uuid(value: str | UUID | None) -> UUID | None: ...


def ensure_uuid(value: str | UUID | None) -> UUID | None:
    """Convert string to UUID if needed, pass through UUID, return None for None.
    
    This function provides type-safe UUID handling that mypy can understand
    through overloaded signatures.
    
    Args:
        value: A UUID, string representation of UUID, or None.
        
    Returns:
        UUID object if value is str or UUID, None if value is None.
        
    Raises:
        ValueError: If string is not a valid UUID format.
        
    Examples:
        >>> ensure_uuid("550e8400-e29b-41d4-a716-446655440000")
        UUID('550e8400-e29b-41d4-a716-446655440000')
        
        >>> ensure_uuid(UUID("550e8400-e29b-41d4-a716-446655440000"))
        UUID('550e8400-e29b-41d4-a716-446655440000')
        
        >>> ensure_uuid(None)
        None
    """
    if value is None:
        return None
    if isinstance(value, UUID):
        return value
    return UUID(value)

