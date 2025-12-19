"""Optimistic locking utilities for concurrent editing safety."""

from typing import Any, Dict


class ConflictError(Exception):
    """Raised when version conflict detected."""

    def __init__(self, entity_type: str, entity_id: str, expected: int, actual: int):
        self.entity_type = entity_type
        self.entity_id = entity_id
        self.expected = expected
        self.actual = actual
        super().__init__(f"Version conflict on {entity_type} {entity_id}: expected {expected}, found {actual}")


def check_version(entity: Any, expected_version: int) -> bool:
    """Check if entity version matches expected."""
    return getattr(entity, 'version', 1) == expected_version


def update_with_lock(entity: Any, updates: Dict[str, Any], expected_version: int, db: Any) -> Any:
    """Update entity with optimistic locking.

    Raises ConflictError if version mismatch.
    Increments version on successful update.
    """
    if not check_version(entity, expected_version):
        raise ConflictError(
            entity.__class__.__name__,
            str(entity.id),
            expected_version,
            entity.version
        )

    for key, value in updates.items():
        setattr(entity, key, value)
    entity.version += 1
    db.commit()
    return entity

