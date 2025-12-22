"""Notes utilities for managing user notes on entities."""
from __future__ import annotations

from typing import List, cast, Any
from uuid import UUID

from amprenta_rag.utils.uuid_utils import ensure_uuid

from amprenta_rag.database.models import Note
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def add_note(entity_type: str, entity_id: str, content: str, user_id: str, db) -> Note:
    """
    Add a note to an entity.

    Args:
        entity_type: Type of entity (experiment, compound, signature, etc.)
        entity_id: UUID of the entity
        content: Note content
        user_id: UUID of the user creating the note
        db: Database session

    Returns:
        Created Note object
    """
    note = Note(
        entity_type=entity_type,
        entity_id=cast(Any, ensure_uuid(entity_id)),
        content=content,
        created_by_id=cast(Any, ensure_uuid(user_id)) if user_id and user_id != "test" else None,
    )

    db.add(note)
    db.commit()
    db.refresh(note)

    logger.info("[NOTES] Added note %s for %s %s", note.id, entity_type, entity_id)
    return note


def get_notes(entity_type: str, entity_id: str, db) -> List[Note]:
    """
    Get all notes for an entity.

    Args:
        entity_type: Type of entity
        entity_id: UUID of the entity
        db: Database session

    Returns:
        List of Note objects, ordered by created_at descending
    """
    notes = (
        db.query(Note)
        .filter(
            Note.entity_type == entity_type,
            Note.entity_id == ensure_uuid(entity_id),
        )
        .order_by(Note.created_at.desc())
        .all()
    )

    return notes


def delete_note(note_id: str, db) -> bool:
    """
    Delete a note.

    Args:
        note_id: UUID of the note
        db: Database session

    Returns:
        True if deleted, False if not found
    """
    note = db.query(Note).filter(
        Note.id == ensure_uuid(note_id)
    ).first()

    if note:
        db.delete(note)
        db.commit()
        logger.info("[NOTES] Deleted note %s", note_id)
        return True

    return False
