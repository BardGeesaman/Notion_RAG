"""Bookmarks service for managing user bookmarks."""
from __future__ import annotations

from typing import List
from uuid import UUID

from amprenta_rag.database.models import Bookmark
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def add_bookmark(user_id: str, entity_type: str, entity_id: str, db) -> Bookmark:
    """
    Add a bookmark for a user.

    Args:
        user_id: UUID of the user
        entity_type: Type of entity (experiment, compound, signature)
        entity_id: UUID of the entity
        db: Database session

    Returns:
        Created Bookmark object
    """
    bookmark = Bookmark(
        user_id=UUID(user_id) if isinstance(user_id, str) else user_id,
        entity_type=entity_type,
        entity_id=UUID(entity_id) if isinstance(entity_id, str) else entity_id,
    )

    db.add(bookmark)
    db.commit()
    db.refresh(bookmark)

    logger.info("[BOOKMARK] Added bookmark %s for user %s", bookmark.id, user_id)
    return bookmark


def remove_bookmark(user_id: str, entity_type: str, entity_id: str, db) -> None:
    """
    Remove a bookmark for a user.

    Args:
        user_id: UUID of the user
        entity_type: Type of entity
        entity_id: UUID of the entity
        db: Database session
    """
    bookmark = (
        db.query(Bookmark)
        .filter(
            Bookmark.user_id == (UUID(user_id) if isinstance(user_id, str) else user_id),
            Bookmark.entity_type == entity_type,
            Bookmark.entity_id == (UUID(entity_id) if isinstance(entity_id, str) else entity_id),
        )
        .first()
    )

    if bookmark:
        db.delete(bookmark)
        db.commit()
        logger.info("[BOOKMARK] Removed bookmark for user %s", user_id)
    else:
        logger.warning("[BOOKMARK] Bookmark not found for removal")


def get_user_bookmarks(user_id: str, db) -> List[Bookmark]:
    """
    Get all bookmarks for a user.

    Args:
        user_id: UUID of the user
        db: Database session

    Returns:
        List of Bookmark objects, ordered by created_at descending
    """
    bookmarks = (
        db.query(Bookmark)
        .filter(Bookmark.user_id == (UUID(user_id) if isinstance(user_id, str) else user_id))
        .order_by(Bookmark.created_at.desc())
        .all()
    )

    return bookmarks


def is_bookmarked(user_id: str, entity_type: str, entity_id: str, db) -> bool:
    """
    Check if an entity is bookmarked by a user.

    Args:
        user_id: UUID of the user
        entity_type: Type of entity
        entity_id: UUID of the entity
        db: Database session

    Returns:
        True if bookmarked, False otherwise
    """
    bookmark = (
        db.query(Bookmark)
        .filter(
            Bookmark.user_id == (UUID(user_id) if isinstance(user_id, str) else user_id),
            Bookmark.entity_type == entity_type,
            Bookmark.entity_id == (UUID(entity_id) if isinstance(entity_id, str) else entity_id),
        )
        .first()
    )

    return bookmark is not None
