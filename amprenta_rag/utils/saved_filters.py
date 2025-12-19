"""Saved filter utilities."""
from __future__ import annotations

from typing import List, Dict, Any
from uuid import UUID

from amprenta_rag.database.models import SavedFilter
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def save_filter(
    name: str,
    entity_type: str,
    filters: Dict[str, Any],
    user_id: str,
    db
) -> SavedFilter:
    """
    Save a filter preset.

    Args:
        name: Name for the saved filter
        entity_type: Type of entity (experiment, compound, etc.)
        filters: Filter criteria as dictionary
        user_id: UUID of the user saving the filter
        db: Database session

    Returns:
        Saved SavedFilter object
    """
    saved_filter = SavedFilter(
        name=name,
        entity_type=entity_type,
        filters=filters,
        user_id=UUID(user_id) if isinstance(user_id, str) else user_id,
    )
    db.add(saved_filter)
    db.commit()
    db.refresh(saved_filter)

    logger.info("[FILTERS] Saved filter '%s' for user %s", name, user_id)
    return saved_filter


def get_user_filters(
    entity_type: str,
    user_id: str,
    db
) -> List[SavedFilter]:
    """
    Get all saved filters for a user and entity type.

    Args:
        entity_type: Type of entity to filter
        user_id: UUID of the user
        db: Database session

    Returns:
        List of SavedFilter objects
    """
    filters = db.query(SavedFilter).filter(
        SavedFilter.entity_type == entity_type,
        SavedFilter.user_id == (UUID(user_id) if isinstance(user_id, str) else user_id)
    ).order_by(SavedFilter.created_at.desc()).all()

    return filters


def delete_filter(filter_id: str, db) -> bool:
    """
    Delete a saved filter.

    Args:
        filter_id: UUID of the filter to delete
        db: Database session

    Returns:
        True if deleted, False if not found
    """
    filter_obj = db.query(SavedFilter).filter(
        SavedFilter.id == (UUID(filter_id) if isinstance(filter_id, str) else filter_id)
    ).first()

    if filter_obj:
        db.delete(filter_obj)
        db.commit()
        logger.info("[FILTERS] Deleted filter %s", filter_id)
        return True

    return False
