"""Workflow action handlers."""
from __future__ import annotations

from typing import Dict, Any

from amprenta_rag.utils.uuid_utils import ensure_uuid

from amprenta_rag.database.models import Notification, Note
from amprenta_rag.utils.validation import validate_experiment, validate_compound
from amprenta_rag.automation.engine import register_action
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def send_notification(config: Dict[str, Any], context: Dict[str, Any], db) -> Dict[str, Any]:
    """
    Send a notification action handler.

    Config:
        - user_id: UUID of user to notify (or "creator" to use context creator)
        - title: Notification title (can use {placeholders})
        - message: Notification message (can use {placeholders})
        - notification_type: info, success, warning, discovery (default: info)

    Context:
        - Any fields that can be used in placeholders

    Returns:
        Result dict with notification_id
    """
    user_id = config.get("user_id")
    if user_id == "creator":
        user_id = context.get("created_by_id") or context.get("user_id")

    if not user_id:
        raise ValueError("No user_id specified in config or context")

    title = config.get("title", "Workflow Notification")
    message = config.get("message", "")
    notification_type = config.get("notification_type", "info")

    # Replace placeholders in title and message
    title = title.format(**context)
    message = message.format(**context) if message else None

    notification = Notification(
        user_id=ensure_uuid(user_id),
        title=title,
        message=message,
        notification_type=notification_type,
    )

    db.add(notification)
    db.commit()
    db.refresh(notification)

    logger.info("[AUTOMATION] Sent notification to user %s", user_id)

    return {"notification_id": str(notification.id), "status": "sent"}


def add_note(config: Dict[str, Any], context: Dict[str, Any], db) -> Dict[str, Any]:
    """
    Add a note to an entity action handler.

    Config:
        - entity_type: Type of entity (experiment, compound, etc.)
        - entity_id: UUID of entity (or use context entity_id)
        - content: Note content (can use {placeholders})
        - user_id: User creating note (or "creator" to use context creator)

    Context:
        - entity_type, entity_id, created_by_id, etc.

    Returns:
        Result dict with note_id
    """
    entity_type = config.get("entity_type") or context.get("entity_type")
    entity_id = config.get("entity_id") or context.get("entity_id")
    content = config.get("content", "")
    user_id = config.get("user_id", "creator")

    if user_id == "creator":
        user_id = context.get("created_by_id") or context.get("user_id")

    if not entity_type or not entity_id:
        raise ValueError("entity_type and entity_id required")

    if not content:
        raise ValueError("content required")

    # Replace placeholders
    content = content.format(**context)

    note = Note(
        entity_type=entity_type,
        entity_id=ensure_uuid(entity_id),
        content=content,
        created_by_id=ensure_uuid(user_id) if user_id else None,
    )

    db.add(note)
    db.commit()
    db.refresh(note)

    logger.info("[AUTOMATION] Added note to %s:%s", entity_type, entity_id)

    return {"note_id": str(note.id), "status": "created"}


def run_validation(config: Dict[str, Any], context: Dict[str, Any], db) -> Dict[str, Any]:
    """
    Run validation on an entity action handler.

    Config:
        - entity_type: Type of entity to validate (experiment, compound)
        - entity_id: UUID of entity (or use context entity_id)

    Context:
        - entity_type, entity_id, etc.

    Returns:
        Result dict with validation issues
    """
    entity_type = config.get("entity_type") or context.get("entity_type")
    entity_id = config.get("entity_id") or context.get("entity_id")

    if not entity_type or not entity_id:
        raise ValueError("entity_type and entity_id required")

    issues = []

    if entity_type == "experiment":
        issues = validate_experiment(str(entity_id), db)
    elif entity_type == "compound":
        issues = validate_compound(str(entity_id), db)
    else:
        raise ValueError(f"Unsupported entity_type for validation: {entity_type}")

    error_count = sum(1 for i in issues if i.severity == "error")
    warning_count = sum(1 for i in issues if i.severity == "warning")

    logger.info("[AUTOMATION] Validated %s:%s - %d errors, %d warnings", entity_type, entity_id, error_count, warning_count)

    return {
        "status": "completed",
        "error_count": error_count,
        "warning_count": warning_count,
        "issues": [
            {
                "field": i.field,
                "issue": i.issue,
                "severity": i.severity,
            }
            for i in issues
        ],
    }


# Register all actions on import
register_action("send_notification", send_notification)
register_action("add_note", add_note)
register_action("run_validation", run_validation)
