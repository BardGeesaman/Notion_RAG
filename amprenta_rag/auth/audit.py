"""Audit logging utilities."""
import logging
from typing import Optional, Any
from datetime import datetime

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import AuditLog

logger = logging.getLogger(__name__)


def log_action(
    action: str,
    user_id: Optional[str] = None,
    username: Optional[str] = None,
    entity_type: Optional[str] = None,
    entity_id: Optional[str] = None,
    details: Optional[dict] = None,
    ip_address: Optional[str] = None
) -> None:
    """
    Log an audit event.

    Args:
        action: The action performed (login, logout, create, update, delete, view)
        user_id: UUID of the user performing the action
        username: Username (cached for when user is deleted)
        entity_type: Type of entity affected (experiment, dataset, etc.)
        entity_id: ID of the entity affected
        details: Additional context as dict
        ip_address: Client IP address
    """
    try:
        db_gen = get_db()
        db = next(db_gen)
        try:
            audit_entry = AuditLog(
                user_id=user_id if user_id and user_id != "test" else None,
                username=username,
                action=action,
                entity_type=entity_type,
                entity_id=entity_id,
                details=details,
                ip_address=ip_address,
                timestamp=datetime.utcnow()
            )
            db.add(audit_entry)
            db.commit()
            logger.info(f"[AUDIT] {username or 'system'}:{action} {entity_type or ''}:{entity_id or ''}")
        finally:
            db_gen.close()
    except Exception as e:
        logger.warning(f"[AUDIT] Failed to log action: {e}")


def log_login(user_id: str, username: str, ip_address: Optional[str] = None) -> None:
    """Log a user login."""
    log_action("login", user_id, username, ip_address=ip_address)


def log_logout(user_id: str, username: str) -> None:
    """Log a user logout."""
    log_action("logout", user_id, username)


def log_create(user_id: str, username: str, entity_type: str, entity_id: str, details: Optional[dict] = None) -> None:
    """Log entity creation."""
    log_action("create", user_id, username, entity_type, entity_id, details)


def log_update(user_id: str, username: str, entity_type: str, entity_id: str, details: Optional[dict] = None) -> None:
    """Log entity update."""
    log_action("update", user_id, username, entity_type, entity_id, details)


def log_delete(user_id: str, username: str, entity_type: str, entity_id: str) -> None:
    """Log entity deletion."""
    log_action("delete", user_id, username, entity_type, entity_id)
