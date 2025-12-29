"""Audit logging utilities."""
import hashlib
import json
import logging
from typing import Any, Dict, List, Optional, cast
from datetime import datetime

from sqlalchemy.orm import Session

from amprenta_rag.utils.uuid_utils import ensure_uuid
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
                user_id=cast(Any, ensure_uuid(user_id)) if user_id and user_id != "test" else None,
                username=username,
                action=action,
                entity_type=entity_type,
                entity_id=str(ensure_uuid(entity_id)) if entity_id else entity_id,
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


def compute_checksum(data: dict) -> str:
    """
    Compute SHA256 checksum of data dictionary.

    Args:
        data: Dictionary to hash

    Returns:
        SHA256 hex digest
    """
    json_str = json.dumps(data, sort_keys=True, default=str)
    return hashlib.sha256(json_str.encode()).hexdigest()


def log_data_change(
    db: Session,
    entity_type: str,
    entity_id: str,
    action: str,
    actor_id: Optional[str],
    old_data: Optional[Dict[str, Any]],
    new_data: Optional[Dict[str, Any]],
) -> None:
    """
    Log data change with checksums and diff.

    Args:
        db: Database session
        entity_type: Type of entity
        entity_id: Entity UUID
        action: Action performed
        actor_id: User ID who made the change
        old_data: Previous data state
        new_data: New data state
    """
    old_checksum = compute_checksum(old_data) if old_data else None
    new_checksum = compute_checksum(new_data) if new_data else None
    
    changes = {}
    if old_data and new_data:
        for key in set(list(old_data.keys()) + list(new_data.keys())):
            old_val = old_data.get(key)
            new_val = new_data.get(key)
            if old_val != new_val:
                changes[key] = {"old": old_val, "new": new_val}
    
    log = AuditLog(
        entity_type=entity_type,
        entity_id=entity_id,
        action=action,
        user_id=ensure_uuid(actor_id) if actor_id else None,
        old_checksum=old_checksum,
        new_checksum=new_checksum,
        changes=changes if changes else None,
    )
    db.add(log)
    db.commit()


def get_audit_trail(db: Session, entity_type: str, entity_id: str) -> List[AuditLog]:
    """
    Get audit trail for an entity.

    Args:
        db: Database session
        entity_type: Type of entity
        entity_id: Entity UUID

    Returns:
        List of AuditLog records
    """
    return (
        db.query(AuditLog)
        .filter(AuditLog.entity_type == entity_type, AuditLog.entity_id == entity_id)
        .order_by(AuditLog.created_at.desc())
        .all()
    )


def verify_integrity(
    db: Session,
    entity_type: str,
    entity_id: str,
    current_data: Dict[str, Any],
) -> bool:
    """
    Verify data integrity against audit log.

    Args:
        db: Database session
        entity_type: Type of entity
        entity_id: Entity UUID
        current_data: Current data state

    Returns:
        True if checksum matches latest audit record
    """
    latest = (
        db.query(AuditLog)
        .filter(AuditLog.entity_type == entity_type, AuditLog.entity_id == entity_id)
        .order_by(AuditLog.created_at.desc())
        .first()
    )
    
    if not latest or not latest.new_checksum:
        return True
    
    current_checksum = compute_checksum(current_data)
    return current_checksum == latest.new_checksum
