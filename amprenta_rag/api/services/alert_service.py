"""
RepositoryNotification service.

This module encapsulates database operations for user-facing alerts:
- listing alerts (optionally unread-only, optionally scoped to a user),
- marking a single alert as read,
- marking all alerts as read for a user.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from amprenta_rag.database.models import RepositoryNotification, RepositorySubscription
from amprenta_rag.database.session import db_session


def _to_dict(alert: RepositoryNotification) -> dict:
    """Convert a `RepositoryNotification` ORM instance into a JSON-serializable dict."""
    return {
        "id": alert.id,
        "subscription_id": alert.subscription_id,
        "dataset_id": alert.dataset_id,
        "is_read": alert.is_read,
        "created_at": alert.created_at,
    }


def list_alerts(user_id: Optional[UUID], unread_only: bool = True) -> List[dict]:
    """
    List alerts visible to a user.

    Args:
        user_id: If provided, filter alerts to subscriptions owned by this user.
        unread_only: If True, only return unread alerts.

    Returns:
        A list of alert dictionaries sorted by newest first.
    """
    with db_session() as db:
        query = db.query(RepositoryNotification).join(RepositorySubscription)
        if unread_only:
            query = query.filter(RepositoryNotification.is_read.is_(False))
        if user_id:
            query = query.filter(RepositorySubscription.user_id == user_id)
        alerts = query.order_by(RepositoryNotification.created_at.desc()).all()
        return [_to_dict(a) for a in alerts]


def mark_read(alert_id: UUID) -> bool:
    """Mark a single alert as read. Returns False if the alert does not exist."""
    with db_session() as db:
        alert = db.query(RepositoryNotification).filter(RepositoryNotification.id == alert_id).first()
        if not alert:
            return False
        alert.is_read = True
        db.commit()
        return True


def mark_all_read(user_id: Optional[UUID]) -> int:
    """Mark all unread alerts as read for the given user. Returns count updated."""
    with db_session() as db:
        query = db.query(RepositoryNotification).join(RepositorySubscription)
        if user_id:
            query = query.filter(RepositorySubscription.user_id == user_id)
        updated = query.filter(RepositoryNotification.is_read.is_(False)).update({RepositoryNotification.is_read: True}, synchronize_session=False)
        db.commit()
        return updated

