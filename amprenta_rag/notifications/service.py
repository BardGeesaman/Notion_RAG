"""Notification service for managing user notifications."""
from __future__ import annotations

from typing import List, Optional, Any, cast

from amprenta_rag.utils.uuid_utils import ensure_uuid

from amprenta_rag.database.models import Notification, User
from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_notification(
    user_id: str,
    title: str,
    message: Optional[str] = None,
    notification_type: str = "info",
    db=None,
) -> Notification:
    """
    Create a new notification for a user.

    Args:
        user_id: UUID of the user
        title: Notification title
        message: Optional notification message
        notification_type: Type of notification (info, success, warning, discovery)
        db: Database session

    Returns:
        Created Notification object
    """
    notification = Notification(
        user_id=cast(Any, ensure_uuid(user_id)),
        title=title,
        message=message,
        notification_type=notification_type,
        is_read=False,
    )

    db.add(notification)
    db.commit()
    db.refresh(notification)

    logger.info("[NOTIFICATION] Created notification %s for user %s", notification.id, user_id)
    return notification


def get_unread_notifications(user_id: str, db, limit: int = 10) -> List[Notification]:
    """
    Get unread notifications for a user.

    Args:
        user_id: UUID of the user
        db: Database session
        limit: Maximum number of notifications to return

    Returns:
        List of Notification objects, ordered by created_at descending
    """
    notifications = (
        db.query(Notification)
        .filter(
            Notification.user_id == ensure_uuid(user_id),
            Notification.is_read.is_(False),
        )
        .order_by(Notification.created_at.desc())
        .limit(limit)
        .all()
    )

    return notifications


def get_unread_count(user_id: str, db) -> int:
    """
    Get count of unread notifications for a user.

    Args:
        user_id: UUID of the user
        db: Database session

    Returns:
        Count of unread notifications
    """
    count = (
        db.query(Notification)
        .filter(
            Notification.user_id == ensure_uuid(user_id),
            Notification.is_read.is_(False),
        )
        .count()
    )

    return count


def mark_as_read(notification_id: str, db) -> None:
    """
    Mark a notification as read.

    Args:
        notification_id: UUID of the notification
        db: Database session
    """
    notification = db.query(Notification).filter(
        Notification.id == ensure_uuid(notification_id)
    ).first()

    if notification:
        notification.is_read = True
        db.commit()
        logger.debug("[NOTIFICATION] Marked notification %s as read", notification_id)
    else:
        logger.warning("[NOTIFICATION] Notification %s not found", notification_id)


def mark_all_read(user_id: str, db) -> None:
    """
    Mark all notifications as read for a user.

    Args:
        user_id: UUID of the user
        db: Database session
    """
    updated = (
        db.query(Notification)
        .filter(
            Notification.user_id == ensure_uuid(user_id),
            Notification.is_read.is_(False),
        )
        .update({"is_read": True})
    )

    db.commit()
    logger.info("[NOTIFICATION] Marked %d notifications as read for user %s", updated, user_id)


def create_admin_notification(
    title: str,
    message: Optional[str] = None,
    notification_type: str = "warning",
) -> List[Notification]:
    """
    Create notifications for all admin users.

    Args:
        title: Notification title
        message: Optional notification message
        notification_type: Type of notification (warning, critical, info)

    Returns:
        List of created Notification objects
    """
    notifications = []
    
    with db_session() as db:
        # Get all admin users
        admin_users = db.query(User).filter(User.role == "admin", User.is_active.is_(True)).all()
        
        if not admin_users:
            logger.warning("[NOTIFICATION] No admin users found for system notification")
            return notifications
        
        # Create notification for each admin
        for admin in admin_users:
            notification = Notification(
                user_id=admin.id,
                title=title,
                message=message,
                notification_type=notification_type,
                is_read=False,
            )
            
            db.add(notification)
            notifications.append(notification)
        
        db.commit()
        
        # Refresh all notifications to get IDs
        for notification in notifications:
            db.refresh(notification)
        
        logger.info(
            "[NOTIFICATION] Created %d admin notifications: %s",
            len(notifications),
            title
        )
    
    return notifications
