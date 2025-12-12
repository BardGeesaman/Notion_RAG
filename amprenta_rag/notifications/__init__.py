"""Notification service module."""
from __future__ import annotations

from amprenta_rag.notifications.service import (
    create_notification,
    get_unread_count,
    get_unread_notifications,
    mark_all_read,
    mark_as_read,
)

__all__ = [
    "create_notification",
    "get_unread_notifications",
    "get_unread_count",
    "mark_as_read",
    "mark_all_read",
]
