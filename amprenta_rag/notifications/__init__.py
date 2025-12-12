"""Notification service module."""
from __future__ import annotations

from amprenta_rag.notifications.service import (
    create_notification,
    get_unread_count,
    get_unread_notifications,
    mark_all_read,
    mark_as_read,
)
from amprenta_rag.notifications.email_service import (
    send_email,
    is_email_configured,
    send_experiment_summary,
    send_share_notification,
)

__all__ = [
    "create_notification",
    "get_unread_notifications",
    "get_unread_count",
    "mark_as_read",
    "mark_all_read",
    "send_email",
    "is_email_configured",
    "send_experiment_summary",
    "send_share_notification",
]
