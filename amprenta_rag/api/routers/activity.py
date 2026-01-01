"""
Activity and notifications API endpoints.

This module provides REST API endpoints for:
- Activity feed management
- Notification management
- Activity event retrieval
"""

from datetime import datetime
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.api.schemas import (
    ActivityEventSchema,
    NotificationSchema,
    NotificationCountSchema,
)
from amprenta_rag.services.activity import (
    get_activity_feed,
    get_user_notifications,
    mark_notification_read,
    mark_all_notifications_read,
    get_unread_count,
)
from amprenta_rag.database.models import ActivityEvent, User
from amprenta_rag.database.session import db_session

# Create routers
activity_router = APIRouter(prefix="/activity", tags=["activity"])
notifications_router = APIRouter(prefix="/notifications", tags=["notifications"])

# Mock user ID for MVP (replace with proper auth later)
MOCK_USER_ID = UUID("00000000-0000-0000-0000-000000000001")


@activity_router.get("/feed", response_model=List[ActivityEventSchema])
def get_activity_feed_endpoint(
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    event_type: Optional[str] = Query(None, description="Filter by event type"),
    since: Optional[datetime] = Query(None, description="Filter events after this datetime"),
    limit: int = Query(50, ge=1, le=100, description="Maximum number of events to return"),
    offset: int = Query(0, ge=0, description="Number of events to skip"),
) -> List[ActivityEventSchema]:
    """
    Get paginated activity feed with optional filters.
    
    Returns a list of activity events ordered by creation date (newest first).
    """
    events = get_activity_feed(
        program_id=program_id,
        event_type=event_type,
        since=since,
        limit=limit,
        offset=offset,
    )
    
    return [
        ActivityEventSchema(
            id=event.id,
            event_type=event.event_type,
            actor_id=event.actor_id,
            target_type=event.target_type,
            target_id=event.target_id,
            target_name=event.target_name,
            program_id=event.program_id,
            metadata=event.event_metadata,
            created_at=event.created_at,
        )
        for event in events
    ]


@activity_router.get("/events/{event_id}", response_model=ActivityEventSchema)
def get_activity_event(event_id: UUID) -> ActivityEventSchema:
    """
    Get a specific activity event by ID.
    
    Returns 404 if the event is not found.
    """
    try:
        with db_session() as db:
            event = db.query(ActivityEvent).filter(ActivityEvent.id == event_id).first()
            
            if not event:
                raise HTTPException(status_code=404, detail="Activity event not found")
            
            return ActivityEventSchema(
                id=event.id,
                event_type=event.event_type,
                actor_id=event.actor_id,
                target_type=event.target_type,
                target_id=event.target_id,
                target_name=event.target_name,
                program_id=event.program_id,
                metadata=event.event_metadata,
                created_at=event.created_at,
            )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to retrieve activity event: {e}")


@notifications_router.get("", response_model=List[NotificationSchema])
def get_notifications_endpoint(
    unread_only: bool = Query(False, description="Return only unread notifications"),
    limit: int = Query(20, ge=1, le=100, description="Maximum number of notifications to return"),
    offset: int = Query(0, ge=0, description="Number of notifications to skip"),
    current_user: User = Depends(get_current_user),
) -> List[NotificationSchema]:
    """
    Get paginated notifications for the current user.
    
    Returns a list of notifications ordered by creation date (newest first).
    Includes related activity event data when available.
    """
    user_id = current_user.id
    
    notifications = get_user_notifications(
        user_id=user_id,
        unread_only=unread_only,
        limit=limit,
        offset=offset,
    )
    
    result = []
    for notification in notifications:
        # Build activity event schema if available
        activity_event_data = None
        if notification.activity_event:
            activity_event_data = ActivityEventSchema(
                id=notification.activity_event.id,
                event_type=notification.activity_event.event_type,
                actor_id=notification.activity_event.actor_id,
                target_type=notification.activity_event.target_type,
                target_id=notification.activity_event.target_id,
                target_name=notification.activity_event.target_name,
                program_id=notification.activity_event.program_id,
                metadata=notification.activity_event.event_metadata,
                created_at=notification.activity_event.created_at,
            )
        
        result.append(
            NotificationSchema(
                id=notification.id,
                user_id=user_id,  # Use the mock user ID
                activity_event_id=notification.activity_event_id,
                notification_type=notification.notification_type,
                title=f"Activity: {notification.activity_event.event_type}" if notification.activity_event else None,
                message=f"{notification.activity_event.target_name}" if notification.activity_event else None,
                is_read=notification.is_read,
                created_at=notification.created_at,
                activity_event=activity_event_data,
            )
        )
    
    return result


@notifications_router.post("/{notification_id}/read")
def mark_notification_read_endpoint(
    notification_id: UUID,
    current_user: User = Depends(get_current_user),
) -> dict:
    """
    Mark a specific notification as read.
    
    Returns 404 if the notification is not found or not owned by the current user.
    """
    user_id = current_user.id
    
    success = mark_notification_read(notification_id, user_id)
    
    if not success:
        raise HTTPException(
            status_code=404, 
            detail="Notification not found or not accessible"
        )
    
    return {
        "status": "read",
        "notification_id": str(notification_id)
    }


@notifications_router.post("/read-all")
def mark_all_notifications_read_endpoint(
    current_user: User = Depends(get_current_user),
) -> dict:
    """
    Mark all notifications as read for the current user.
    
    Returns the count of notifications that were marked as read.
    """
    user_id = current_user.id
    
    count = mark_all_notifications_read(user_id)
    
    return {
        "status": "read",
        "count": count
    }


@notifications_router.get("/count", response_model=NotificationCountSchema)
def get_notifications_count(
    current_user: User = Depends(get_current_user),
) -> NotificationCountSchema:
    """
    Get the count of unread notifications for the current user.
    """
    user_id = current_user.id
    
    unread_count = get_unread_count(user_id)
    
    return NotificationCountSchema(unread_count=unread_count)


# Export both routers for registration in main.py
__all__ = ["activity_router", "notifications_router"]
