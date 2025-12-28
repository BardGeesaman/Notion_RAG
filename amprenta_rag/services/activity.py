"""
Activity service for logging events and managing notifications.

This service handles:
- Logging activity events (compound added, model trained, etc.)
- Routing notifications to relevant users
- Managing activity feeds and notification states
"""

from datetime import datetime
from typing import Dict, List, Optional
from uuid import UUID

from amprenta_rag.database.models import (
    ActivityEvent,
    ActivityEventType,
    RepositoryNotification,
    Program,
    User,
)
from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def log_activity(
    event_type: ActivityEventType,
    target_type: str,
    target_id: UUID,
    target_name: str,
    actor_id: Optional[UUID] = None,
    program_id: Optional[UUID] = None,
    metadata: Optional[Dict] = None,
) -> Optional[ActivityEvent]:
    """
    Log an activity event and route notifications.
    
    Args:
        event_type: Type of activity event
        target_type: Type of target object (e.g., "compound", "experiment")
        target_id: UUID of target object
        target_name: Human-readable name of target (for orphan handling)
        actor_id: UUID of user who performed the action
        program_id: UUID of associated program
        metadata: Additional event metadata
        
    Returns:
        Created ActivityEvent or None on failure
    """
    try:
        with db_session() as db:
            # Handle both enum and string event types
            if isinstance(event_type, ActivityEventType):
                event_type_str = event_type.value
            else:
                event_type_str = str(event_type)
            
            # Create activity event
            event = ActivityEvent(
                event_type=event_type_str,
                actor_id=actor_id,
                target_type=target_type,
                target_id=target_id,
                target_name=target_name,
                program_id=program_id,
                event_metadata=metadata or {},
            )
            
            db.add(event)
            db.flush()  # Get the event ID
            
            # Route notifications
            notified_users = route_notifications(event, db)
            logger.info(
                f"Activity logged: {event_type_str} for {target_type}:{target_id}, "
                f"notified {len(notified_users)} users"
            )
            
            db.commit()
            
            # Load all attributes before returning (to avoid DetachedInstanceError)
            event_id = event.id
            event_type_val = event.event_type
            actor_id_val = event.actor_id
            target_type_val = event.target_type
            target_id_val = event.target_id
            target_name_val = event.target_name
            program_id_val = event.program_id
            metadata_val = event.event_metadata
            created_at_val = event.created_at
            
            db.expunge(event)  # Detach from session so it can be used after context
            return event
            
    except Exception as e:
        # Use event_type_str if available, otherwise convert event_type
        event_type_display = event_type_str if 'event_type_str' in locals() else (
            event_type.value if isinstance(event_type, ActivityEventType) else str(event_type)
        )
        logger.error(f"Failed to log activity {event_type_display}: {e}")
        return None


def route_notifications(event: ActivityEvent, db) -> List[UUID]:
    """
    Route notifications for an activity event.
    
    MVP implementation: Notify program lead only.
    
    Args:
        event: ActivityEvent to route
        db: Database session
        
    Returns:
        List of user IDs that were notified
    """
    notified_users = []
    
    try:
        if not event.program_id:
            logger.debug(f"No program_id for event {event.id}, skipping notifications")
            return notified_users
            
        # Get program creator (as lead)
        program = db.query(Program).filter(Program.id == event.program_id).first()
        if not program or not program.created_by_id:
            logger.debug(f"No program creator found for program {event.program_id}")
            return notified_users
            
        # Don't notify the actor about their own actions
        if event.actor_id == program.created_by_id:
            logger.debug(f"Actor {event.actor_id} is program creator, skipping self-notification")
            return notified_users
            
        # Create notification for program lead
        notification = RepositoryNotification(
            subscription_id=program.id,  # Use program_id as subscription_id for now
            dataset_id=event.target_id,  # Use target_id as dataset_id for now
            activity_event_id=event.id,
            notification_type="activity",
            is_read=False,
        )
        
        db.add(notification)
        notified_users.append(program.created_by_id)
        
        logger.debug(f"Created notification for program creator {program.created_by_id}")
        
    except Exception as e:
        logger.error(f"Failed to route notifications for event {event.id}: {e}")
        
    return notified_users


def get_activity_feed(
    program_id: Optional[UUID] = None,
    event_type: Optional[str] = None,
    since: Optional[datetime] = None,
    limit: int = 50,
    offset: int = 0,
) -> List[ActivityEvent]:
    """
    Get paginated activity feed with optional filters.
    
    Args:
        program_id: Filter by program
        event_type: Filter by event type
        since: Filter events after this datetime
        limit: Maximum number of events to return
        offset: Number of events to skip
        
    Returns:
        List of ActivityEvent objects
    """
    try:
        with db_session() as db:
            query = db.query(ActivityEvent)
            
            # Apply filters
            if program_id:
                query = query.filter(ActivityEvent.program_id == program_id)
            
            if event_type:
                query = query.filter(ActivityEvent.event_type == event_type)
                
            if since:
                query = query.filter(ActivityEvent.created_at >= since)
            
            # Order and paginate
            events = (
                query.order_by(ActivityEvent.created_at.desc())
                .offset(offset)
                .limit(limit)
                .all()
            )
            
            # Load all attributes and expunge to avoid DetachedInstanceError
            for event in events:
                # Access all attributes to load them
                _ = (event.id, event.event_type, event.actor_id, event.target_type,
                     event.target_id, event.target_name, event.program_id,
                     event.event_metadata, event.created_at)
                db.expunge(event)
            
            return events
            
    except Exception as e:
        logger.error(f"Failed to get activity feed: {e}")
        return []


def get_user_notifications(
    user_id: UUID,
    unread_only: bool = False,
    limit: int = 20,
    offset: int = 0,
) -> List[RepositoryNotification]:
    """
    Get paginated notifications for a user.
    
    Args:
        user_id: UUID of user
        unread_only: If True, only return unread notifications
        limit: Maximum number of notifications to return
        offset: Number of notifications to skip
        
    Returns:
        List of RepositoryNotification objects with activity_event loaded
    """
    try:
        with db_session() as db:
            # Note: This is a simplified implementation
            # In a real system, we'd need proper user-notification relationships
            query = db.query(RepositoryNotification).join(
                ActivityEvent, 
                RepositoryNotification.activity_event_id == ActivityEvent.id
            ).join(
                Program,
                ActivityEvent.program_id == Program.id
            ).filter(
                Program.created_by_id == user_id
            )
            
            if unread_only:
                query = query.filter(RepositoryNotification.is_read == False)
            
            notifications = (
                query.order_by(RepositoryNotification.created_at.desc())
                .offset(offset)
                .limit(limit)
                .all()
            )
            
            return notifications
            
    except Exception as e:
        logger.error(f"Failed to get notifications for user {user_id}: {e}")
        return []


def mark_notification_read(notification_id: UUID, user_id: UUID) -> bool:
    """
    Mark a notification as read.
    
    Args:
        notification_id: UUID of notification
        user_id: UUID of user (for security check)
        
    Returns:
        True if successful, False otherwise
    """
    try:
        with db_session() as db:
            # Find notification that belongs to this user
            notification = (
                db.query(RepositoryNotification)
                .join(ActivityEvent, RepositoryNotification.activity_event_id == ActivityEvent.id)
                .join(Program, ActivityEvent.program_id == Program.id)
                .filter(
                    RepositoryNotification.id == notification_id,
                    Program.created_by_id == user_id
                )
                .first()
            )
            
            if not notification:
                logger.warning(f"Notification {notification_id} not found for user {user_id}")
                return False
                
            notification.is_read = True
            db.commit()
            
            logger.debug(f"Marked notification {notification_id} as read for user {user_id}")
            return True
            
    except Exception as e:
        logger.error(f"Failed to mark notification {notification_id} as read: {e}")
        return False


def mark_all_notifications_read(user_id: UUID) -> int:
    """
    Mark all notifications as read for a user.
    
    Args:
        user_id: UUID of user
        
    Returns:
        Number of notifications marked as read
    """
    try:
        with db_session() as db:
            # Update all unread notifications for this user
            count = (
                db.query(RepositoryNotification)
                .join(ActivityEvent, RepositoryNotification.activity_event_id == ActivityEvent.id)
                .join(Program, ActivityEvent.program_id == Program.id)
                .filter(
                    Program.created_by_id == user_id,
                    RepositoryNotification.is_read == False
                )
                .update({"is_read": True})
            )
            
            db.commit()
            
            logger.info(f"Marked {count} notifications as read for user {user_id}")
            return count
            
    except Exception as e:
        logger.error(f"Failed to mark all notifications as read for user {user_id}: {e}")
        return 0


def get_unread_count(user_id: UUID) -> int:
    """
    Get count of unread notifications for a user.
    
    Args:
        user_id: UUID of user
        
    Returns:
        Number of unread notifications
    """
    try:
        with db_session() as db:
            count = (
                db.query(RepositoryNotification)
                .join(ActivityEvent, RepositoryNotification.activity_event_id == ActivityEvent.id)
                .join(Program, ActivityEvent.program_id == Program.id)
                .filter(
                    Program.created_by_id == user_id,
                    RepositoryNotification.is_read == False
                )
                .count()
            )
            
            return count
            
    except Exception as e:
        logger.error(f"Failed to get unread count for user {user_id}: {e}")
        return 0
