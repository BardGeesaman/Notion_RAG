"""Comment service for contextual commenting on entities."""

import re
from datetime import datetime, timezone
from typing import Any, List, Optional, cast
from uuid import UUID
from sqlalchemy.orm import Session
from amprenta_rag.database.models import Comment, User
from amprenta_rag.utils.uuid_utils import ensure_uuid
from amprenta_rag.services.activity import log_activity


def add_comment(
    entity_type: str,
    entity_id: UUID,
    content: str,
    user_id: UUID,
    db: Session,
    parent_id: Optional[UUID] = None,
) -> Comment:
    """Add a comment to an entity."""
    comment = Comment(
        entity_type=entity_type,
        entity_id=cast(Any, ensure_uuid(entity_id)),
        content=content,
        created_by_id=cast(Any, ensure_uuid(user_id)),
        parent_id=cast(Any, ensure_uuid(parent_id)),
    )
    db.add(comment)
    db.commit()
    db.refresh(comment)
    return comment


def get_comments(entity_type: str, entity_id: UUID, db: Session) -> List[dict]:
    """Get threaded comments for an entity."""
    comments = (
        db.query(Comment)
        .filter(
            Comment.entity_type == entity_type,
            Comment.entity_id == ensure_uuid(entity_id),
            Comment.parent_id.is_(None),  # Top-level only
        )
        .order_by(Comment.created_at.desc())
        .all()
    )

    result = []
    for c in comments:
        user = db.query(User).filter(User.id == c.created_by_id).first()
        parent_uuid = cast(UUID, ensure_uuid(cast(UUID | None, c.id)))
        result.append(
            {
                "id": str(c.id),
                "content": c.content,
                "author": user.username if user else "Unknown",
                "author_id": str(c.created_by_id) if c.created_by_id else None,
                "created_at": c.created_at,
                "replies": _get_replies(parent_uuid, db),
            }
        )
    return result


def _get_replies(parent_id: UUID, db: Session) -> List[dict]:
    """Get replies to a comment."""
    replies = (
        db.query(Comment)
        .filter(Comment.parent_id == parent_id)
        .order_by(Comment.created_at)
        .all()
    )
    result = []
    for r in replies:
        user = db.query(User).filter(User.id == r.created_by_id).first()
        result.append(
            {
                "id": str(r.id),
                "content": r.content,
                "author": user.username if user else "Unknown",
                "author_id": str(r.created_by_id) if r.created_by_id else None,
                "created_at": r.created_at,
            }
        )
    return result


def delete_comment(comment_id: UUID, user_id: UUID, db: Session) -> bool:
    """Delete a comment (only if user is the author)."""
    comment = (
        db.query(Comment)
        .filter(Comment.id == comment_id, Comment.created_by_id == user_id)
        .first()
    )
    if comment:
        db.delete(comment)
        db.commit()
        return True
    return False


def update_comment(comment_id: UUID, content: str, user_id: UUID, db: Session) -> Optional[Comment]:
    """Update a comment (only if user is the author).
    
    Args:
        comment_id: Comment UUID
        content: New content
        user_id: User UUID attempting the update
        db: Database session
        
    Returns:
        Updated Comment or None if unauthorized/not found
    """
    comment = (
        db.query(Comment)
        .filter(Comment.id == comment_id, Comment.created_by_id == user_id)
        .first()
    )
    if comment:
        comment.content = content
        comment.updated_at = datetime.now(timezone.utc)
        db.commit()
        db.refresh(comment)
        return comment
    return None


def parse_mentions(content: str) -> List[str]:
    """Parse @mentions from comment content.
    
    Args:
        content: Comment text
        
    Returns:
        List of usernames (without @ prefix)
    """
    pattern = r'@([a-zA-Z0-9_-]+)'
    matches = re.findall(pattern, content)
    return matches


def resolve_mentions(usernames: List[str], db: Session) -> List[UUID]:
    """Resolve usernames to user IDs.
    
    Args:
        usernames: List of usernames
        db: Database session
        
    Returns:
        List of user UUIDs
    """
    if not usernames:
        return []
    
    users = db.query(User.id).filter(User.username.in_(usernames)).all()
    return [user.id for user in users]


def notify_mentions(comment: Comment, user_ids: List[UUID], db: Session) -> None:
    """Create activity events for mentioned users.
    
    Args:
        comment: The comment containing mentions
        user_ids: List of mentioned user IDs
        db: Database session
    """
    author_id = comment.created_by_id
    
    for user_id in user_ids:
        # Don't notify the comment author
        if user_id == author_id:
            continue
        
        try:
            log_activity(
                event_type="mention_received",
                target_type=comment.entity_type,
                target_id=comment.entity_id,
                target_name=f"Comment on {comment.entity_type}",
                actor_id=author_id,
                program_id=None,
                metadata={
                    "comment_id": str(comment.id),
                    "comment_preview": comment.content[:100],
                },
            )
        except Exception as e:
            # Log but don't fail the comment operation
            import logging
            logging.getLogger(__name__).error(f"Failed to notify mention: {e}")

