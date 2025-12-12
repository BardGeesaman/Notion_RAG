"""Comment service for contextual commenting on entities."""

from typing import List, Optional
from uuid import UUID
from sqlalchemy.orm import Session
from amprenta_rag.database.models import Comment, User


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
        entity_id=entity_id,
        content=content,
        created_by_id=user_id,
        parent_id=parent_id,
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
            Comment.entity_id == entity_id,
            Comment.parent_id == None,  # Top-level only
        )
        .order_by(Comment.created_at.desc())
        .all()
    )

    result = []
    for c in comments:
        user = db.query(User).filter(User.id == c.created_by_id).first()
        result.append(
            {
                "id": str(c.id),
                "content": c.content,
                "author": user.username if user else "Unknown",
                "author_id": str(c.created_by_id) if c.created_by_id else None,
                "created_at": c.created_at,
                "replies": _get_replies(c.id, db),
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

