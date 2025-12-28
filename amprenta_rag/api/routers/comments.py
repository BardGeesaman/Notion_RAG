"""Comments API endpoints."""
from __future__ import annotations

from typing import List
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db
from amprenta_rag.api import schemas
from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.utils.comments import (
    add_comment,
    get_comments,
    update_comment,
    delete_comment,
    parse_mentions,
    resolve_mentions,
    notify_mentions,
)
from amprenta_rag.services.activity import log_activity
from amprenta_rag.database.models import ActivityEventType, User

router = APIRouter(prefix="/comments", tags=["comments"])


@router.post("", response_model=schemas.CommentResponse)
def create_comment(
    comment: schemas.CommentCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> schemas.CommentResponse:
    """Create a new comment and notify mentioned users.
    
    Args:
        comment: Comment creation data
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        Created comment with metadata
    """
    user_id = current_user.id
    
    try:
        # Create comment
        new_comment = add_comment(
            entity_type=comment.entity_type,
            entity_id=comment.entity_id,
            content=comment.content,
            user_id=user_id,
            db=db,
            parent_id=comment.parent_id,
        )
        
        # Parse and notify mentions
        usernames = parse_mentions(comment.content)
        if usernames:
            mentioned_user_ids = resolve_mentions(usernames, db)
            notify_mentions(new_comment, mentioned_user_ids, db)
        
        # Log activity
        try:
            log_activity(
                event_type=ActivityEventType.COMMENT_ADDED,
                target_type=comment.entity_type,
                target_id=comment.entity_id,
                target_name=f"{comment.entity_type} comment",
                actor_id=user_id,
                program_id=None,
                metadata={"comment_id": str(new_comment.id)},
            )
        except Exception as e:
            # Don't fail comment creation if activity logging fails
            import logging
            logging.getLogger(__name__).error(f"Failed to log comment activity: {e}")
        
        # Get author name from User table
        user = db.query(User).filter(User.id == user_id).first()
        author_name = user.username if user else "Unknown"
        
        # Build response
        return schemas.CommentResponse(
            id=new_comment.id,
            entity_type=new_comment.entity_type,
            entity_id=new_comment.entity_id,
            content=new_comment.content,
            author=author_name,
            author_id=user_id,
            created_at=new_comment.created_at,
            updated_at=new_comment.updated_at,
            replies=[],
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to create comment: {e}")


@router.get("", response_model=List[dict])
def list_comments(
    entity_type: str,
    entity_id: UUID,
    db: Session = Depends(get_db),
) -> List[dict]:
    """Get comments for an entity.
    
    Args:
        entity_type: Type of entity (e.g., 'dataset', 'compound')
        entity_id: UUID of the entity
        db: Database session
        
    Returns:
        List of comments with replies
    """
    try:
        return get_comments(entity_type, entity_id, db)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch comments: {e}")


@router.put("/{comment_id}", response_model=schemas.CommentResponse)
def edit_comment(
    comment_id: UUID,
    comment_update: schemas.CommentUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> schemas.CommentResponse:
    """Update a comment (author only).
    
    Args:
        comment_id: Comment UUID
        comment_update: Updated content
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        Updated comment
        
    Raises:
        HTTPException: 403 if unauthorized, 404 if not found
    """
    user_id = current_user.id
    
    updated = update_comment(comment_id, comment_update.content, user_id, db)
    
    if updated is None:
        raise HTTPException(
            status_code=403,
            detail="Comment not found or you are not the author"
        )
    
    # Get author name from User table
    user = db.query(User).filter(User.id == user_id).first()
    author_name = user.username if user else "Unknown"
    
    return schemas.CommentResponse(
        id=updated.id,
        entity_type=updated.entity_type,
        entity_id=updated.entity_id,
        content=updated.content,
        author=author_name,
        author_id=user_id,
        created_at=updated.created_at,
        updated_at=updated.updated_at,
        replies=[],
    )


@router.delete("/{comment_id}")
def remove_comment(
    comment_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> dict:
    """Delete a comment (author only).
    
    Args:
        comment_id: Comment UUID
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        Success message
        
    Raises:
        HTTPException: 403 if unauthorized, 404 if not found
    """
    user_id = current_user.id
    
    deleted = delete_comment(comment_id, user_id, db)
    
    if not deleted:
        raise HTTPException(
            status_code=403,
            detail="Comment not found or you are not the author"
        )
    
    return {"status": "deleted", "comment_id": str(comment_id)}


__all__ = ["router"]

