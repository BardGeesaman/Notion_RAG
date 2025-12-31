"""Service layer for review thread and comment operations."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session, selectinload

from amprenta_rag.database.models import ReviewThread, ReviewComment
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_thread(
    db: Session,
    review_id: UUID,
    title: str,
    created_by_id: UUID,
    cell_index: Optional[int] = None
) -> ReviewThread:
    """
    Create a new review thread.
    
    Args:
        db: Database session
        review_id: ID of the notebook review
        title: Thread title/summary
        created_by_id: ID of user creating the thread
        cell_index: Optional cell index for cell-specific threads
    
    Returns:
        Created ReviewThread instance
    """
    logger.info(f"Creating thread for review {review_id}: {title}")
    
    thread = ReviewThread(
        review_id=review_id,
        title=title,
        created_by_id=created_by_id,
        cell_index=cell_index,
        status="open"
    )
    
    db.add(thread)
    db.commit()
    db.refresh(thread)
    
    logger.info(f"Created thread {thread.id} for review {review_id}")
    return thread


def get_threads(db: Session, review_id: UUID) -> List[ReviewThread]:
    """
    Get all threads for a review with comments eager loaded.
    
    Args:
        db: Database session
        review_id: ID of the notebook review
    
    Returns:
        List of ReviewThread instances with comments loaded
    """
    logger.debug(f"Fetching threads for review {review_id}")
    
    threads = db.query(ReviewThread).options(
        selectinload(ReviewThread.comments).selectinload(ReviewComment.created_by),
        selectinload(ReviewThread.created_by)
    ).filter(
        ReviewThread.review_id == review_id
    ).order_by(
        ReviewThread.created_at.asc()
    ).all()
    
    logger.debug(f"Found {len(threads)} threads for review {review_id}")
    return threads


def add_comment(
    db: Session,
    thread_id: UUID,
    content: str,
    created_by_id: UUID,
    parent_id: Optional[UUID] = None
) -> ReviewComment:
    """
    Add a comment to a thread.
    
    Args:
        db: Database session
        thread_id: ID of the thread
        content: Comment content (markdown)
        created_by_id: ID of user creating the comment
        parent_id: Optional parent comment ID for nested replies
    
    Returns:
        Created ReviewComment instance
    """
    logger.info(f"Adding comment to thread {thread_id}")
    
    comment = ReviewComment(
        thread_id=thread_id,
        content=content,
        created_by_id=created_by_id,
        parent_id=parent_id
    )
    
    db.add(comment)
    db.commit()
    db.refresh(comment)
    
    logger.info(f"Created comment {comment.id} in thread {thread_id}")
    return comment


def resolve_thread(db: Session, thread_id: UUID, status: str = "resolved") -> ReviewThread:
    """
    Mark a thread as resolved or set custom status.
    
    Args:
        db: Database session
        thread_id: ID of the thread
        status: New status (resolved, wontfix, etc.)
    
    Returns:
        Updated ReviewThread instance
    
    Raises:
        ValueError: If thread not found
    """
    logger.info(f"Resolving thread {thread_id} with status: {status}")
    
    thread = db.query(ReviewThread).filter(ReviewThread.id == thread_id).first()
    if not thread:
        raise ValueError(f"Thread {thread_id} not found")
    
    thread.status = status
    db.commit()
    db.refresh(thread)
    
    logger.info(f"Thread {thread_id} marked as {status}")
    return thread


def reopen_thread(db: Session, thread_id: UUID) -> ReviewThread:
    """
    Reopen a resolved thread.
    
    Args:
        db: Database session
        thread_id: ID of the thread
    
    Returns:
        Updated ReviewThread instance
    
    Raises:
        ValueError: If thread not found
    """
    logger.info(f"Reopening thread {thread_id}")
    
    thread = db.query(ReviewThread).filter(ReviewThread.id == thread_id).first()
    if not thread:
        raise ValueError(f"Thread {thread_id} not found")
    
    thread.status = "open"
    db.commit()
    db.refresh(thread)
    
    logger.info(f"Thread {thread_id} reopened")
    return thread


def delete_thread(db: Session, thread_id: UUID) -> bool:
    """
    Delete a thread and all its comments.
    
    Args:
        db: Database session
        thread_id: ID of the thread to delete
    
    Returns:
        True if thread was deleted, False if not found
    """
    logger.info(f"Deleting thread {thread_id}")
    
    thread = db.query(ReviewThread).filter(ReviewThread.id == thread_id).first()
    if not thread:
        logger.warning(f"Thread {thread_id} not found for deletion")
        return False
    
    db.delete(thread)
    db.commit()
    
    logger.info(f"Thread {thread_id} deleted")
    return True
