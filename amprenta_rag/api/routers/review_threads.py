"""REST API endpoints for notebook review threads and diff functionality."""

from __future__ import annotations

from datetime import datetime
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, ConfigDict
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session, get_current_user
from amprenta_rag.database.models import User
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.services.review_threads import (
    create_thread,
    get_threads,
    add_comment,
    resolve_thread,
    reopen_thread,
)
from amprenta_rag.services.notebook_diff import get_review_diff

logger = get_logger(__name__)
router = APIRouter(prefix="/reviews", tags=["Review Threads"])


# Pydantic Schemas
class ThreadCreateRequest(BaseModel):
    """Request schema for creating a new review thread."""
    title: str
    cell_index: Optional[int] = None


class ThreadUpdateRequest(BaseModel):
    """Request schema for updating thread status."""
    status: str  # "open", "resolved", "wontfix"


class CommentCreateRequest(BaseModel):
    """Request schema for creating a comment."""
    content: str
    parent_id: Optional[UUID] = None


class CommentResponse(BaseModel):
    """Response schema for a comment."""
    id: UUID
    thread_id: UUID
    parent_id: Optional[UUID]
    content: str
    created_by_id: Optional[UUID]
    created_at: datetime

    model_config = ConfigDict(from_attributes=True)


class ThreadResponse(BaseModel):
    """Response schema for a review thread."""
    id: UUID
    review_id: UUID
    title: str
    cell_index: Optional[int]
    status: str
    created_by_id: Optional[UUID]
    created_at: datetime
    comments: List[CommentResponse]

    model_config = ConfigDict(from_attributes=True)


class CellSummary(BaseModel):
    """Summary of a notebook cell for diff display."""
    index: int
    cell_type: str
    source_preview: str


class CellModification(BaseModel):
    """Summary of a modified cell showing before/after."""
    index: int
    cell_type: str
    old_source_preview: str
    new_source_preview: str


class DiffResponse(BaseModel):
    """Response schema for notebook diff."""
    added: List[CellSummary]
    removed: List[CellSummary]
    modified: List[CellModification]
    unchanged: List[CellSummary]


# API Endpoints
@router.post("/{review_id}/threads", response_model=ThreadResponse, status_code=status.HTTP_201_CREATED)
def create_review_thread(
    review_id: UUID,
    request: ThreadCreateRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Create a new discussion thread for a notebook review."""
    logger.info(f"User {current_user.id} creating thread for review {review_id}")
    
    try:
        thread = create_thread(
            db=db,
            review_id=review_id,
            title=request.title,
            created_by_id=current_user.id,
            cell_index=request.cell_index,
        )
        
        # Convert to response format with empty comments list initially
        return ThreadResponse(
            id=thread.id,
            review_id=thread.review_id,
            title=thread.title,
            cell_index=thread.cell_index,
            status=thread.status,
            created_by_id=thread.created_by_id,
            created_at=thread.created_at,
            comments=[],
        )
    except ValueError as e:
        logger.warning(f"Invalid review ID {review_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Review {review_id} not found",
        )
    except Exception as e:
        logger.error(f"Failed to create thread for review {review_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to create thread",
        )


@router.get("/{review_id}/threads", response_model=List[ThreadResponse])
def list_review_threads(
    review_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """List all discussion threads for a notebook review."""
    logger.info(f"User {current_user.id} listing threads for review {review_id}")
    
    try:
        threads = get_threads(db=db, review_id=review_id)
        
        # Convert threads with their comments to response format
        thread_responses = []
        for thread in threads:
            comments = [
                CommentResponse(
                    id=comment.id,
                    thread_id=comment.thread_id,
                    parent_id=comment.parent_id,
                    content=comment.content,
                    created_by_id=comment.created_by_id,
                    created_at=comment.created_at,
                )
                for comment in thread.comments
            ]
            
            thread_responses.append(ThreadResponse(
                id=thread.id,
                review_id=thread.review_id,
                title=thread.title,
                cell_index=thread.cell_index,
                status=thread.status,
                created_by_id=thread.created_by_id,
                created_at=thread.created_at,
                comments=comments,
            ))
        
        return thread_responses
    except Exception as e:
        logger.error(f"Failed to list threads for review {review_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to list threads",
        )


# Thread-specific endpoints (not nested under reviews)
threads_router = APIRouter(prefix="/threads", tags=["Review Threads"])


@threads_router.post("/{thread_id}/comments", response_model=CommentResponse, status_code=status.HTTP_201_CREATED)
def add_thread_comment(
    thread_id: UUID,
    request: CommentCreateRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Add a comment to a discussion thread."""
    logger.info(f"User {current_user.id} adding comment to thread {thread_id}")
    
    try:
        comment = add_comment(
            db=db,
            thread_id=thread_id,
            content=request.content,
            created_by_id=current_user.id,
            parent_id=request.parent_id,
        )
        
        return CommentResponse(
            id=comment.id,
            thread_id=comment.thread_id,
            parent_id=comment.parent_id,
            content=comment.content,
            created_by_id=comment.created_by_id,
            created_at=comment.created_at,
        )
    except ValueError as e:
        logger.warning(f"Invalid thread ID {thread_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Thread {thread_id} not found",
        )
    except Exception as e:
        logger.error(f"Failed to add comment to thread {thread_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to add comment",
        )


@threads_router.patch("/{thread_id}", response_model=ThreadResponse)
def update_thread_status(
    thread_id: UUID,
    request: ThreadUpdateRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Update the status of a discussion thread (resolve/reopen)."""
    logger.info(f"User {current_user.id} updating thread {thread_id} status to {request.status}")
    
    if request.status not in ["open", "resolved", "wontfix"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Status must be 'open', 'resolved', or 'wontfix'",
        )
    
    try:
        if request.status == "open":
            thread = reopen_thread(db=db, thread_id=thread_id)
        else:
            thread = resolve_thread(db=db, thread_id=thread_id, status=request.status)
        
        return ThreadResponse(
            id=thread.id,
            review_id=thread.review_id,
            title=thread.title,
            cell_index=thread.cell_index,
            status=thread.status,
            created_by_id=thread.created_by_id,
            created_at=thread.created_at,
            comments=[],  # Don't fetch comments for status update
        )
    except ValueError as e:
        logger.warning(f"Invalid thread ID {thread_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Thread {thread_id} not found",
        )
    except Exception as e:
        logger.error(f"Failed to update thread {thread_id} status: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to update thread status",
        )


@router.get("/{review_id}/diff", response_model=DiffResponse)
def get_notebook_diff(
    review_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Get the notebook diff between snapshot and current version."""
    logger.info(f"User {current_user.id} requesting diff for review {review_id}")
    
    try:
        # We need to get the notebook path from the review to compute the diff
        from amprenta_rag.database.models import NotebookReview
        
        review = db.query(NotebookReview).filter(NotebookReview.id == review_id).first()
        if not review:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Review {review_id} not found",
            )
        
        diff_result = get_review_diff(
            db=db,
            review_id=review_id,
            current_notebook_path=review.notebook_path,
        )
        
        # Convert diff result to response format
        def convert_cell_summary(cell_data):
            return CellSummary(
                index=cell_data["index"],
                cell_type=cell_data["cell_type"],
                source_preview=cell_data["source_preview"],
            )
        
        def convert_cell_modification(cell_data):
            return CellModification(
                index=cell_data["index"],
                cell_type=cell_data["cell_type"],
                old_source_preview=cell_data["old_source_preview"],
                new_source_preview=cell_data["new_source_preview"],
            )
        
        return DiffResponse(
            added=[convert_cell_summary(cell) for cell in diff_result["added"]],
            removed=[convert_cell_summary(cell) for cell in diff_result["removed"]],
            modified=[convert_cell_modification(cell) for cell in diff_result["modified"]],
            unchanged=[convert_cell_summary(cell) for cell in diff_result["unchanged"]],
        )
    except ValueError as e:
        logger.warning(f"Invalid review or diff computation failed for {review_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e),
        )
    except Exception as e:
        logger.error(f"Failed to get diff for review {review_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to compute notebook diff",
        )


# Export both routers for registration
review_threads_router = router
thread_actions_router = threads_router
