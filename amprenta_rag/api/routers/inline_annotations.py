"""Inline annotations API endpoints."""

from __future__ import annotations

from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, status
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.api.schemas import (
    InlineAnnotationCreate,
    InlineAnnotationListResponse,
    InlineAnnotationReplyCreate,
    InlineAnnotationResponse,
)
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import User
from amprenta_rag.services.inline_annotations import (
    create_annotation,
    delete_annotation,
    get_annotation,
    get_annotation_count,
    get_annotations,
    reopen_annotation,
    reply_to_annotation,
    resolve_annotation,
)

router = APIRouter(prefix="/annotations", tags=["inline-annotations"])


@router.post("", response_model=InlineAnnotationResponse, status_code=status.HTTP_201_CREATED)
def create_inline_annotation(
    annotation: InlineAnnotationCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> InlineAnnotationResponse:
    """Create a new inline annotation.
    
    Args:
        annotation: Annotation creation data
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        Created annotation
        
    Raises:
        HTTPException: 400 if position data is invalid
    """
    try:
        new_annotation = create_annotation(
            entity_type=annotation.entity_type,
            entity_id=annotation.entity_id,
            position_type=annotation.position_type,
            position_data=annotation.position_data,
            content=annotation.content,
            user_id=current_user.id,
            db=db,
            parent_id=annotation.parent_id,
        )
        
        return InlineAnnotationResponse.model_validate(new_annotation)
        
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to create annotation: {e}")


@router.get("", response_model=InlineAnnotationListResponse)
def list_inline_annotations(
    entity_type: str = Query(..., description="Type of entity"),
    entity_id: UUID = Query(..., description="ID of the entity"),
    status_filter: Optional[str] = Query(None, alias="status", description="Filter by status (open, resolved)"),
    position_type: Optional[str] = Query(None, description="Filter by position type"),
    db: Session = Depends(get_db),
) -> InlineAnnotationListResponse:
    """Get inline annotations for an entity with optional filters.
    
    Args:
        entity_type: Type of entity (notebook, dataset, experiment)
        entity_id: ID of the entity
        status_filter: Optional status filter
        position_type: Optional position type filter
        db: Database session
        
    Returns:
        List of annotations with counts
    """
    try:
        # Get annotations with filters
        annotations = get_annotations(
            entity_type=entity_type,
            entity_id=entity_id,
            db=db,
            status=status_filter,
            position_type=position_type,
        )
        
        # Get counts for the entity
        counts = get_annotation_count(entity_type, entity_id, db)
        
        return InlineAnnotationListResponse(
            annotations=[InlineAnnotationResponse.model_validate(ann) for ann in annotations],
            total=counts["total"],
            open_count=counts["open"],
            resolved_count=counts["resolved"],
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch annotations: {e}")


@router.get("/{annotation_id}", response_model=InlineAnnotationResponse)
def get_inline_annotation(
    annotation_id: UUID,
    db: Session = Depends(get_db),
) -> InlineAnnotationResponse:
    """Get a single inline annotation by ID.
    
    Args:
        annotation_id: Annotation UUID
        db: Database session
        
    Returns:
        Annotation details
        
    Raises:
        HTTPException: 404 if annotation not found
    """
    annotation = get_annotation(annotation_id, db)
    
    if annotation is None:
        raise HTTPException(status_code=404, detail="Annotation not found")
        
    return InlineAnnotationResponse.model_validate(annotation)


@router.patch("/{annotation_id}/resolve", response_model=InlineAnnotationResponse)
def resolve_inline_annotation(
    annotation_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> InlineAnnotationResponse:
    """Mark an annotation as resolved.
    
    Args:
        annotation_id: Annotation UUID
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        Updated annotation
        
    Raises:
        HTTPException: 404 if annotation not found
    """
    resolved_annotation = resolve_annotation(annotation_id, current_user.id, db)
    
    if resolved_annotation is None:
        raise HTTPException(status_code=404, detail="Annotation not found")
        
    return InlineAnnotationResponse.model_validate(resolved_annotation)


@router.patch("/{annotation_id}/reopen", response_model=InlineAnnotationResponse)
def reopen_inline_annotation(
    annotation_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> InlineAnnotationResponse:
    """Reopen a resolved annotation.
    
    Args:
        annotation_id: Annotation UUID
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        Updated annotation
        
    Raises:
        HTTPException: 404 if annotation not found
    """
    reopened_annotation = reopen_annotation(annotation_id, current_user.id, db)
    
    if reopened_annotation is None:
        raise HTTPException(status_code=404, detail="Annotation not found")
        
    return InlineAnnotationResponse.model_validate(reopened_annotation)


@router.post("/{annotation_id}/reply", response_model=InlineAnnotationResponse, status_code=status.HTTP_201_CREATED)
def reply_to_inline_annotation(
    annotation_id: UUID,
    reply: InlineAnnotationReplyCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> InlineAnnotationResponse:
    """Add a reply to an annotation.
    
    Args:
        annotation_id: Parent annotation UUID
        reply: Reply content
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        Created reply annotation
        
    Raises:
        HTTPException: 404 if parent annotation not found
    """
    reply_annotation = reply_to_annotation(
        annotation_id=annotation_id,
        content=reply.content,
        user_id=current_user.id,
        db=db,
    )
    
    if reply_annotation is None:
        raise HTTPException(status_code=404, detail="Parent annotation not found")
        
    return InlineAnnotationResponse.model_validate(reply_annotation)


@router.delete("/{annotation_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_inline_annotation(
    annotation_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> None:
    """Delete an annotation (author only).
    
    Args:
        annotation_id: Annotation UUID
        db: Database session
        current_user: Current authenticated user
        
    Raises:
        HTTPException: 403 if unauthorized, 404 if not found
    """
    deleted = delete_annotation(annotation_id, current_user.id, db)
    
    if not deleted:
        raise HTTPException(
            status_code=403,
            detail="Annotation not found or you are not authorized to delete it"
        )


__all__ = ["router"]
