"""Inline annotations service for position-anchored comments."""

from datetime import datetime, timezone
from typing import Any, Dict, List, Optional
from uuid import UUID

from sqlalchemy import and_, func
from sqlalchemy.orm import Session

from amprenta_rag.models.user_prefs import InlineAnnotation


# Position schema validation
POSITION_SCHEMAS = {
    "cell": {"required": ["cell_index"]},
    "cell_range": {"required": ["start_cell", "end_cell"]},
    "column": {"required": ["column"]},
    "row": {"required": ["row_index"]},
    "field": {"required": ["field"]},
    "range": {"required": ["start_cell", "end_cell"]},  # Alternative name for cell_range
}


def validate_position_data(position_type: str, position_data: Dict[str, Any]) -> bool:
    """Validate position_data matches expected schema for position_type.
    
    Args:
        position_type: The type of position (cell, column, row, field, etc.)
        position_data: Dictionary containing position-specific data
        
    Returns:
        True if valid, False otherwise
        
    Examples:
        >>> validate_position_data("cell", {"cell_index": 5})
        True
        >>> validate_position_data("column", {"column": "gene_name"})
        True
        >>> validate_position_data("cell", {"invalid_key": 5})
        False
    """
    if position_type not in POSITION_SCHEMAS:
        return False
        
    schema = POSITION_SCHEMAS[position_type]
    required_fields = schema["required"]
    
    # Check all required fields are present
    for field in required_fields:
        if field not in position_data:
            return False
            
    return True


def create_annotation(
    entity_type: str,
    entity_id: UUID,
    position_type: str,
    position_data: Dict[str, Any],
    content: str,
    user_id: Optional[UUID],
    db: Session,
    parent_id: Optional[UUID] = None,
) -> InlineAnnotation:
    """Create a new inline annotation.
    
    Args:
        entity_type: Type of entity being annotated (notebook, dataset, etc.)
        entity_id: ID of the entity being annotated
        position_type: Type of position (cell, column, row, field, etc.)
        position_data: Dictionary containing position-specific data
        content: Annotation content/text
        user_id: ID of user creating annotation (None for system annotations)
        db: Database session
        parent_id: Optional parent annotation ID for replies
        
    Returns:
        Created InlineAnnotation instance
        
    Raises:
        ValueError: If position_data is invalid for the position_type
    """
    # Validate position data
    if not validate_position_data(position_type, position_data):
        raise ValueError(f"Invalid position_data for position_type '{position_type}': {position_data}")
    
    annotation = InlineAnnotation(
        entity_type=entity_type,
        entity_id=entity_id,
        position_type=position_type,
        position_data=position_data,
        content=content,
        created_by_id=user_id,
        parent_id=parent_id,
    )
    
    db.add(annotation)
    db.commit()
    db.refresh(annotation)
    db.expunge(annotation)  # Detach from session for safe return
    
    return annotation


def get_annotations(
    entity_type: str,
    entity_id: UUID,
    db: Session,
    status: Optional[str] = None,
    position_type: Optional[str] = None,
) -> List[InlineAnnotation]:
    """Get annotations for an entity with optional filters.
    
    Args:
        entity_type: Type of entity
        entity_id: ID of the entity
        db: Database session
        status: Optional status filter (open, resolved)
        position_type: Optional position type filter
        
    Returns:
        List of InlineAnnotation instances
    """
    query = db.query(InlineAnnotation).filter(
        and_(
            InlineAnnotation.entity_type == entity_type,
            InlineAnnotation.entity_id == entity_id,
        )
    )
    
    if status is not None:
        query = query.filter(InlineAnnotation.status == status)
        
    if position_type is not None:
        query = query.filter(InlineAnnotation.position_type == position_type)
    
    # Order by created_at for consistent results
    annotations = query.order_by(InlineAnnotation.created_at).all()
    
    # Detach from session for safe return
    for annotation in annotations:
        db.expunge(annotation)
        
    return annotations


def get_annotation(annotation_id: UUID, db: Session) -> Optional[InlineAnnotation]:
    """Get single annotation by ID.
    
    Args:
        annotation_id: ID of the annotation
        db: Database session
        
    Returns:
        InlineAnnotation instance if found, None otherwise
    """
    annotation = db.query(InlineAnnotation).filter(InlineAnnotation.id == annotation_id).first()
    
    if annotation:
        db.expunge(annotation)
        
    return annotation


def resolve_annotation(annotation_id: UUID, user_id: Optional[UUID], db: Session) -> Optional[InlineAnnotation]:
    """Mark annotation as resolved.
    
    Args:
        annotation_id: ID of the annotation to resolve
        user_id: ID of user resolving the annotation (None for system)
        db: Database session
        
    Returns:
        Updated InlineAnnotation instance if found, None otherwise
    """
    annotation = db.query(InlineAnnotation).filter(InlineAnnotation.id == annotation_id).first()
    
    if not annotation:
        return None
        
    annotation.status = "resolved"
    annotation.resolved_by_id = user_id
    annotation.resolved_at = datetime.now(timezone.utc)
    
    db.commit()
    db.refresh(annotation)
    db.expunge(annotation)
    
    return annotation


def reopen_annotation(annotation_id: UUID, user_id: Optional[UUID], db: Session) -> Optional[InlineAnnotation]:
    """Reopen a resolved annotation.
    
    Args:
        annotation_id: ID of the annotation to reopen
        user_id: ID of user reopening the annotation (None for system)
        db: Database session
        
    Returns:
        Updated InlineAnnotation instance if found, None otherwise
    """
    annotation = db.query(InlineAnnotation).filter(InlineAnnotation.id == annotation_id).first()
    
    if not annotation:
        return None
        
    annotation.status = "open"
    annotation.resolved_by_id = None
    annotation.resolved_at = None
    
    db.commit()
    db.refresh(annotation)
    db.expunge(annotation)
    
    return annotation


def reply_to_annotation(
    annotation_id: UUID,
    content: str,
    user_id: Optional[UUID],
    db: Session,
) -> Optional[InlineAnnotation]:
    """Add reply to an annotation.
    
    Args:
        annotation_id: ID of the parent annotation
        content: Reply content
        user_id: ID of user creating reply (None for system)
        db: Database session
        
    Returns:
        Created reply InlineAnnotation instance if parent found, None otherwise
    """
    parent = db.query(InlineAnnotation).filter(InlineAnnotation.id == annotation_id).first()
    
    if not parent:
        return None
        
    reply = InlineAnnotation(
        entity_type=parent.entity_type,
        entity_id=parent.entity_id,
        position_type=parent.position_type,
        position_data=parent.position_data,
        content=content,
        created_by_id=user_id,
        parent_id=annotation_id,
    )
    
    db.add(reply)
    db.commit()
    db.refresh(reply)
    db.expunge(reply)
    
    return reply


def delete_annotation(annotation_id: UUID, user_id: Optional[UUID], db: Session) -> bool:
    """Delete annotation (only if user is author or user_id is None for system).
    
    Args:
        annotation_id: ID of the annotation to delete
        user_id: ID of user requesting deletion (None for system/admin)
        db: Database session
        
    Returns:
        True if deleted successfully, False if not found or unauthorized
    """
    annotation = db.query(InlineAnnotation).filter(InlineAnnotation.id == annotation_id).first()
    
    if not annotation:
        return False
        
    # Allow deletion if:
    # 1. user_id is None (system/admin can delete anything)
    # 2. user is the author (user_id == annotation.created_by_id)
    # 3. annotation was created by system (annotation.created_by_id is None) and user_id is None
    if user_id is not None and annotation.created_by_id != user_id:
        return False
        
    # Delete all replies first (cascade)
    replies = db.query(InlineAnnotation).filter(InlineAnnotation.parent_id == annotation_id).all()
    for reply in replies:
        db.delete(reply)
        
    # Delete the annotation
    db.delete(annotation)
    db.commit()
    
    return True


def get_annotation_count(entity_type: str, entity_id: UUID, db: Session) -> Dict[str, int]:
    """Get annotation counts by status for an entity.
    
    Args:
        entity_type: Type of entity
        entity_id: ID of the entity
        db: Database session
        
    Returns:
        Dictionary with status counts: {"open": 5, "resolved": 2, "total": 7}
    """
    # Count by status
    status_counts = (
        db.query(
            InlineAnnotation.status,
            func.count(InlineAnnotation.id).label('count')
        )
        .filter(
            and_(
                InlineAnnotation.entity_type == entity_type,
                InlineAnnotation.entity_id == entity_id,
            )
        )
        .group_by(InlineAnnotation.status)
        .all()
    )
    
    # Convert to dictionary
    counts = {"open": 0, "resolved": 0}
    for status, count in status_counts:
        counts[status] = count
        
    counts["total"] = counts["open"] + counts["resolved"]
    
    return counts
