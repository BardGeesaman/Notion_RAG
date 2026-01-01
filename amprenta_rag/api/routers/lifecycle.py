"""Data lifecycle management API endpoints.

Provides REST API for:
- Cascade impact preview
- Status updates
- Bulk operations
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.database.models import User, LifecycleStatus
from amprenta_rag.services.lifecycle import (
    calculate_deletion_impact,
    update_lifecycle_status,
    bulk_update_status,
    bulk_delete_preview,
    execute_bulk_archive,
)

router = APIRouter(prefix="/lifecycle", tags=["lifecycle"])


# --- Schemas ---

class DeletionImpactResponse(BaseModel):
    """Response for deletion impact preview."""
    entity: dict
    impact: dict
    blocking_references: List[dict]
    can_delete: bool
    error: Optional[str] = None
    blocking_reason: Optional[str] = None


class StatusUpdateRequest(BaseModel):
    """Request to update entity lifecycle status."""
    entity_type: str = Field(..., description="Type: dataset, experiment, compound, signature")
    entity_id: UUID
    new_status: str = Field(..., description="Status: active, quarantined, invalid, archived")
    reason: Optional[str] = Field(None, max_length=500)


class StatusUpdateResponse(BaseModel):
    """Response for status update."""
    success: bool
    message: str


class BulkStatusRequest(BaseModel):
    """Request for bulk status update."""
    entity_type: str
    entity_ids: List[UUID]
    new_status: str
    reason: Optional[str] = None


class BulkStatusResponse(BaseModel):
    """Response for bulk status update."""
    total: int
    success: int
    failed: int
    errors: List[dict]


class BulkDeletePreviewRequest(BaseModel):
    """Request for bulk delete preview (dry run)."""
    entity_type: str
    entity_ids: List[UUID]


class BulkDeletePreviewResponse(BaseModel):
    """Response for bulk delete preview."""
    entity_type: str
    entity_count: int
    total_impact: dict
    blocking_entities: List[dict]
    can_proceed: bool


class BulkArchiveRequest(BaseModel):
    """Request for bulk archive operation."""
    entity_type: str
    entity_ids: List[UUID]
    reason: str = Field(..., min_length=1, max_length=500)
    confirmed: bool = Field(False, description="Must be True to execute")


# --- Endpoints ---

@router.get(
    "/impact/{entity_type}/{entity_id}",
    response_model=DeletionImpactResponse,
    summary="Preview deletion impact",
)
def get_deletion_impact(
    entity_type: str,
    entity_id: UUID,
    current_user: User = Depends(get_current_user),
) -> DeletionImpactResponse:
    """
    Preview cascade impact of deleting an entity.
    
    Returns counts of related entities and any blocking references
    that would prevent deletion.
    """
    if entity_type not in ["dataset", "experiment", "compound", "signature"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type: {entity_type}",
        )
    
    impact = calculate_deletion_impact(entity_type, entity_id)
    return DeletionImpactResponse(**impact)


@router.post(
    "/status",
    response_model=StatusUpdateResponse,
    summary="Update entity lifecycle status",
)
def update_status(
    request: StatusUpdateRequest,
    current_user: User = Depends(get_current_user),
) -> StatusUpdateResponse:
    """
    Update lifecycle status for a single entity.
    
    Creates an audit trail entry for the change.
    """
    if request.entity_type not in ["dataset", "experiment", "compound", "signature"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type: {request.entity_type}",
        )
    
    if request.new_status not in [s.value for s in LifecycleStatus]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid status: {request.new_status}",
        )
    
    success, message = update_lifecycle_status(
        entity_type=request.entity_type,
        entity_id=request.entity_id,
        new_status=request.new_status,
        reason=request.reason,
        actor_id=current_user.id,
    )
    
    if not success and "not found" in message.lower():
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=message,
        )
    
    return StatusUpdateResponse(success=success, message=message)


@router.post(
    "/bulk/status",
    response_model=BulkStatusResponse,
    summary="Bulk update lifecycle status",
)
def bulk_status_update(
    request: BulkStatusRequest,
    current_user: User = Depends(get_current_user),
) -> BulkStatusResponse:
    """
    Update lifecycle status for multiple entities at once.
    """
    if request.entity_type not in ["dataset", "experiment", "compound", "signature"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type: {request.entity_type}",
        )
    
    if request.new_status not in [s.value for s in LifecycleStatus]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid status: {request.new_status}",
        )
    
    if len(request.entity_ids) > 100:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Maximum 100 entities per bulk operation",
        )
    
    results = bulk_update_status(
        entity_type=request.entity_type,
        entity_ids=request.entity_ids,
        new_status=request.new_status,
        reason=request.reason,
        actor_id=current_user.id,
    )
    
    return BulkStatusResponse(**results)


@router.post(
    "/bulk/preview",
    response_model=BulkDeletePreviewResponse,
    summary="Preview bulk delete impact (dry run)",
)
def preview_bulk_delete(
    request: BulkDeletePreviewRequest,
    current_user: User = Depends(get_current_user),
) -> BulkDeletePreviewResponse:
    """
    Preview cascade impact for bulk deletion (dry run).
    
    Does not modify any data. Returns aggregate impact counts
    and identifies any entities with blocking references.
    """
    if request.entity_type not in ["dataset", "experiment", "compound", "signature"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type: {request.entity_type}",
        )
    
    if len(request.entity_ids) > 100:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Maximum 100 entities per preview",
        )
    
    preview = bulk_delete_preview(request.entity_type, request.entity_ids)
    return BulkDeletePreviewResponse(**preview)


@router.post(
    "/bulk/archive",
    response_model=BulkStatusResponse,
    summary="Bulk archive entities (soft delete)",
)
def bulk_archive(
    request: BulkArchiveRequest,
    current_user: User = Depends(get_current_user),
) -> BulkStatusResponse:
    """
    Bulk archive (soft delete) entities.
    
    Requires `confirmed: true` to execute. Always run preview first.
    """
    if not request.confirmed:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Set confirmed=true to execute bulk archive. Run preview first.",
        )
    
    if request.entity_type not in ["dataset", "experiment", "compound", "signature"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type: {request.entity_type}",
        )
    
    if len(request.entity_ids) > 100:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Maximum 100 entities per bulk archive",
        )
    
    results = execute_bulk_archive(
        entity_type=request.entity_type,
        entity_ids=request.entity_ids,
        reason=request.reason,
        actor_id=current_user.id,
    )
    
    return BulkStatusResponse(**results)
