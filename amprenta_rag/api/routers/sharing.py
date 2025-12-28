"""API router for entity sharing functionality."""

from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.api.schemas import EntityShareCreate, EntityShareResponse, EntityShareList
from amprenta_rag.models.auth import User
from amprenta_rag.services.sharing import (
    share_entity,
    unshare_entity,
    list_shares_for_entity,
    get_my_shares,
    check_share_permission,
)
from amprenta_rag.auth.permissions import can_share_entity

router = APIRouter(prefix="/sharing", tags=["sharing"])


@router.post("/entities/{entity_type}/{entity_id}/share", status_code=status.HTTP_201_CREATED)
def share_entity_endpoint(
    entity_type: str,
    entity_id: UUID,
    share_request: EntityShareCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityShareResponse:
    """
    Share an entity with a user or team.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity to share
        share_request: Share details (user_id/team_id, permission)
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        EntityShare details
    """
    # Validate entity type
    valid_types = ["dataset", "experiment", "compound", "signature"]
    if entity_type not in valid_types:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type. Must be one of: {', '.join(valid_types)}"
        )
    
    # Check authorization - user must have share permission on entity
    if not can_share_entity(str(current_user.id), entity_type, str(entity_id), db):
        raise HTTPException(status_code=403, detail="Not authorized to share this entity")
    
    share = share_entity(
        entity_type=entity_type,
        entity_id=entity_id,
        user_id=share_request.shared_with_user_id,
        team_id=share_request.shared_with_team_id,
        permission=share_request.permission,
        shared_by=current_user.id,
        db=db,
    )
    
    if not share:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Failed to create share"
        )
    
    return EntityShareResponse(
        id=share.id,
        entity_type=share.entity_type,
        entity_id=share.entity_id,
        shared_with_user_id=share.shared_with_user_id,
        shared_with_team_id=share.shared_with_team_id,
        permission=share.permission,
        shared_by_id=share.shared_by_id,
        created_at=share.created_at,
    )


@router.delete("/entities/{entity_type}/{entity_id}/share")
def unshare_entity_endpoint(
    entity_type: str,
    entity_id: UUID,
    user_id: Optional[UUID] = None,
    team_id: Optional[UUID] = None,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> dict:
    """
    Remove sharing for an entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        user_id: UUID of user to unshare from (query param)
        team_id: UUID of team to unshare from (query param)
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        Success message
    """
    # Validate entity type
    valid_types = ["dataset", "experiment", "compound", "signature"]
    if entity_type not in valid_types:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type. Must be one of: {', '.join(valid_types)}"
        )
    
    if not user_id and not team_id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Either user_id or team_id query parameter must be provided"
        )
    
    # Check authorization - user must have share permission on entity
    if not can_share_entity(str(current_user.id), entity_type, str(entity_id), db):
        raise HTTPException(status_code=403, detail="Not authorized to share this entity")
    
    success = unshare_entity(
        entity_type=entity_type,
        entity_id=entity_id,
        user_id=user_id,
        team_id=team_id,
        db=db,
    )
    
    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Share not found or failed to remove"
        )
    
    return {"message": "Share removed successfully"}


@router.get("/entities/{entity_type}/{entity_id}/shares")
def list_entity_shares(
    entity_type: str,
    entity_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityShareList:
    """
    List all shares for a specific entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        List of shares for the entity
    """
    # Validate entity type
    valid_types = ["dataset", "experiment", "compound", "signature"]
    if entity_type not in valid_types:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type. Must be one of: {', '.join(valid_types)}"
        )
    
    # Check authorization - user must have share permission on entity
    if not can_share_entity(str(current_user.id), entity_type, str(entity_id), db):
        raise HTTPException(status_code=403, detail="Not authorized to share this entity")
    
    shares = list_shares_for_entity(
        entity_type=entity_type,
        entity_id=entity_id,
        db=db,
    )
    
    share_responses = [
        EntityShareResponse(
            id=share.id,
            entity_type=share.entity_type,
            entity_id=share.entity_id,
            shared_with_user_id=share.shared_with_user_id,
            shared_with_team_id=share.shared_with_team_id,
            permission=share.permission,
            shared_by_id=share.shared_by_id,
            created_at=share.created_at,
        )
        for share in shares
    ]
    
    return EntityShareList(shares=share_responses)


@router.get("/my-shares")
def get_my_shares_endpoint(
    entity_type: Optional[str] = None,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> EntityShareList:
    """
    Get all entities shared with the current user.
    
    Args:
        entity_type: Optional filter by entity type (query param)
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        List of entities shared with the user
    """
    if entity_type:
        valid_types = ["dataset", "experiment", "compound", "signature"]
        if entity_type not in valid_types:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Invalid entity_type. Must be one of: {', '.join(valid_types)}"
            )
    
    shares = get_my_shares(
        user_id=current_user.id,
        entity_type=entity_type,
        db=db,
    )
    
    share_responses = [
        EntityShareResponse(
            id=share.id,
            entity_type=share.entity_type,
            entity_id=share.entity_id,
            shared_with_user_id=share.shared_with_user_id,
            shared_with_team_id=share.shared_with_team_id,
            permission=share.permission,
            shared_by_id=share.shared_by_id,
            created_at=share.created_at,
        )
        for share in shares
    ]
    
    return EntityShareList(shares=share_responses)


@router.get("/entities/{entity_type}/{entity_id}/check-permission")
def check_entity_permission(
    entity_type: str,
    entity_id: UUID,
    permission: str = "view",
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> dict:
    """
    Check if current user has permission to access an entity.
    
    Args:
        entity_type: Type of entity (dataset, experiment, compound, signature)
        entity_id: UUID of the entity
        permission: Required permission level (view, edit, admin) - query param
        current_user: Current authenticated user
        db: Database session
    
    Returns:
        Permission check result
    """
    # Validate entity type
    valid_types = ["dataset", "experiment", "compound", "signature"]
    if entity_type not in valid_types:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid entity_type. Must be one of: {', '.join(valid_types)}"
        )
    
    # Validate permission level
    valid_permissions = ["view", "edit", "admin"]
    if permission not in valid_permissions:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid permission. Must be one of: {', '.join(valid_permissions)}"
        )
    
    has_permission = check_share_permission(
        user_id=current_user.id,
        entity_type=entity_type,
        entity_id=entity_id,
        required_permission=permission,
        db=db,
    )
    
    return {
        "entity_type": entity_type,
        "entity_id": str(entity_id),
        "permission": permission,
        "has_permission": has_permission,
        "user_id": str(current_user.id),
    }
