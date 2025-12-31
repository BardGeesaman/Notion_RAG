"""Version management API endpoints."""

from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.auth.versioning import (
    create_version,
    get_versions,
    get_version,
    compare_versions,
    rollback_to_version,
    VERSIONABLE_ENTITIES,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)
router = APIRouter(prefix="/versions", tags=["Versions"])


# Request/Response Schemas
class VersionResponse(BaseModel):
    id: UUID
    entity_type: str
    entity_id: UUID
    version_number: int
    checksum_sha256: str
    created_by: Optional[UUID]
    created_at: str
    change_summary: Optional[str]


class VersionDetailResponse(VersionResponse):
    data_snapshot: Dict[str, Any]


class CreateVersionRequest(BaseModel):
    data: Dict[str, Any]
    change_summary: Optional[str] = None


class CompareVersionsRequest(BaseModel):
    version_a: int
    version_b: int


class VersionDiffResponse(BaseModel):
    version_a: int
    version_b: int
    added: Dict[str, Any]
    removed: Dict[str, Any]
    changed: Dict[str, Any]


class RestoreVersionRequest(BaseModel):
    confirm: bool  # P2: Must be True to proceed
    reason: Optional[str] = None  # Optional restore reason


class RestoreVersionResponse(BaseModel):
    success: bool
    restored_from_version: int
    new_version_number: int
    entity_type: str
    entity_id: UUID
    message: str


# Endpoints
# Note: More specific routes must come before generic ones
@router.post("/restore/{version_id}", response_model=RestoreVersionResponse)
def restore_version(
    version_id: UUID,
    request: RestoreVersionRequest,
    db: Session = Depends(get_database_session),
    current_user = Depends(get_current_user),
) -> RestoreVersionResponse:
    """
    Restore entity to a previous version.
    
    P2 Safeguards:
    - Requires confirm=True in request body
    - Admin role required
    - Creates audit log entry
    """
    # P2: Require confirmation
    if not request.confirm:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Restore operation requires confirm=True"
        )
    
    # P2: Admin-only access
    if current_user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Admin access required for version restore"
        )
    
    # Get version to restore
    version = get_version(db, version_id)
    if not version:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Version {version_id} not found"
        )
    
    # Create new version from old snapshot using rollback_to_version
    change_summary = f"Restored from version {version.version_number}"
    if request.reason:
        change_summary += f": {request.reason}"
    
    try:
        new_version = rollback_to_version(
            db=db,
            entity_type=version.entity_type,
            entity_id=version.entity_id,
            target_version_number=version.version_number,
            user_id=current_user.id,
            reason=request.reason,
        )
        
        # P2: Log audit entry
        from amprenta_rag.auth.audit import log_action
        log_action(
            action="restore_version",
            user_id=str(current_user.id),
            username=getattr(current_user, 'username', 'unknown'),
            entity_type=version.entity_type,
            entity_id=str(version.entity_id),
            details={
                "restored_from_version": version.version_number,
                "new_version_number": new_version.version_number,
                "reason": request.reason,
            }
        )
        
        logger.info(
            f"User {getattr(current_user, 'username', current_user.id)} restored "
            f"{version.entity_type}/{version.entity_id} "
            f"from v{version.version_number} to v{new_version.version_number}"
        )
        
        return RestoreVersionResponse(
            success=True,
            restored_from_version=version.version_number,
            new_version_number=new_version.version_number,
            entity_type=version.entity_type,
            entity_id=version.entity_id,
            message=change_summary,
        )
        
    except Exception as e:
        logger.error(f"Failed to restore version {version_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to restore version: {str(e)}"
        )


@router.get("/{entity_type}/{entity_id}", response_model=List[VersionResponse])
def list_versions(
    entity_type: str,
    entity_id: UUID,
    db: Session = Depends(get_database_session),
    current_user = Depends(get_current_user),
) -> List[VersionResponse]:
    """List all versions for an entity."""
    if entity_type not in VERSIONABLE_ENTITIES:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Entity type '{entity_type}' is not versionable"
        )
    
    versions = get_versions(db, entity_type, entity_id)
    
    return [
        VersionResponse(
            id=v.id,
            entity_type=v.entity_type,
            entity_id=v.entity_id,
            version_number=v.version_number,
            checksum_sha256=v.checksum_sha256,
            created_by=v.created_by,
            created_at=v.created_at.isoformat(),
            change_summary=v.change_summary,
        )
        for v in versions
    ]


@router.get("/{version_id}", response_model=VersionDetailResponse)
def get_version_detail(
    version_id: UUID,
    db: Session = Depends(get_database_session),
    current_user = Depends(get_current_user),
) -> VersionDetailResponse:
    """Get a specific version with full snapshot."""
    version = get_version(db, version_id)
    
    if not version:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Version {version_id} not found"
        )
    
    return VersionDetailResponse(
        id=version.id,
        entity_type=version.entity_type,
        entity_id=version.entity_id,
        version_number=version.version_number,
        checksum_sha256=version.checksum_sha256,
        created_by=version.created_by,
        created_at=version.created_at.isoformat(),
        change_summary=version.change_summary,
        data_snapshot=version.data_snapshot,
    )


@router.post("/{entity_type}/{entity_id}", response_model=VersionResponse, status_code=201)
def create_manual_version(
    entity_type: str,
    entity_id: UUID,
    request: CreateVersionRequest,
    db: Session = Depends(get_database_session),
    current_user = Depends(get_current_user),
) -> VersionResponse:
    """Create a manual version snapshot."""
    if entity_type not in VERSIONABLE_ENTITIES:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Entity type '{entity_type}' is not versionable"
        )
    
    try:
        version = create_version(
            db=db,
            entity_type=entity_type,
            entity_id=entity_id,
            data=request.data,
            user_id=current_user.id,
            change_summary=request.change_summary,
        )
        
        return VersionResponse(
            id=version.id,
            entity_type=version.entity_type,
            entity_id=version.entity_id,
            version_number=version.version_number,
            checksum_sha256=version.checksum_sha256,
            created_by=version.created_by,
            created_at=version.created_at.isoformat(),
            change_summary=version.change_summary,
        )
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        logger.error(f"Failed to create version: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to create version"
        )


@router.post("/{entity_type}/{entity_id}/compare", response_model=VersionDiffResponse)
def compare_entity_versions(
    entity_type: str,
    entity_id: UUID,
    request: CompareVersionsRequest,
    db: Session = Depends(get_database_session),
    current_user = Depends(get_current_user),
) -> VersionDiffResponse:
    """Compare two versions of an entity."""
    if entity_type not in VERSIONABLE_ENTITIES:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Entity type '{entity_type}' is not versionable"
        )
    
    # Get both versions
    versions = get_versions(db, entity_type, entity_id)
    version_map = {v.version_number: v for v in versions}
    
    version_a = version_map.get(request.version_a)
    version_b = version_map.get(request.version_b)
    
    if not version_a:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Version {request.version_a} not found"
        )
    
    if not version_b:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Version {request.version_b} not found"
        )
    
    # Compare versions
    diff = compare_versions(version_a.data_snapshot, version_b.data_snapshot)
    
    return VersionDiffResponse(
        version_a=request.version_a,
        version_b=request.version_b,
        added=diff["added"],
        removed=diff["removed"],
        changed=diff["changed"],
    )
