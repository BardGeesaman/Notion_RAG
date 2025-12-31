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


# Endpoints
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
