"""ID Mapping API endpoints for lookup and management."""

from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional

from fastapi import APIRouter, Depends, HTTPException, Query, status
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.database.models import User, IDMapping, MappingRefreshLog
from amprenta_rag.database.session import db_session
from amprenta_rag.services.id_mapping_service import (
    get_mapping,
    get_mappings_batch,
    get_mapping_stats,
    get_last_successful_refresh,
)
from amprenta_rag.jobs.tasks.mapping_refresh import refresh_uniprot_mappings_task
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)
router = APIRouter(prefix="/mappings", tags=["Mappings"])


def _require_admin(current_user: User) -> None:
    """Helper to ensure current user is admin."""
    if not current_user or current_user.role != "admin":
        raise HTTPException(status_code=403, detail="Admin access required")


# Pydantic Schemas
class MappingRefreshRequest(BaseModel):
    source: str = "uniprot"  # "uniprot" | "all"


class MappingRefreshResponse(BaseModel):
    status: str
    job_id: Optional[str] = None
    message: str


class MappingStatusResponse(BaseModel):
    uniprot_last_refresh: Optional[datetime]
    total_mappings: int
    mappings_by_type: Dict[str, int]
    expired_count: int


class MappingStatsResponse(BaseModel):
    total: int
    by_source_type: Dict[str, int]
    by_target_type: Dict[str, int]
    expired: int
    permanent: int


class MappingLookupResponse(BaseModel):
    source_type: str
    source_id: str
    mappings: Dict[str, str]  # {target_type: target_id}


class BatchMappingRequest(BaseModel):
    ids: List[str]
    source_type: str
    target_type: str
    organism: str = "human"


class BatchMappingResponse(BaseModel):
    results: Dict[str, Optional[str]]  # {source_id: target_id or None}
    found: int
    not_found: int


class KEGGCacheStatusResponse(BaseModel):
    """KEGG cache status information."""
    expiring_7_days: int
    expiring_30_days: int
    total_kegg_mappings: int
    last_refresh: Optional[str]
    last_refresh_count: Optional[int]


# API Endpoints
@router.post("/refresh", response_model=MappingRefreshResponse)
def trigger_refresh(
    request: MappingRefreshRequest,
    current_user: User = Depends(get_current_user),
):
    """Trigger manual ID mapping refresh (admin only)."""
    _require_admin(current_user)
    
    logger.info(f"Admin user {current_user.id} triggered {request.source} mapping refresh")
    
    if request.source == "uniprot":
        task = refresh_uniprot_mappings_task.delay()
        return MappingRefreshResponse(
            status="queued",
            job_id=str(task.id),
            message="UniProt refresh job queued"
        )
    elif request.source == "all":
        task = refresh_uniprot_mappings_task.delay()
        return MappingRefreshResponse(
            status="queued",
            job_id=str(task.id),
            message="All mapping refresh jobs queued"
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Unknown source: {request.source}"
        )


@router.get("/status", response_model=MappingStatusResponse)
def get_status(
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
):
    """Get ID mapping status and last refresh times."""
    logger.info(f"User {current_user.id} requested mapping status")
    
    stats = get_mapping_stats()
    uniprot_last_refresh = get_last_successful_refresh("uniprot")
    
    return MappingStatusResponse(
        uniprot_last_refresh=uniprot_last_refresh,
        total_mappings=stats.get("total_mappings", 0),
        mappings_by_type=stats.get("by_source_type", {}),
        expired_count=stats.get("expired_mappings", 0)
    )


@router.get("/stats", response_model=MappingStatsResponse)
def get_stats(
    current_user: User = Depends(get_current_user),
):
    """Get detailed mapping statistics."""
    logger.info(f"User {current_user.id} requested mapping stats")
    
    stats = get_mapping_stats()
    
    return MappingStatsResponse(
        total=stats.get("total_mappings", 0),
        by_source_type=stats.get("by_source_type", {}),
        by_target_type=stats.get("by_target_type", {}),
        expired=stats.get("expired_mappings", 0),
        permanent=stats.get("permanent_mappings", 0)
    )


@router.get("/{source_type}/{source_id}", response_model=MappingLookupResponse)
def lookup_mapping(
    source_type: str,
    source_id: str,
    fallback: bool = Query(True, description="Fall back to API if not in DB"),
    current_user: User = Depends(get_current_user),
):
    """Look up all mappings for a source ID."""
    logger.info(f"User {current_user.id} looking up {source_type}:{source_id}")
    
    # Get all target types for this source
    target_types = ["uniprot", "kegg_gene", "kegg_compound"]
    mappings = {}
    
    for target_type in target_types:
        result = get_mapping(source_type, source_id, target_type, fallback=fallback)
        if result:
            mappings[target_type] = result
    
    return MappingLookupResponse(
        source_type=source_type,
        source_id=source_id,
        mappings=mappings
    )


@router.post("/batch", response_model=BatchMappingResponse)
def batch_lookup(
    request: BatchMappingRequest,
    current_user: User = Depends(get_current_user),
):
    """Batch lookup mappings (max 1000 IDs)."""
    logger.info(
        f"User {current_user.id} batch lookup {len(request.ids)} IDs "
        f"{request.source_type} -> {request.target_type}"
    )
    
    if len(request.ids) > 1000:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Maximum 1000 IDs per request"
        )
    
    results = get_mappings_batch(
        request.ids,
        request.source_type,
        request.target_type,
        request.organism
    )
    
    found = sum(1 for v in results.values() if v is not None)
    
    return BatchMappingResponse(
        results=results,
        found=found,
        not_found=len(results) - found
    )


@router.get("/kegg/status", response_model=KEGGCacheStatusResponse)
def get_kegg_cache_status(
    current_user: User = Depends(get_current_user),
) -> KEGGCacheStatusResponse:
    """Get KEGG cache status: expiring counts, last refresh timestamp."""
    now = datetime.now(timezone.utc)
    seven_days = now + timedelta(days=7)
    thirty_days = now + timedelta(days=30)
    
    with db_session() as db:
        # Count expiring mappings
        expiring_7 = db.query(IDMapping).filter(
            IDMapping.target_type.like("kegg_%"),
            IDMapping.expires_at.isnot(None),
            IDMapping.expires_at <= seven_days,
            IDMapping.expires_at > now,
        ).count()
        
        expiring_30 = db.query(IDMapping).filter(
            IDMapping.target_type.like("kegg_%"),
            IDMapping.expires_at.isnot(None),
            IDMapping.expires_at <= thirty_days,
            IDMapping.expires_at > now,
        ).count()
        
        # Total KEGG mappings
        total = db.query(IDMapping).filter(
            IDMapping.target_type.like("kegg_%"),
        ).count()
        
        # Last refresh info
        last_refresh_log = db.query(MappingRefreshLog).filter(
            MappingRefreshLog.source == "kegg_refresh",
        ).order_by(MappingRefreshLog.completed_at.desc()).first()
        
        return KEGGCacheStatusResponse(
            expiring_7_days=expiring_7,
            expiring_30_days=expiring_30,
            total_kegg_mappings=total,
            last_refresh=last_refresh_log.completed_at.isoformat() if last_refresh_log and last_refresh_log.completed_at else None,
            last_refresh_count=last_refresh_log.records_processed if last_refresh_log else None,
        )
