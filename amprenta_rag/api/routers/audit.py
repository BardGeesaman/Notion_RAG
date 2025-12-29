"""
Audit trail API endpoints.

Provides audit log querying and integrity verification.
"""

from __future__ import annotations

from typing import Any, Dict, List

from fastapi import APIRouter, Depends, Query
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.auth.audit import get_audit_trail, verify_integrity
from amprenta_rag.database.models import AuditLog

router = APIRouter(prefix="/audit", tags=["Audit"])


class AuditLogResponse(BaseModel):
    """Audit log entry response."""

    entity_type: str | None
    entity_id: str | None
    action: str
    username: str | None
    timestamp: str
    old_checksum: str | None
    new_checksum: str | None
    changes: Dict[str, Any] | None


class VerifyIntegrityRequest(BaseModel):
    """Request for integrity verification."""

    current_data: Dict[str, Any]


class VerifyIntegrityResponse(BaseModel):
    """Response for integrity verification."""

    verified: bool
    message: str


@router.get("/{entity_type}/{entity_id}", response_model=List[AuditLogResponse])
async def get_audit_trail_endpoint(
    entity_type: str,
    entity_id: str,
    db: Session = Depends(get_database_session),
) -> List[AuditLogResponse]:
    """Get audit trail for an entity."""
    logs = get_audit_trail(db, entity_type, entity_id)
    
    return [
        AuditLogResponse(
            entity_type=log.entity_type,
            entity_id=log.entity_id,
            action=log.action,
            username=log.username,
            timestamp=log.timestamp.isoformat() if log.timestamp else "",
            old_checksum=log.old_checksum,
            new_checksum=log.new_checksum,
            changes=log.changes,
        )
        for log in logs
    ]


@router.post("/{entity_type}/{entity_id}/verify", response_model=VerifyIntegrityResponse)
async def verify_integrity_endpoint(
    entity_type: str,
    entity_id: str,
    request: VerifyIntegrityRequest,
    db: Session = Depends(get_database_session),
) -> VerifyIntegrityResponse:
    """Verify data integrity against audit log."""
    verified = verify_integrity(db, entity_type, entity_id, request.current_data)
    
    message = "Data integrity verified" if verified else "Data integrity check failed - checksum mismatch"
    
    return VerifyIntegrityResponse(verified=verified, message=message)


@router.get("/recent", response_model=List[AuditLogResponse])
async def get_recent_changes(
    limit: int = Query(50, ge=1, le=200),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_database_session),
) -> List[AuditLogResponse]:
    """Get recent audit log entries."""
    logs = (
        db.query(AuditLog)
        .order_by(AuditLog.timestamp.desc())
        .offset(offset)
        .limit(limit)
        .all()
    )
    
    return [
        AuditLogResponse(
            entity_type=log.entity_type,
            entity_id=log.entity_id,
            action=log.action,
            username=log.username,
            timestamp=log.timestamp.isoformat() if log.timestamp else "",
            old_checksum=log.old_checksum,
            new_checksum=log.new_checksum,
            changes=log.changes,
        )
        for log in logs
    ]


@router.get("/user/{user_id}", response_model=List[AuditLogResponse])
async def get_user_changes(
    user_id: str,
    limit: int = Query(50, ge=1, le=200),
    db: Session = Depends(get_database_session),
) -> List[AuditLogResponse]:
    """Get changes by specific user."""
    logs = (
        db.query(AuditLog)
        .filter(AuditLog.user_id == user_id)
        .order_by(AuditLog.timestamp.desc())
        .limit(limit)
        .all()
    )
    
    return [
        AuditLogResponse(
            entity_type=log.entity_type,
            entity_id=log.entity_id,
            action=log.action,
            username=log.username,
            timestamp=log.timestamp.isoformat() if log.timestamp else "",
            old_checksum=log.old_checksum,
            new_checksum=log.new_checksum,
            changes=log.changes,
        )
        for log in logs
    ]

