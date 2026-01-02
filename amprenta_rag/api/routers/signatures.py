"""
Electronic signatures API endpoints.

Provides signature creation, verification, and querying.
"""

from __future__ import annotations

from typing import Any, Dict, List
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Request
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.auth.signatures import create_signature, get_signatures, verify_signature
from amprenta_rag.api.rate_limit import limiter
from amprenta_rag.auth.lockout import record_failed_attempt, is_locked_out, clear_failed_attempts

router = APIRouter()


class SignRequest(BaseModel):
    """Request for creating electronic signature."""

    user_id: UUID
    action: str
    entity_type: str
    entity_id: UUID
    meaning: str
    password: str


class SignResponse(BaseModel):
    """Response for signature creation."""

    signature_id: UUID
    signature_hash: str
    timestamp: str
    valid: bool


class SignatureInfoResponse(BaseModel):
    """Signature information response."""

    id: UUID
    user_id: UUID
    action: str
    entity_type: str
    entity_id: UUID
    signature_hash: str
    meaning: str
    timestamp: str
    verified_at: str | None


@router.post("/sign", response_model=SignResponse, status_code=201)
@limiter.limit("5/minute")
async def sign_document(
    request: Request,
    sign_request: SignRequest,
    db: Session = Depends(get_database_session),
) -> SignResponse:
    """Create electronic signature with brute force protection."""
    # Build lockout identifier from IP
    client_ip = request.client.host if request.client else "unknown"
    identifier = f"sign:{client_ip}"
    
    # Check lockout status
    locked, remaining = is_locked_out(identifier, db)
    if locked:
        raise HTTPException(
            status_code=429,
            detail=f"Too many failed attempts. Try again in {remaining} seconds.",
            headers={"Retry-After": str(remaining)},
        )
    
    # Attempt signature creation
    signature = create_signature(
        user_id=str(sign_request.user_id),
        action=sign_request.action,
        entity_type=sign_request.entity_type,
        entity_id=str(sign_request.entity_id),
        password=sign_request.password,
        meaning=sign_request.meaning,
        db=db,
    )
    
    if not signature:
        record_failed_attempt(identifier, db)
        raise HTTPException(status_code=401, detail="Authentication failed - incorrect password")
    
    # Success - clear any failed attempts
    clear_failed_attempts(identifier, db)
    
    return SignResponse(
        signature_id=signature.id,
        signature_hash=signature.signature_hash,
        timestamp=signature.timestamp.isoformat() if signature.timestamp else "",
        valid=True,
    )


@router.get("/{signature_id}/verify", response_model=Dict[str, Any])
async def verify_signature_endpoint(
    signature_id: UUID,
    db: Session = Depends(get_database_session),
) -> Dict[str, Any]:
    """Verify signature integrity."""
    is_valid = verify_signature(str(signature_id), db)
    
    return {
        "signature_id": str(signature_id),
        "valid": is_valid,
        "message": "Signature is valid" if is_valid else "Signature verification failed",
    }


@router.get("/{entity_type}/{entity_id}", response_model=List[SignatureInfoResponse])
async def list_entity_signatures(
    entity_type: str,
    entity_id: UUID,
    db: Session = Depends(get_database_session),
) -> List[SignatureInfoResponse]:
    """Get all signatures for an entity."""
    signatures = get_signatures(entity_type, str(entity_id), db)
    
    return [
        SignatureInfoResponse(
            id=sig.id,
            user_id=sig.user_id,
            action=sig.action,
            entity_type=sig.entity_type,
            entity_id=sig.entity_id,
            signature_hash=sig.signature_hash,
            meaning=sig.meaning,
            timestamp=sig.timestamp.isoformat() if sig.timestamp else "",
            verified_at=sig.verified_at.isoformat() if sig.verified_at else None,
        )
        for sig in signatures
    ]
