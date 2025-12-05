"""
API router for Signatures.
"""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import Signature, SignatureCreate, SignatureUpdate
from amprenta_rag.api.services import signatures as signature_service

router = APIRouter()


@router.post("/", response_model=Signature, status_code=201)
async def create_signature(
    signature: SignatureCreate,
    db: Session = Depends(get_database_session),
):
    """Create a new signature."""
    return signature_service.create_signature(db, signature)


@router.get("/", response_model=List[Signature])
async def list_signatures(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    db: Session = Depends(get_database_session),
):
    """List all signatures."""
    return signature_service.get_signatures(
        db,
        skip=skip,
        limit=limit,
        name_filter=name,
        program_id=program_id,
    )


@router.get("/{signature_id}", response_model=Signature)
async def get_signature(
    signature_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Get a signature by ID."""
    signature = signature_service.get_signature(db, signature_id)
    if not signature:
        raise HTTPException(status_code=404, detail="Signature not found")
    return signature


@router.patch("/{signature_id}", response_model=Signature)
async def update_signature(
    signature_id: UUID,
    signature: SignatureUpdate,
    db: Session = Depends(get_database_session),
):
    """Update a signature."""
    updated = signature_service.update_signature(db, signature_id, signature)
    if not updated:
        raise HTTPException(status_code=404, detail="Signature not found")
    return updated


@router.delete("/{signature_id}", status_code=204)
async def delete_signature(
    signature_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Delete a signature."""
    success = signature_service.delete_signature(db, signature_id)
    if not success:
        raise HTTPException(status_code=404, detail="Signature not found")

