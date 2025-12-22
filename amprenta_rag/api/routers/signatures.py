from datetime import datetime
"""
API router for Signatures.
"""

from typing import List, Optional, cast
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import (
    AnnotationCreate,
    Signature,
    SignatureCreate,
    SignatureStatusUpdate,
    SignatureUpdate,
)
from amprenta_rag.api.services import signatures as signature_service
from amprenta_rag.database.models import Note
from typing import Any, Dict
from amprenta_rag.utils.uuid_utils import ensure_uuid

router = APIRouter()


@router.post("/", response_model=Signature, status_code=201)
async def create_signature(
    signature: SignatureCreate,
    db: Session = Depends(get_database_session),
) -> Signature:
    """Create a new signature."""
    return cast(Signature, signature_service.create_signature(db, signature))


@router.get("/", response_model=List[Signature])
async def list_signatures(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    db: Session = Depends(get_database_session),
) -> List[Signature]:
    """List all signatures."""
    return cast(
        List[Signature],
        signature_service.get_signatures(
        db,
        skip=skip,
        limit=limit,
        name_filter=name,
        program_id=program_id,
        ),
    )


@router.get("/{signature_id}", response_model=Signature)
async def get_signature(
    signature_id: UUID,
    db: Session = Depends(get_database_session),
) -> Signature:
    """Get a signature by ID."""
    signature = signature_service.get_signature(db, signature_id)
    if not signature:
        raise HTTPException(status_code=404, detail="Signature not found")
    return cast(Signature, signature)


@router.post("/{signature_id}/annotations", summary="Add annotation to signature")
async def add_signature_annotation(
    signature_id: UUID,
    annotation: AnnotationCreate,
    db: Session = Depends(get_database_session),
) -> Dict[str, Any]:
    """Add a note/annotation to a signature."""
    signature = signature_service.get_signature(db, signature_id)
    if not signature:
        raise HTTPException(status_code=404, detail="Signature not found")

    note = Note(
        entity_type="signature",
        entity_id=cast(Any, ensure_uuid(signature_id)),
        annotation_type=annotation.annotation_type,
        content=annotation.text,
    )
    db.add(note)
    db.commit()
    db.refresh(note)

    created_at_val = getattr(note, "created_at", None)
    return {
        "id": str(note.id),
        "entity_type": note.entity_type,
        "entity_id": str(note.entity_id),
        "text": note.content,
        "annotation_type": note.annotation_type,
        "created_at": created_at_val.isoformat() if isinstance(created_at_val, datetime) else None,
    }


@router.patch("/{signature_id}", response_model=Signature)
async def update_signature(
    signature_id: UUID,
    signature: SignatureUpdate,
    db: Session = Depends(get_database_session),
) -> Signature:
    """Update a signature."""
    updated = signature_service.update_signature(db, signature_id, signature)
    if not updated:
        raise HTTPException(status_code=404, detail="Signature not found")
    return cast(Signature, updated)


@router.delete("/{signature_id}", status_code=204)
async def delete_signature(
    signature_id: UUID,
    db: Session = Depends(get_database_session),
) -> None:
    """Delete a signature."""
    success = signature_service.delete_signature(db, signature_id)
    if not success:
        raise HTTPException(status_code=404, detail="Signature not found")


@router.patch("/{signature_id}/status", response_model=Signature)
async def update_signature_status(
    signature_id: UUID,
    status_update: SignatureStatusUpdate,
    db: Session = Depends(get_database_session),
) -> Signature:
    """Update signature validation status and record an audit note."""
    signature = signature_service.get_signature(db, signature_id)
    if not signature:
        raise HTTPException(status_code=404, detail="Signature not found")

    signature.validation_status = status_update.status

    note_content = status_update.reviewer_notes or f"Status changed to {status_update.status}"
    note = Note(
        entity_type="signature",
        entity_id=cast(Any, ensure_uuid(signature_id)),
        annotation_type="validation",
        content=note_content,
    )
    db.add(note)
    db.add(signature)
    db.commit()
    db.refresh(signature)

    return Signature.model_validate(signature)

