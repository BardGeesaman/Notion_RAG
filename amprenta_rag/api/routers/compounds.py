from __future__ import annotations

from typing import List

from fastapi import APIRouter, HTTPException

from amprenta_rag.api import schemas
from amprenta_rag.api.services import compounds as service
from amprenta_rag.database.models import Compound, Note
from amprenta_rag.database.session import db_session

router = APIRouter()


@router.get("/", summary="List compounds", response_model=List[schemas.CompoundResponse])
def list_compounds() -> List[schemas.CompoundResponse]:
    return [schemas.CompoundResponse(**c) for c in service.list_compounds()]


@router.get(
    "/{compound_id}",
    summary="Get compound by ID",
    response_model=schemas.CompoundResponse,
)
def get_compound(compound_id: str) -> schemas.CompoundResponse:
    compound = service.get_compound_by_id(compound_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    return schemas.CompoundResponse(**compound)


@router.get(
    "/{compound_id}/programs",
    summary="Get programs linked to compound",
    response_model=List[schemas.ProgramLinkResponse],
)
def get_compound_programs(compound_id: str) -> List[schemas.ProgramLinkResponse]:
    return [schemas.ProgramLinkResponse(**p) for p in service.get_compound_programs(compound_id)]


@router.post(
    "/{compound_id}/annotations",
    summary="Add annotation to compound",
)
def add_compound_annotation(compound_id: str, annotation: schemas.AnnotationCreate) -> dict:
    """Add a note/annotation to a compound (by compound_id string)."""
    with db_session() as db:
        compound = db.query(Compound).filter(Compound.compound_id == compound_id).first()
        if not compound:
            raise HTTPException(status_code=404, detail="Compound not found")

        note = Note(
            entity_type="compound",
            entity_id=compound.id,
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
            "created_at": created_at_val.isoformat() if created_at_val else None,
        }

