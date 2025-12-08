from __future__ import annotations

from typing import List

from fastapi import APIRouter, HTTPException

from amprenta_rag.api import schemas
from amprenta_rag.api.services import compounds as service

router = APIRouter()


@router.get("/", summary="List compounds", response_model=List[schemas.CompoundResponse])
def list_compounds():
    return service.list_compounds()


@router.get(
    "/{compound_id}",
    summary="Get compound by ID",
    response_model=schemas.CompoundResponse,
)
def get_compound(compound_id: str):
    compound = service.get_compound_by_id(compound_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    return compound


@router.get(
    "/{compound_id}/programs",
    summary="Get programs linked to compound",
    response_model=List[schemas.ProgramLinkResponse],
)
def get_compound_programs(compound_id: str):
    return service.get_compound_programs(compound_id)

