from __future__ import annotations

from fastapi import APIRouter, HTTPException

from amprenta_rag.api.services import compounds as service

router = APIRouter()


@router.get("/", summary="List compounds")
def list_compounds():
    return service.list_compounds()


@router.get("/{compound_id}", summary="Get compound by ID")
def get_compound(compound_id: str):
    compound = service.get_compound_by_id(compound_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    return compound


@router.get("/{compound_id}/programs", summary="Get programs linked to compound")
def get_compound_programs(compound_id: str):
    return service.get_compound_programs(compound_id)

