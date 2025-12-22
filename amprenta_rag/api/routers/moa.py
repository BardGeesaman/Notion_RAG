from __future__ import annotations

from typing import List, Dict, Any, cast
from uuid import UUID

from fastapi import APIRouter, HTTPException

from amprenta_rag.analysis.moa_inference import infer_moa
from amprenta_rag.api import schemas
from amprenta_rag.database.models import Compound
from amprenta_rag.database.session import db_session

router = APIRouter()


def _validate_compound(compound_id: UUID) -> None:
    with db_session() as db:
        exists = db.query(Compound.id).filter(Compound.id == compound_id).first()
    if not exists:
        raise HTTPException(status_code=404, detail="Compound not found")


@router.post(
    "/moa/infer",
    summary="Infer MOA for compound across datasets",
    response_model=schemas.MOAInferenceResult,
)
def infer_moa_endpoint(request: schemas.MOAInferenceRequest) -> schemas.MOAInferenceResult:
    _validate_compound(request.compound_id)
    candidates = infer_moa(request.compound_id, request.dataset_ids)
    return schemas.MOAInferenceResult(
        compound_id=request.compound_id,
        candidates=[schemas.MOACandidate(**cast(Dict[str, Any], c.asdict())) for c in candidates],
    )


@router.get(
    "/compounds/{compound_id}/moa",
    summary="Get inferred MOA for a compound",
    response_model=schemas.MOAInferenceResult,
)
def compound_moa(compound_id: UUID, dataset_ids: List[UUID] = []) -> schemas.MOAInferenceResult:
    _validate_compound(compound_id)
    candidates = infer_moa(compound_id, dataset_ids)
    return schemas.MOAInferenceResult(
        compound_id=compound_id,
        candidates=[schemas.MOACandidate(**cast(Dict[str, Any], c.asdict())) for c in candidates],
    )

