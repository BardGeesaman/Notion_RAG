from __future__ import annotations

from typing import List

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field

from amprenta_rag.api import schemas
from amprenta_rag.api.services import sar_data as service
from amprenta_rag.analysis.sar_whatif import TRANSFORMATIONS, compare_properties, scaffold_hop, validate_smiles

router = APIRouter()


@router.get(
    "/targets",
    summary="List SAR targets",
    response_model=List[schemas.TargetResponse],
)
def list_targets(limit: int = Query(200, ge=1, le=2000)):
    return service.list_targets(limit=limit)


@router.get(
    "/targets/{target}/compounds",
    summary="List compounds (with activity) for a target",
    response_model=List[schemas.CompoundActivityResponse],
)
def get_compounds_by_target(target: str, limit: int = Query(2000, ge=1, le=20000)):
    return service.get_compounds_by_target(target=target, limit=limit)


@router.get(
    "/targets/{target}/cliffs",
    summary="Detect activity cliffs for a target",
    response_model=List[schemas.ActivityCliffResponse],
)
def get_activity_cliffs(
    target: str,
    similarity_threshold: float = Query(0.6, ge=0.0, le=1.0),
    fold_change: float = Query(10.0, ge=1.0, le=1e6),
    limit: int = Query(50, ge=1, le=500),
):
    return service.get_activity_cliffs_for_target(
        target=target,
        similarity_threshold=similarity_threshold,
        fold_change=fold_change,
        limit=limit,
    )


class ValidateSmilesRequest(BaseModel):
    smiles: str = Field(..., min_length=1)


class PredictRequest(BaseModel):
    smiles_list: List[str]


class ScaffoldHopRequest(BaseModel):
    smiles: str = Field(..., min_length=1)
    transformation: str = Field(..., min_length=1)


@router.post("/validate")
def sar_validate(payload: ValidateSmilesRequest) -> dict:
    try:
        ok = validate_smiles(payload.smiles)
        return {"smiles": payload.smiles, "valid": bool(ok)}
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))


@router.post("/predict")
def sar_predict(payload: PredictRequest) -> List[dict]:
    try:
        df = compare_properties(payload.smiles_list)
        return df.to_dict(orient="records")
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/transformations")
def list_transformations() -> List[dict]:
    out: List[dict] = []
    for k, v in TRANSFORMATIONS.items():
        out.append({"id": k, **v})
    return out


@router.post("/scaffold-hop")
def sar_scaffold_hop(payload: ScaffoldHopRequest) -> dict:
    try:
        products = scaffold_hop(payload.smiles, payload.transformation)
        return {"smiles": payload.smiles, "transformation": payload.transformation, "products": products}
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=str(e))


