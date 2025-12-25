"""QSAR API endpoints (per-target activity prediction)."""

from __future__ import annotations

from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException

from amprenta_rag.api.schemas import (
    QSARCompoundResult,
    QSARPredictRequest,
    QSARPredictResponse,
    QSARTargetInfo,
    QSARTargetPrediction,
)
from amprenta_rag.ml.qsar.predictor import TargetQSARPredictor


router = APIRouter(prefix="/qsar", tags=["qsar"])


@router.get("/targets", response_model=List[QSARTargetInfo])
def list_targets() -> List[QSARTargetInfo]:
    pred = TargetQSARPredictor()
    out = pred.list_available_targets()
    return [QSARTargetInfo.model_validate(x) for x in out]


@router.post("/predict", response_model=QSARPredictResponse)
def predict_activity(request: QSARPredictRequest) -> QSARPredictResponse:
    if len(request.smiles_list) > 100:
        raise HTTPException(status_code=400, detail="Max 100 SMILES per request")

    pred = TargetQSARPredictor()
    raw = pred.predict(request.smiles_list, targets=request.targets)
    results: List[QSARCompoundResult] = []

    for item in raw:
        preds: Dict[str, QSARTargetPrediction] = {}
        for tgt, v in (item.get("predictions") or {}).items():
            if not isinstance(v, dict):
                continue
            preds[str(tgt)] = QSARTargetPrediction.model_validate(v)
        results.append(
            QSARCompoundResult(
                smiles=str(item.get("smiles") or ""),
                predictions=preds,
                error=item.get("error"),
            )
        )

    return QSARPredictResponse(results=results)


@router.get("/targets/{target}/info")
def get_target_info(target: str) -> Dict[str, Any]:
    pred = TargetQSARPredictor()
    info = pred.get_model_info(target)
    return {"target": target, "info": info}


__all__ = ["router"]


