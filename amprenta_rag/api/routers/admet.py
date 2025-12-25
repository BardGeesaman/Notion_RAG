"""ADMET prediction endpoints (with optional uncertainty)."""

from __future__ import annotations

from typing import Any, Dict, List

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from amprenta_rag.api.schemas import (
    ADMETCompoundPrediction,
    ADMETEndpointPrediction,
    ADMETPredictRequest,
    ADMETPredictResponse,
)
from amprenta_rag.database.base import get_db
from amprenta_rag.ml.admet.predictor import ADMET_MODELS, get_admet_predictor
from amprenta_rag.ml.registry import get_registry


router = APIRouter(prefix="/admet", tags=["ADMET"])


def _model_info(endpoints: List[str]) -> Dict[str, Any]:
    info: Dict[str, Any] = {"endpoints": endpoints, "models": {}}

    for ep in endpoints:
        name = ADMET_MODELS.get(ep)
        if not name:
            continue
        try:
            registry = get_registry()
            rec = registry.get_active_model(name)
        except Exception as e:  # noqa: BLE001
            info["models"][ep] = {"name": name, "status": f"registry_unavailable:{type(e).__name__}"}
            continue
        if rec is None:
            info["models"][ep] = {"name": name, "status": "not_registered"}
            continue
        try:
            obj = registry.load_model(rec.id)
        except Exception:  # noqa: BLE001
            obj = None
        trained_at = obj.get("trained_at") if isinstance(obj, dict) else None
        metrics = obj.get("metrics") if isinstance(obj, dict) else None
        ensemble_size = None
        if isinstance(obj, dict):
            ensemble_size = (obj.get("ensemble") or {}).get("n_models")

        info["models"][ep] = {
            "name": name,
            "model_id": str(rec.id),
            "version": getattr(rec, "version", None),
            "trained_at": trained_at,
            "ensemble_size": ensemble_size,
            "metrics": metrics,
            "status": "active",
        }
    return info


@router.post("/predict", response_model=ADMETPredictResponse)
def predict_admet(
    request: ADMETPredictRequest,
    db: Session = Depends(get_db),  # noqa: ARG001
) -> ADMETPredictResponse:
    predictor = get_admet_predictor()
    endpoints = request.endpoints or list(ADMET_MODELS.keys())

    if request.include_uncertainty:
        raw = predictor.predict_with_uncertainty(request.smiles, endpoints=endpoints)
        results: List[ADMETCompoundPrediction] = []
        for item in raw:
            preds: Dict[str, ADMETEndpointPrediction] = {}
            for ep, v in (item.get("predictions") or {}).items():
                if not isinstance(v, dict):
                    continue
                # allow NaN mean if model missing
                if "mean" not in v:
                    continue
                preds[ep] = ADMETEndpointPrediction.model_validate(v)
            results.append(
                ADMETCompoundPrediction(
                    smiles=str(item.get("smiles") or ""),
                    predictions=preds,
                    error=item.get("error"),
                )
            )
        return ADMETPredictResponse(results=results, model_info=_model_info(endpoints))

    # Backward-compatible simple predictions -> normalize into endpoint prediction schema.
    raw2 = predictor.predict(request.smiles, endpoints=endpoints, include_shap=False)
    results2: List[ADMETCompoundPrediction] = []
    for item in raw2:
        smi = str(item.get("smiles") or "")
        err = item.get("error")
        preds: Dict[str, ADMETEndpointPrediction] = {}
        for ep in endpoints:
            v = item.get(ep)
            if not isinstance(v, dict):
                continue
            if "probability" in v:
                preds[ep] = ADMETEndpointPrediction(mean=float(v["probability"]), calibrated=False)
            elif "value" in v:
                preds[ep] = ADMETEndpointPrediction(mean=float(v["value"]), calibrated=False)
        results2.append(ADMETCompoundPrediction(smiles=smi, predictions=preds, error=err))

    return ADMETPredictResponse(results=results2, model_info=_model_info(endpoints))


__all__ = ["router"]


