"""ADMET prediction endpoints (with optional uncertainty)."""

from __future__ import annotations

import asyncio
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from amprenta_rag.api.schemas import (
    ADMETCompoundPrediction,
    ADMETExplainRequest,
    ADMETExplainResponse,
    ADMETEndpointPrediction,
    ADMETPredictRequest,
    ADMETPredictResponse,
)
from amprenta_rag.database.base import get_db
from amprenta_rag.ml.admet.predictor import ADMET_MODELS, get_admet_predictor
from amprenta_rag.ml.registry import get_registry
from amprenta_rag.ml.gnn.predictor import get_gnn_predictor
from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.models.auth import User
from fastapi import HTTPException, Query


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


# Sync helper functions for ML model inference
def _sync_predict_admet(request: ADMETPredictRequest, db: Session) -> ADMETPredictResponse:
    """Sync helper for ADMET ML model inference."""
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


def _sync_explain_admet(request: ADMETExplainRequest, db: Session) -> ADMETExplainResponse:
    """Sync helper for ADMET prediction with SHAP analysis."""
    predictor = get_admet_predictor()
    endpoint = str(request.endpoint or "herg")

    raw = predictor.predict_with_uncertainty(
        [request.smiles],
        endpoints=[endpoint],
        include_shap=True,
        shap_top_k=int(request.top_k),
    )
    item = (raw or [{}])[0] if isinstance(raw, list) else {}
    err: Optional[str] = item.get("error") if isinstance(item, dict) else "prediction_failed"

    preds = item.get("predictions") if isinstance(item, dict) else None
    pred_ep = preds.get(endpoint) if isinstance(preds, dict) else None

    shap: Optional[Dict[str, Any]] = None
    prediction: Dict[str, Any] = {}

    if isinstance(pred_ep, dict):
        shap = pred_ep.get("shap")
        prediction = dict(pred_ep)
        prediction.pop("shap", None)
        if "error" in prediction and not err:
            err = str(prediction.get("error"))
    else:
        if err is None:
            err = "prediction_unavailable"

    return ADMETExplainResponse(
        smiles=str(request.smiles or ""),
        endpoint=endpoint,
        prediction=prediction,
        shap=shap,
        error=err,
    )


@router.post("/predict", response_model=ADMETPredictResponse)
async def predict_admet(
    request: ADMETPredictRequest,
    db: Session = Depends(get_db),
) -> ADMETPredictResponse:
    """Predict ADMET properties using async thread pool for ML inference."""
    return await asyncio.to_thread(_sync_predict_admet, request, db)


@router.post("/explain", response_model=ADMETExplainResponse)
async def explain_admet_prediction(
    request: ADMETExplainRequest,
    db: Session = Depends(get_db),
) -> ADMETExplainResponse:
    """Explain ADMET prediction using async thread pool for ML + SHAP analysis."""
    return await asyncio.to_thread(_sync_explain_admet, request, db)


# ============ GNN Toxicity Endpoints ============

GNN_ENDPOINTS = ["herg", "ames", "dili", "ld50", "clintox"]


@router.post(
    "/predict-gnn",
    summary="Predict using GNN models",
    description="Deep learning toxicity predictions using Graph Neural Networks."
)
def predict_gnn(
    request: ADMETPredictRequest,
    endpoints: List[str] = Query(default=["herg"], description="GNN endpoints to run"),
    with_uncertainty: bool = Query(default=True, description="Include MC Dropout uncertainty"),
    current_user: User = Depends(get_current_user)
):
    """Predict toxicity using GNN models.
    
    Uses Graph Neural Networks trained on TDC datasets for molecular
    toxicity prediction with uncertainty quantification.
    """
    results = []
    
    for smiles in request.smiles_list:
        pred = {"smiles": smiles, "endpoints": {}}
        
        for endpoint in endpoints:
            if endpoint not in GNN_ENDPOINTS:
                pred["endpoints"][endpoint] = {"error": f"Unknown endpoint: {endpoint}"}
                continue
            
            predictor = get_gnn_predictor(endpoint)
            if predictor is None:
                pred["endpoints"][endpoint] = {"error": "Model not available"}
                continue
            
            try:
                result = predictor.predict([smiles], with_uncertainty=with_uncertainty)[0]
                pred["endpoints"][endpoint] = {
                    k: v for k, v in result.items() if k != "smiles"
                }
            except Exception as e:
                pred["endpoints"][endpoint] = {"error": str(e)}
        
        results.append(pred)
    
    return {"predictions": results}


@router.get(
    "/gnn-models",
    summary="List GNN models",
    description="List available GNN toxicity models."
)
def list_gnn_models(current_user: User = Depends(get_current_user)):
    """List available GNN models.
    
    Returns status and performance metrics for each GNN toxicity endpoint.
    """
    models = []
    for endpoint in GNN_ENDPOINTS:
        model_path = Path(f"models/gnn/gnn_{endpoint}.pt")
        if model_path.exists():
            try:
                checkpoint = torch.load(model_path, map_location="cpu")
                models.append({
                    "endpoint": endpoint,
                    "available": True,
                    "metrics": checkpoint.get("metrics", {}),
                    "trained_at": checkpoint.get("trained_at"),
                    "task_type": checkpoint.get("config", {}).get("task_type", "unknown")
                })
            except Exception as e:
                models.append({
                    "endpoint": endpoint,
                    "available": False,
                    "error": str(e)
                })
        else:
            models.append({
                "endpoint": endpoint,
                "available": False,
                "error": "Model not trained"
            })
    
    return {"models": models}


@router.get(
    "/gnn-models/{endpoint}/info",
    summary="Get GNN model info",
    description="Get detailed information about a specific GNN model."
)
def get_gnn_model_info(
    endpoint: str,
    current_user: User = Depends(get_current_user)
):
    """Get GNN model metadata.
    
    Returns detailed information about model architecture,
    training metrics, and configuration.
    """
    if endpoint not in GNN_ENDPOINTS:
        raise HTTPException(status_code=404, detail=f"Unknown endpoint: {endpoint}")
    
    model_path = Path(f"models/gnn/gnn_{endpoint}.pt")
    if not model_path.exists():
        raise HTTPException(status_code=404, detail=f"Model not trained: {endpoint}")
    
    try:
        checkpoint = torch.load(model_path, map_location="cpu")
        
        return {
            "endpoint": endpoint,
            "config": checkpoint.get("config", {}),
            "metrics": checkpoint.get("metrics", {}),
            "trained_at": checkpoint.get("trained_at"),
            "model_path": str(model_path),
            "file_size_mb": model_path.stat().st_size / (1024 * 1024)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to load model info: {str(e)}")


@router.post(
    "/compare",
    summary="Compare XGBoost vs GNN",
    description="Compare predictions from traditional and deep learning models."
)
def compare_models(
    request: ADMETPredictRequest,
    endpoint: str = Query(default="herg", description="Endpoint to compare"),
    current_user: User = Depends(get_current_user)
):
    """Compare XGBoost ensemble vs GNN predictions.
    
    Provides side-by-side comparison of traditional ML (XGBoost ensemble)
    vs deep learning (GNN) approaches for toxicity prediction.
    """
    from amprenta_rag.ml.admet.predictor import ADMETPredictor
    
    # Get XGBoost predictions
    try:
        xgb_predictor = ADMETPredictor()
        xgb_results = xgb_predictor.predict(request.smiles_list, endpoints=[endpoint])
    except Exception as e:
        xgb_results = [{"error": f"XGBoost prediction failed: {str(e)}"}] * len(request.smiles_list)
    
    # Get GNN predictions
    gnn_predictor = get_gnn_predictor(endpoint)
    
    comparisons = []
    for i, smiles in enumerate(request.smiles_list):
        comparison = {"smiles": smiles}
        
        # XGBoost result
        if i < len(xgb_results) and endpoint in xgb_results[i]:
            comparison["xgboost"] = xgb_results[i][endpoint]
        else:
            comparison["xgboost"] = {"error": "Not available"}
        
        # GNN result
        if gnn_predictor:
            try:
                gnn_result = gnn_predictor.predict([smiles])[0]
                comparison["gnn"] = {k: v for k, v in gnn_result.items() if k != "smiles"}
            except Exception as e:
                comparison["gnn"] = {"error": str(e)}
        else:
            comparison["gnn"] = {"error": "Model not available"}
        
        comparisons.append(comparison)
    
    return {"endpoint": endpoint, "comparisons": comparisons}


__all__ = ["router"]


