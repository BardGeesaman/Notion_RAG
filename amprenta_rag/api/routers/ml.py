"""ML Model Registry API endpoints."""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, ConfigDict

from amprenta_rag.ml.registry import get_registry

router = APIRouter(prefix="/ml", tags=["ml"])


class MLModelResponse(BaseModel):
    id: UUID
    name: str
    version: str
    model_type: str
    framework: str
    features: Optional[List[str]] = None
    metrics: Optional[Dict[str, Any]] = None
    status: str
    description: Optional[str] = None

    model_config = ConfigDict(from_attributes=True)


class PredictRequest(BaseModel):
    inputs: List[Dict[str, Any]]


@router.get("/models", response_model=List[MLModelResponse])
def list_models(
    model_type: Optional[str] = None,
    status: Optional[str] = None,
) -> List[MLModelResponse]:
    """List registered ML models."""
    registry = get_registry()
    models = registry.list_models(model_type=model_type, status=status)
    return [MLModelResponse.model_validate(m) for m in models]


@router.get("/models/{model_id}", response_model=MLModelResponse)
def get_model(model_id: UUID) -> MLModelResponse:
    """Get model details."""
    from amprenta_rag.database.models import MLModel
    from amprenta_rag.database.session import db_session

    with db_session() as db:
        model = db.query(MLModel).filter(MLModel.id == model_id).first()
        if not model:
            raise HTTPException(status_code=404, detail="Model not found")
        return MLModelResponse.model_validate(model)


@router.post("/models/{model_id}/archive", response_model=MLModelResponse)
def archive_model(model_id: UUID) -> MLModelResponse:
    """Archive a model."""
    registry = get_registry()
    try:
        model = registry.archive_model(model_id)
        return MLModelResponse.model_validate(model)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.post("/models/{model_id}/predict")
def predict(model_id: UUID, request: PredictRequest) -> Dict[str, Any]:
    """Run prediction using a registered model."""
    registry = get_registry()
    try:
        model = registry.load_model(model_id)
        if not hasattr(model, "predict"):
            raise HTTPException(status_code=400, detail="Model does not support predict()")
        predictions = model.predict(request.inputs)
        return {"predictions": predictions.tolist() if hasattr(predictions, "tolist") else list(predictions)}
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except HTTPException:
        raise
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Prediction failed: {e}")


__all__ = ["router"]


