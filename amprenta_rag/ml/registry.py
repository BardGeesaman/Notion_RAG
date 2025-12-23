"""ML Model Registry service."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import UUID

import joblib

from amprenta_rag.database.models import MLModel
from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Default artifact storage path
ARTIFACT_BASE_PATH = Path(os.getenv("ML_ARTIFACT_PATH", "models/"))


class MLModelRegistry:
    """Service for managing ML models."""

    def __init__(self, artifact_base: Optional[Path] = None):
        self.artifact_base = artifact_base or ARTIFACT_BASE_PATH
        self.artifact_base.mkdir(parents=True, exist_ok=True)
        self._model_cache: Dict[str, Any] = {}

    def register_model(
        self,
        name: str,
        version: str,
        model_type: str,
        framework: str,
        model_object: Any,
        features: Optional[List[str]] = None,
        hyperparameters: Optional[Dict[str, Any]] = None,
        metrics: Optional[Dict[str, float]] = None,
        description: Optional[str] = None,
        training_dataset_id: Optional[UUID] = None,
    ) -> MLModel:
        """Register a new model with artifact storage."""
        # Save artifact
        artifact_filename = f"{name}_{version}.joblib"
        artifact_path = self.artifact_base / artifact_filename
        joblib.dump(model_object, artifact_path)

        with db_session() as db:
            ml_model = MLModel(
                name=name,
                version=version,
                model_type=model_type,
                framework=framework,
                artifact_path=str(artifact_path),
                training_dataset_id=training_dataset_id,
                features=features,
                hyperparameters=hyperparameters,
                metrics=metrics,
                description=description,
                status="active",
            )
            db.add(ml_model)
            db.commit()
            db.refresh(ml_model)
            logger.info("[ML] Registered model %s:%s (%s)", name, version, ml_model.id)
            return ml_model

    def load_model(self, model_id: UUID) -> Any:
        """Load model from artifact storage (with caching)."""
        cache_key = str(model_id)
        if cache_key in self._model_cache:
            return self._model_cache[cache_key]

        with db_session() as db:
            ml_model = db.query(MLModel).filter(MLModel.id == model_id).first()
            if not ml_model:
                raise ValueError(f"Model {model_id} not found")

            model_obj = joblib.load(ml_model.artifact_path)
            self._model_cache[cache_key] = model_obj
            return model_obj

    def get_active_model(self, name: str) -> Optional[MLModel]:
        """Get the latest active model by name."""
        with db_session() as db:
            return (
                db.query(MLModel)
                .filter(MLModel.name == name, MLModel.status == "active")
                .order_by(MLModel.created_at.desc())
                .first()
            )

    def list_models(
        self,
        model_type: Optional[str] = None,
        status: Optional[str] = None,
    ) -> List[MLModel]:
        """List models with optional filtering."""
        with db_session() as db:
            query = db.query(MLModel)
            if model_type:
                query = query.filter(MLModel.model_type == model_type)
            if status:
                query = query.filter(MLModel.status == status)
            return query.order_by(MLModel.created_at.desc()).all()

    def archive_model(self, model_id: UUID) -> MLModel:
        """Archive a model (soft delete)."""
        with db_session() as db:
            ml_model = db.query(MLModel).filter(MLModel.id == model_id).first()
            if not ml_model:
                raise ValueError(f"Model {model_id} not found")
            ml_model.status = "archived"
            db.commit()
            db.refresh(ml_model)
            self._model_cache.pop(str(model_id), None)
            return ml_model


# Singleton instance
_registry: Optional[MLModelRegistry] = None


def get_registry() -> MLModelRegistry:
    global _registry
    if _registry is None:
        _registry = MLModelRegistry()
    return _registry


__all__ = ["MLModelRegistry", "get_registry"]


