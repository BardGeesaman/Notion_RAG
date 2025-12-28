"""Multi-omics latent factor analysis API endpoints."""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, ConfigDict, Field

from amprenta_rag.database.models import (
    FactorLoading,
    FactorScore,
    Feature,
    LatentFactor,
    MultiOmicsExperiment,
)
from amprenta_rag.database.session import db_session
from amprenta_rag.multi_omics.analysis_service import run_experiment


router = APIRouter(prefix="/multi-omics", tags=["Multi-Omics"])


class CreateExperimentRequest(BaseModel):
    name: str
    dataset_ids: Any = Field(..., description="JSON list/dict describing datasets/views for the experiment")
    sample_mapping: Optional[Any] = Field(None, description="JSON list describing sample alignment across views")
    n_factors: int = Field(10, ge=1, le=200)
    convergence_mode: str = Field("fast")
    description: Optional[str] = None


class ExperimentResponse(BaseModel):
    id: UUID
    name: str
    description: Optional[str] = None
    dataset_ids: Optional[Any] = None
    sample_mapping: Optional[Any] = None
    n_factors: Optional[int] = None
    convergence_mode: Optional[str] = None
    status: Optional[str] = None
    processing_log: Optional[str] = None
    created_at: Optional[datetime] = None
    processed_at: Optional[datetime] = None

    model_config = ConfigDict(from_attributes=True)


class RunResponse(BaseModel):
    experiment_id: UUID
    factors_created: int
    status: str


class FactorResponse(BaseModel):
    factor_index: int
    variance_explained: Optional[Dict[str, Any]] = None
    description: Optional[str] = None


class LoadingResponse(BaseModel):
    feature_id: UUID
    feature_name: str
    loading: float
    abs_loading: float
    omics_type: Optional[str] = None


class ScoreResponse(BaseModel):
    sample_id: str
    score: float


@router.post("/experiments", response_model=ExperimentResponse)
def create_experiment(payload: CreateExperimentRequest) -> ExperimentResponse:
    with db_session() as db:
        exp = MultiOmicsExperiment(
            name=payload.name,
            description=payload.description,
            dataset_ids=payload.dataset_ids,
            sample_mapping=payload.sample_mapping,
            n_factors=int(payload.n_factors),
            convergence_mode=str(payload.convergence_mode),
            status="pending",
            processing_log=None,
            processed_at=None,
        )
        db.add(exp)
        db.commit()
        db.refresh(exp)
        return ExperimentResponse.model_validate(exp)


@router.get("/experiments", response_model=List[ExperimentResponse])
def list_experiments(limit: int = 100) -> List[ExperimentResponse]:
    with db_session() as db:
        rows = db.query(MultiOmicsExperiment).order_by(MultiOmicsExperiment.created_at.desc()).limit(limit).all()
        return [ExperimentResponse.model_validate(r) for r in rows]


@router.get("/experiments/{experiment_id}", response_model=ExperimentResponse)
def get_experiment(experiment_id: UUID) -> ExperimentResponse:
    with db_session() as db:
        exp = db.query(MultiOmicsExperiment).filter(MultiOmicsExperiment.id == experiment_id).first()
        if not exp:
            raise HTTPException(status_code=404, detail="MultiOmicsExperiment not found")
        return ExperimentResponse.model_validate(exp)


@router.post("/experiments/{experiment_id}/run", response_model=RunResponse)
def run_multi_omics_experiment(experiment_id: UUID) -> RunResponse:
    with db_session() as db:
        try:
            factors = run_experiment(experiment_id, db)
        except Exception as e:  # noqa: BLE001
            raise HTTPException(status_code=400, detail=str(e))
        return RunResponse(experiment_id=experiment_id, factors_created=len(factors), status="completed")


def _get_factor(db, experiment_id: UUID, factor_index: int) -> LatentFactor:
    lf = (
        db.query(LatentFactor)
        .filter(LatentFactor.experiment_id == experiment_id)
        .filter(LatentFactor.factor_index == factor_index)
        .first()
    )
    if not lf:
        raise HTTPException(status_code=404, detail="LatentFactor not found")
    return lf


@router.get("/experiments/{experiment_id}/factors", response_model=List[FactorResponse])
def list_factors(experiment_id: UUID) -> List[FactorResponse]:
    with db_session() as db:
        rows = (
            db.query(LatentFactor)
            .filter(LatentFactor.experiment_id == experiment_id)
            .order_by(LatentFactor.factor_index.asc())
            .all()
        )
        return [
            FactorResponse(
                factor_index=r.factor_index,
                variance_explained=r.variance_explained,
                description=r.description,
            )
            for r in rows
        ]


@router.get(
    "/experiments/{experiment_id}/factors/{factor_index}/loadings",
    response_model=List[LoadingResponse],
)
def get_factor_loadings(experiment_id: UUID, factor_index: int, limit: int = 50) -> List[LoadingResponse]:
    with db_session() as db:
        lf = _get_factor(db, experiment_id, factor_index)
        rows = (
            db.query(FactorLoading, Feature)
            .join(Feature, Feature.id == FactorLoading.feature_id)
            .filter(FactorLoading.factor_id == lf.id)
            .order_by(FactorLoading.abs_loading.desc())
            .limit(limit)
            .all()
        )
        out: List[LoadingResponse] = []
        for fl, feat in rows:
            out.append(
                LoadingResponse(
                    feature_id=fl.feature_id,
                    feature_name=feat.name,
                    loading=float(fl.loading),
                    abs_loading=float(fl.abs_loading),
                    omics_type=fl.omics_type,
                )
            )
        return out


@router.get(
    "/experiments/{experiment_id}/factors/{factor_index}/scores",
    response_model=List[ScoreResponse],
)
def get_factor_scores(experiment_id: UUID, factor_index: int, limit: int = 5000) -> List[ScoreResponse]:
    with db_session() as db:
        lf = _get_factor(db, experiment_id, factor_index)
        rows = (
            db.query(FactorScore)
            .filter(FactorScore.factor_id == lf.id)
            .order_by(FactorScore.sample_id.asc())
            .limit(limit)
            .all()
        )
        return [ScoreResponse(sample_id=str(r.sample_id), score=float(r.score)) for r in rows]


__all__ = ["router"]


