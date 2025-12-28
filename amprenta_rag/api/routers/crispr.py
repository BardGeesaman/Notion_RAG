"""CRISPR screen analysis API endpoints."""

from __future__ import annotations

import math
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, ConfigDict, Field

from amprenta_rag.crispr.analysis_service import run_screen_analysis
from amprenta_rag.database.models import CRISPRResult, CRISPRScreen
from amprenta_rag.database.session import db_session


router = APIRouter(prefix="/crispr", tags=["CRISPR"])


class CreateScreenRequest(BaseModel):
    name: str
    dataset_id: UUID
    library_type: Optional[str] = None
    cell_line: Optional[str] = None
    treatment: Optional[str] = None
    control_label: str
    treatment_label: str


class ScreenResponse(BaseModel):
    id: UUID
    dataset_id: UUID
    name: str
    library_type: Optional[str] = None
    cell_line: Optional[str] = None
    treatment: Optional[str] = None
    control_label: Optional[str] = None
    treatment_label: Optional[str] = None
    status: Optional[str] = None

    model_config = ConfigDict(from_attributes=True)


class AnalyzeRequest(BaseModel):
    method: str = Field("test", description="MAGeCK method (MVP: test)")


class ResultResponse(BaseModel):
    id: UUID
    screen_id: UUID
    gene_symbol: Optional[str] = None
    feature_id: Optional[UUID] = None
    beta_score: Optional[float] = None
    p_value: Optional[float] = None
    fdr: Optional[float] = None
    neg_lfc: Optional[float] = None
    pos_lfc: Optional[float] = None
    rank: Optional[int] = None
    is_hit: bool

    model_config = ConfigDict(from_attributes=True)


@router.post("/screens", response_model=ScreenResponse)
def create_screen(payload: CreateScreenRequest) -> ScreenResponse:
    with db_session() as db:
        screen = CRISPRScreen(
            dataset_id=payload.dataset_id,
            name=payload.name,
            library_type=payload.library_type,
            cell_line=payload.cell_line,
            treatment=payload.treatment,
            control_label=payload.control_label,
            treatment_label=payload.treatment_label,
            status="pending",
        )
        db.add(screen)
        db.commit()
        db.refresh(screen)
        return ScreenResponse.model_validate(screen)


@router.get("/screens", response_model=List[ScreenResponse])
def list_screens(limit: int = 100) -> List[ScreenResponse]:
    with db_session() as db:
        rows = db.query(CRISPRScreen).order_by(CRISPRScreen.created_at.desc()).limit(limit).all()
        return [ScreenResponse.model_validate(r) for r in rows]


@router.get("/screens/{screen_id}", response_model=ScreenResponse)
def get_screen(screen_id: UUID) -> ScreenResponse:
    with db_session() as db:
        screen = db.query(CRISPRScreen).filter(CRISPRScreen.id == screen_id).first()
        if not screen:
            raise HTTPException(status_code=404, detail="CRISPRScreen not found")
        return ScreenResponse.model_validate(screen)


@router.post("/screens/{screen_id}/analyze")
def analyze_screen(screen_id: UUID, payload: AnalyzeRequest) -> dict:
    with db_session() as db:
        try:
            rows = run_screen_analysis(screen_id, db, method=payload.method)
        except Exception as e:  # noqa: BLE001
            raise HTTPException(status_code=400, detail=str(e))
        return {"screen_id": str(screen_id), "results": len(rows), "status": "completed"}


@router.get("/screens/{screen_id}/results", response_model=List[ResultResponse])
def get_results(screen_id: UUID, fdr_threshold: float = 0.05, limit: int = 100) -> List[ResultResponse]:
    with db_session() as db:
        q = db.query(CRISPRResult).filter(CRISPRResult.screen_id == screen_id)
        if fdr_threshold is not None:
            q = q.filter((CRISPRResult.fdr.is_(None)) | (CRISPRResult.fdr <= float(fdr_threshold)))
        rows = q.order_by(CRISPRResult.rank.asc().nullslast()).limit(limit).all()
        return [ResultResponse.model_validate(r) for r in rows]


@router.get("/screens/{screen_id}/volcano-data")
def volcano_data(screen_id: UUID, limit: int = 5000) -> dict:
    with db_session() as db:
        rows = (
            db.query(CRISPRResult)
            .filter(CRISPRResult.screen_id == screen_id)
            .order_by(CRISPRResult.rank.asc().nullslast())
            .limit(limit)
            .all()
        )
        pts = []
        for r in rows:
            fdr = r.fdr
            y = None
            if fdr is not None and fdr > 0:
                y = -math.log10(float(fdr))
            pts.append(
                {
                    "gene": r.gene_symbol,
                    "neg_lfc": r.neg_lfc,
                    "pos_lfc": r.pos_lfc,
                    "fdr": r.fdr,
                    "neglog10_fdr": y,
                    "is_hit": bool(r.is_hit),
                }
            )
        return {"screen_id": str(screen_id), "points": pts}


__all__ = ["router"]


