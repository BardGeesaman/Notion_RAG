"""Connectivity Mapping (LINCS / CMap) API endpoints."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, ConfigDict, Field

from amprenta_rag.connectivity.ingest_service import ingest_lincs_data
from amprenta_rag.connectivity.connectivity_service import ConnectivityService
from amprenta_rag.database.models import ConnectivityScore, LINCSSignature
from amprenta_rag.database.session import db_session


router = APIRouter(prefix="/connectivity", tags=["Connectivity"])

SERVICE = ConnectivityService()


class IngestRequest(BaseModel):
    max_compounds: int = Field(default=1000, ge=1, le=200000)


class IngestResponse(BaseModel):
    ingested_signatures: int


class LINCSSignatureResponse(BaseModel):
    id: UUID
    sig_id: str
    pert_iname: Optional[str] = None
    pert_id: Optional[str] = None
    pert_type: Optional[str] = None
    cell_id: Optional[str] = None
    tas: Optional[float] = None

    model_config = ConfigDict(from_attributes=True)


class ConnectivityScoreResponse(BaseModel):
    lincs_signature_id: UUID
    sig_id: str
    pert_iname: Optional[str] = None
    pert_id: Optional[str] = None
    cell_id: Optional[str] = None
    score: float


@router.post("/ingest", response_model=IngestResponse)
def ingest(payload: IngestRequest) -> IngestResponse:
    n = ingest_lincs_data(max_compounds=int(payload.max_compounds))
    return IngestResponse(ingested_signatures=int(n))


@router.get("/lincs-signatures", response_model=List[LINCSSignatureResponse])
def list_lincs_signatures(
    pert_id: Optional[str] = Query(None),
    pert_iname: Optional[str] = Query(None),
    cell_id: Optional[str] = Query(None),
    limit: int = Query(100, ge=1, le=1000),
) -> List[LINCSSignatureResponse]:
    with db_session() as db:
        q = db.query(LINCSSignature).order_by(LINCSSignature.ingested_at.desc())
        if pert_id:
            q = q.filter(LINCSSignature.pert_id == pert_id)
        if pert_iname:
            q = q.filter(LINCSSignature.pert_iname == pert_iname)
        if cell_id:
            q = q.filter(LINCSSignature.cell_id == cell_id)
        rows = q.limit(limit).all()
        return [LINCSSignatureResponse.model_validate(r) for r in rows]


@router.post("/compute/{signature_id}")
def compute(signature_id: UUID) -> dict:
    try:
        n = SERVICE.compute_scores(signature_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"computed": int(n)}


@router.get("/signatures/{signature_id}/top-reversals", response_model=List[ConnectivityScoreResponse])
def top_reversals(signature_id: UUID, n: int = Query(50, ge=1, le=500)) -> List[ConnectivityScoreResponse]:
    with db_session() as db:
        rows = (
            db.query(ConnectivityScore, LINCSSignature)
            .join(LINCSSignature, LINCSSignature.id == ConnectivityScore.lincs_signature_id)
            .filter(ConnectivityScore.query_signature_id == signature_id)
            .filter(ConnectivityScore.score < 0)
            .order_by(ConnectivityScore.score.asc())
            .limit(n)
            .all()
        )
        out: List[ConnectivityScoreResponse] = []
        for s, lincs_sig in rows:
            out.append(
                ConnectivityScoreResponse(
                    lincs_signature_id=lincs_sig.id,
                    sig_id=lincs_sig.sig_id,
                    pert_iname=lincs_sig.pert_iname,
                    pert_id=lincs_sig.pert_id,
                    cell_id=lincs_sig.cell_id,
                    score=float(s.score),
                )
            )
        return out


@router.get("/signatures/{signature_id}/top-mimics", response_model=List[ConnectivityScoreResponse])
def top_mimics(signature_id: UUID, n: int = Query(50, ge=1, le=500)) -> List[ConnectivityScoreResponse]:
    with db_session() as db:
        rows = (
            db.query(ConnectivityScore, LINCSSignature)
            .join(LINCSSignature, LINCSSignature.id == ConnectivityScore.lincs_signature_id)
            .filter(ConnectivityScore.query_signature_id == signature_id)
            .filter(ConnectivityScore.score > 0)
            .order_by(ConnectivityScore.score.desc())
            .limit(n)
            .all()
        )
        out: List[ConnectivityScoreResponse] = []
        for s, lincs_sig in rows:
            out.append(
                ConnectivityScoreResponse(
                    lincs_signature_id=lincs_sig.id,
                    sig_id=lincs_sig.sig_id,
                    pert_iname=lincs_sig.pert_iname,
                    pert_id=lincs_sig.pert_id,
                    cell_id=lincs_sig.cell_id,
                    score=float(s.score),
                )
            )
        return out


@router.get("/compounds/{pert_id}/mechanisms", response_model=List[LINCSSignatureResponse])
def compound_mechanisms(pert_id: str, limit: int = Query(200, ge=1, le=2000)) -> List[LINCSSignatureResponse]:
    """Return LINCS signatures for a given compound pert_id (proxy for mechanism exploration)."""
    with db_session() as db:
        rows = (
            db.query(LINCSSignature)
            .filter(LINCSSignature.pert_id == pert_id)
            .order_by(LINCSSignature.ingested_at.desc())
            .limit(limit)
            .all()
        )
        return [LINCSSignatureResponse.model_validate(r) for r in rows]


__all__ = ["router"]


