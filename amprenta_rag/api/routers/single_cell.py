"""Single-cell omics API endpoints."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field

from amprenta_rag.database.models import (
    CellAnnotation,
    CellCluster,
    CellTypeMarker,
    SingleCellDataset,
)
from amprenta_rag.database.session import db_session
from amprenta_rag.single_cell.h5ad_parser import load_h5ad, validate_h5ad
from amprenta_rag.single_cell.ingest_service import ingest_h5ad, recluster_dataset


router = APIRouter(prefix="/single-cell", tags=["SingleCell"])


class IngestRequest(BaseModel):
    h5ad_path: str
    dataset_id: Optional[UUID] = None


class ReclusterRequest(BaseModel):
    resolution: float = Field(default=0.5, ge=0.05, le=5.0)


class SingleCellDatasetResponse(BaseModel):
    id: UUID
    dataset_id: UUID
    h5ad_path: str
    n_cells: Optional[int] = None
    n_genes: Optional[int] = None
    clustering_resolution: Optional[float] = None
    processing_status: str
    processing_log: Optional[str] = None

    class Config:
        from_attributes = True


class CellClusterResponse(BaseModel):
    cluster_id: int
    cell_type: Optional[str] = None
    n_cells: Optional[int] = None


class CellTypeMarkerResponse(BaseModel):
    gene_symbol: Optional[str] = None
    cluster_id: int
    log2_fold_change: Optional[float] = None
    pval_adj: Optional[float] = None
    feature_id: Optional[UUID] = None


class UMAPPoint(BaseModel):
    barcode: str
    umap_1: Optional[float] = None
    umap_2: Optional[float] = None
    cluster_id: Optional[int] = None


@router.post("/ingest", response_model=SingleCellDatasetResponse)
def ingest(payload: IngestRequest) -> SingleCellDatasetResponse:
    scd = ingest_h5ad(payload.h5ad_path, dataset_id=payload.dataset_id)
    return SingleCellDatasetResponse.model_validate(scd)


@router.get("/datasets", response_model=List[SingleCellDatasetResponse])
def list_single_cell_datasets(limit: int = Query(200, ge=1, le=2000)) -> List[SingleCellDatasetResponse]:
    with db_session() as db:
        rows = db.query(SingleCellDataset).order_by(SingleCellDataset.ingested_at.desc()).limit(limit).all()
        return [SingleCellDatasetResponse.model_validate(r) for r in rows]


@router.get("/datasets/{scd_id}", response_model=SingleCellDatasetResponse)
def get_single_cell_dataset(scd_id: UUID) -> SingleCellDatasetResponse:
    with db_session() as db:
        r = db.query(SingleCellDataset).filter_by(id=scd_id).first()
        if not r:
            raise HTTPException(status_code=404, detail="SingleCellDataset not found")
        return SingleCellDatasetResponse.model_validate(r)


@router.get("/datasets/{scd_id}/clusters", response_model=List[CellClusterResponse])
def get_clusters(scd_id: UUID) -> List[CellClusterResponse]:
    with db_session() as db:
        rows = (
            db.query(CellCluster)
            .filter_by(single_cell_dataset_id=scd_id)
            .order_by(CellCluster.cluster_id.asc())
            .all()
        )
        return [CellClusterResponse(cluster_id=r.cluster_id, cell_type=r.cell_type, n_cells=r.n_cells) for r in rows]


@router.post("/datasets/{scd_id}/recluster")
def recluster(scd_id: UUID, payload: ReclusterRequest) -> dict:
    try:
        recluster_dataset(scd_id, resolution=float(payload.resolution))
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    return {"ok": True, "resolution": float(payload.resolution), "status": "pending"}


@router.get("/datasets/{scd_id}/markers", response_model=List[CellTypeMarkerResponse])
def get_markers(
    scd_id: UUID,
    cluster_id: Optional[int] = Query(None),
    limit: int = Query(200, ge=1, le=5000),
) -> List[CellTypeMarkerResponse]:
    with db_session() as db:
        q = db.query(CellTypeMarker).filter_by(single_cell_dataset_id=scd_id)
        if cluster_id is not None:
            q = q.filter(CellTypeMarker.cluster_id == int(cluster_id))
        rows = q.order_by(CellTypeMarker.pval_adj.asc().nullslast()).limit(limit).all()
        return [
            CellTypeMarkerResponse(
                gene_symbol=r.gene_symbol,
                cluster_id=r.cluster_id,
                log2_fold_change=r.log2_fold_change,
                pval_adj=r.pval_adj,
                feature_id=r.feature_id,
            )
            for r in rows
        ]


@router.get("/datasets/{scd_id}/umap", response_model=List[UMAPPoint])
def get_umap(scd_id: UUID, limit: int = Query(5000, ge=1, le=50000)) -> List[UMAPPoint]:
    with db_session() as db:
        rows = (
            db.query(CellAnnotation)
            .filter_by(single_cell_dataset_id=scd_id)
            .order_by(CellAnnotation.id.asc())
            .limit(limit)
            .all()
        )
        return [
            UMAPPoint(barcode=r.barcode, umap_1=r.umap_1, umap_2=r.umap_2, cluster_id=r.cluster_id)
            for r in rows
        ]


@router.get("/datasets/{scd_id}/expression/{gene_symbol}")
def get_expression(scd_id: UUID, gene_symbol: str) -> dict:
    with db_session() as db:
        scd = db.query(SingleCellDataset).filter_by(id=scd_id).first()
        if not scd:
            raise HTTPException(status_code=404, detail="SingleCellDataset not found")

    try:
        adata = load_h5ad(scd.h5ad_path)
        validate_h5ad(adata)
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))

    gene = str(gene_symbol)
    var_names = [str(x) for x in adata.var_names]
    if gene not in var_names:
        raise HTTPException(status_code=404, detail="Gene not found in var_names")
    idx = var_names.index(gene)
    x = adata.X[:, idx]
    try:
        arr = x.toarray().ravel().tolist()
    except Exception:
        import numpy as np  # type: ignore

        arr = np.asarray(x).ravel().tolist()

    # Return a compact payload; client can join by obs_names
    obs_names = [str(x) for x in adata.obs_names]
    return {"gene_symbol": gene, "obs_names": obs_names, "expression": arr}


__all__ = ["router"]


