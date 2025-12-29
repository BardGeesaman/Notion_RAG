"""
High-dimensional projection API endpoints.

Provides UMAP, t-SNE, and PCA projections for datasets.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

import numpy as np
from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from amprenta_rag.analysis.projection_engine import ProjectorEngine
from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.database.models import Dataset
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

router = APIRouter(prefix="/projector", tags=["Projector"])


class ProjectionRequest(BaseModel):
    """Request for computing projection."""

    dataset_id: Optional[UUID] = Field(None, description="Dataset ID to project")
    matrix: Optional[List[List[float]]] = Field(None, description="Raw matrix to project (n_samples x n_features)")
    algorithm: str = Field(..., description="umap, tsne, or pca")
    n_components: int = Field(2, ge=2, le=3, description="2D or 3D projection")
    n_neighbors: int = Field(15, ge=2, le=100, description="UMAP n_neighbors")
    min_dist: float = Field(0.1, ge=0.0, le=1.0, description="UMAP min_dist")
    perplexity: int = Field(30, ge=5, le=50, description="t-SNE perplexity")
    random_state: int = Field(42, description="Random seed")


class ProjectionResponse(BaseModel):
    """Response for projection computation."""

    coordinates: List[List[float]]
    algorithm_used: str
    n_samples: int
    n_features: int
    cached: bool


class DatasetListResponse(BaseModel):
    """Response for listing available datasets."""

    datasets: List[Dict[str, Any]]
    total: int


class ExportRequest(BaseModel):
    """Request for exporting projection."""

    coordinates: List[List[float]]


class ExportResponse(BaseModel):
    """Response for projection export."""

    csv_data: str
    filename: str


@router.post("/compute", response_model=ProjectionResponse)
async def compute_projection(
    request: ProjectionRequest,
    db: Session = Depends(get_database_session),
) -> ProjectionResponse:
    """
    Compute dimensionality reduction projection.

    Accepts either a dataset_id (fetches data from database) or raw matrix.
    """
    # Validate input
    if not request.dataset_id and not request.matrix:
        raise HTTPException(
            status_code=400,
            detail="Either dataset_id or matrix must be provided"
        )
    
    # Get data matrix
    if request.dataset_id:
        # Fetch dataset and convert to matrix
        dataset = db.query(Dataset).filter(Dataset.id == request.dataset_id).first()
        if not dataset:
            raise HTTPException(status_code=404, detail="Dataset not found")
        
        # Placeholder: would need to load actual dataset matrix
        # For now, create a random matrix
        logger.warning("[PROJECTOR] Dataset loading not implemented - using placeholder")
        np.random.seed(42)
        data_matrix = np.random.randn(100, 20)
    else:
        # Use provided matrix
        data_matrix = np.array(request.matrix)
    
    # Validate matrix
    if data_matrix.size == 0:
        raise HTTPException(status_code=400, detail="Empty data matrix")
    
    if len(data_matrix.shape) != 2:
        raise HTTPException(status_code=400, detail="Matrix must be 2-dimensional")
    
    # Compute projection
    engine = ProjectorEngine()
    
    try:
        if request.algorithm.lower() == "umap":
            result = engine.compute_umap(
                data_matrix,
                n_components=request.n_components,
                n_neighbors=request.n_neighbors,
                min_dist=request.min_dist,
                random_state=request.random_state,
            )
        elif request.algorithm.lower() == "tsne":
            result = engine.compute_tsne(
                data_matrix,
                n_components=request.n_components,
                perplexity=request.perplexity,
                random_state=request.random_state,
            )
        elif request.algorithm.lower() == "pca":
            result = engine.compute_pca(
                data_matrix,
                n_components=request.n_components,
            )
        else:
            raise HTTPException(
                status_code=400,
                detail=f"Unsupported algorithm: {request.algorithm}. Use umap, tsne, or pca."
            )
        
        logger.info(
            "[PROJECTOR] Computed %s projection: %d samples → %dD",
            request.algorithm,
            result.n_samples,
            request.n_components,
        )
        
        return ProjectionResponse(
            coordinates=result.coordinates,
            algorithm_used=result.algorithm_used,
            n_samples=result.n_samples,
            n_features=result.n_features,
            cached=result.cached,
        )
    
    except Exception as e:
        logger.error("[PROJECTOR] Projection failed: %r", e)
        raise HTTPException(status_code=500, detail=f"Projection failed: {str(e)}")


@router.get("/datasets", response_model=DatasetListResponse)
async def list_available_datasets(
    limit: int = 50,
    db: Session = Depends(get_database_session),
) -> DatasetListResponse:
    """
    List datasets available for projection.

    Returns datasets that could be projected (numeric data).
    """
    datasets = db.query(Dataset).order_by(Dataset.created_at.desc()).limit(limit).all()
    
    dataset_list = [
        {
            "id": str(ds.id),
            "name": ds.name,
            "omics_type": ds.omics_type,
            "created_at": ds.created_at.isoformat() if ds.created_at else None,
        }
        for ds in datasets
    ]
    
    return DatasetListResponse(datasets=dataset_list, total=len(dataset_list))


@router.post("/export", response_model=ExportResponse)
async def export_projection(
    request: ExportRequest,
) -> ExportResponse:
    """
    Export projection coordinates as CSV.

    Args:
        request: Export request with coordinates

    Returns:
        CSV data and suggested filename
    """
    coordinates = request.coordinates
    
    if not coordinates:
        raise HTTPException(status_code=400, detail="No coordinates provided")
    
    # Convert to CSV
    import io
    
    output = io.StringIO()
    
    # Determine dimensionality
    n_dims = len(coordinates[0]) if coordinates else 0
    
    if n_dims == 2:
        output.write("Sample,X,Y\n")
        for i, (x, y) in enumerate(coordinates):
            output.write(f"sample_{i},{x},{y}\n")
    elif n_dims == 3:
        output.write("Sample,X,Y,Z\n")
        for i, (x, y, z) in enumerate(coordinates):
            output.write(f"sample_{i},{x},{y},{z}\n")
    else:
        raise HTTPException(status_code=400, detail="Coordinates must be 2D or 3D")
    
    csv_data = output.getvalue()
    
    return ExportResponse(
        csv_data=csv_data,
        filename=f"projection_{n_dims}d.csv",
    )

