"""Multi-omics visualization API endpoints."""

from typing import List, Dict, Any
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset
from amprenta_rag.analysis.multi_omics_viz import compute_alluvial_data, compute_upset_data
from amprenta_rag.api.schemas import (
    AlluvialRequest,
    AlluvialResponse,
    UpSetRequest,
    UpSetResponse,
)

router = APIRouter(prefix="/multi-omics-viz", tags=["multi-omics-viz"])


@router.get("/datasets")
async def get_datasets(db: Session = Depends(get_db)) -> List[Dict[str, Any]]:
    """Get all datasets for multi-omics visualization selection."""
    datasets = db.query(Dataset).all()
    return [
        {
            "id": str(d.id),
            "name": d.name,
            "omics_type": d.omics_type,
        }
        for d in datasets
    ]


@router.post("/alluvial", response_model=AlluvialResponse)
async def get_alluvial_data(
    request: AlluvialRequest,
    db: Session = Depends(get_db),
) -> AlluvialResponse:
    """Get alluvial/Sankey diagram data for multi-omics integration."""
    try:
        dataset_ids = [UUID(id_str) for id_str in request.dataset_ids]
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid UUID format")
    data = compute_alluvial_data(dataset_ids, db)
    return AlluvialResponse(**data)


@router.post("/upset", response_model=UpSetResponse)
async def get_upset_data(
    request: UpSetRequest,
    db: Session = Depends(get_db),
) -> UpSetResponse:
    """Get UpSet plot data for feature intersections across datasets."""
    try:
        dataset_ids = [UUID(id_str) for id_str in request.dataset_ids]
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid UUID format")
    data = compute_upset_data(dataset_ids, db)
    return UpSetResponse(**data)