"""
Data export API endpoints.

Provides export endpoints for datasets, experiments, and compounds.
"""

from __future__ import annotations

from typing import List
from uuid import UUID

from fastapi import APIRouter, Depends, Query
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.export.export_engine import (
    create_export_package,
    export_compounds,
    export_dataset,
    export_experiment,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

router = APIRouter(prefix="/export", tags=["Export"])


class ExportCompoundsRequest(BaseModel):
    """Request for compound export."""

    compound_ids: List[UUID]
    format: str = "csv"


class ExportPackageRequest(BaseModel):
    """Request for export package."""

    items: List[str]  # List of export identifiers
    format: str = "zip"


@router.get("/dataset/{dataset_id}")
async def export_dataset_endpoint(
    dataset_id: UUID,
    format: str = Query("csv", description="csv, excel, or json"),
    include_metadata: bool = Query(True, description="Include metadata header"),
    db: Session = Depends(get_database_session),
) -> StreamingResponse:
    """Export dataset to specified format."""
    try:
        data = export_dataset(db, str(dataset_id), format=format, include_metadata=include_metadata)
        
        # Determine content type
        content_types = {
            "csv": "text/csv",
            "excel": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            "json": "application/json",
        }
        
        content_type = content_types.get(format, "application/octet-stream")
        filename = f"dataset_{dataset_id}.{format if format != 'excel' else 'xlsx'}"
        
        return StreamingResponse(
            iter([data]),
            media_type=content_type,
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )
    
    except Exception as e:
        logger.error("[EXPORT] Dataset export failed: %r", e)
        raise


@router.get("/experiment/{experiment_id}")
async def export_experiment_endpoint(
    experiment_id: UUID,
    format: str = Query("csv", description="csv or json"),
    include_datasets: bool = Query(True, description="Include related datasets"),
    db: Session = Depends(get_database_session),
) -> StreamingResponse:
    """Export experiment to specified format."""
    try:
        data = export_experiment(db, str(experiment_id), format=format, include_datasets=include_datasets)
        
        content_type = "text/csv" if format == "csv" else "application/json"
        filename = f"experiment_{experiment_id}.{format}"
        
        return StreamingResponse(
            iter([data]),
            media_type=content_type,
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )
    
    except Exception as e:
        logger.error("[EXPORT] Experiment export failed: %r", e)
        raise


@router.post("/compounds")
async def export_compounds_endpoint(
    request: ExportCompoundsRequest,
    db: Session = Depends(get_database_session),
) -> StreamingResponse:
    """Export compounds list."""
    try:
        compound_ids_str = [str(cid) for cid in request.compound_ids]
        data = export_compounds(db, compound_ids_str, format=request.format)
        
        content_type = "text/csv" if request.format == "csv" else "application/json"
        filename = f"compounds.{request.format}"
        
        return StreamingResponse(
            iter([data]),
            media_type=content_type,
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )
    
    except Exception as e:
        logger.error("[EXPORT] Compounds export failed: %r", e)
        raise


@router.post("/package")
async def create_package_endpoint(
    request: ExportPackageRequest,
    db: Session = Depends(get_database_session),
) -> StreamingResponse:
    """Create export package with multiple files."""
    try:
        # Placeholder: would fetch multiple exports
        exports = {}
        for item_id in request.items:
            # Simplified: just create a placeholder file
            exports[f"{item_id}.csv"] = b"placeholder,data\n1,2"
        
        data = create_export_package(exports, format=request.format)
        
        return StreamingResponse(
            iter([data]),
            media_type="application/zip",
            headers={"Content-Disposition": 'attachment; filename="export_package.zip"'},
        )
    
    except Exception as e:
        logger.error("[EXPORT] Package creation failed: %r", e)
        raise

