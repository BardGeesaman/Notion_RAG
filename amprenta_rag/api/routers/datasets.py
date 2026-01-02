"""
API router for Datasets.
"""

import json
import re
from datetime import datetime
from typing import Any, Dict, List, Optional, cast
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, Response
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select

from amprenta_rag.api.async_dependencies import get_async_database_session
from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.database.models import User
from amprenta_rag.api.schemas import (
    AnnotationCreate,
    Dataset,
    DatasetCreate,
    DatasetFinderRequest,
    DatasetFinderResponse,
    DatasetUpdate,
    EnrichmentStatusResponse,
    MetadataEnrichmentResponse,
)
from amprenta_rag.api.services import datasets as dataset_service
from amprenta_rag.database.models import Dataset as DatasetModel, Note
from amprenta_rag.models.domain import OmicsType
from amprenta_rag.notebooks import generate_dataset_notebook
from amprenta_rag.utils.optimistic_lock import ConflictError, update_with_lock
from amprenta_rag.utils.uuid_utils import ensure_uuid

router = APIRouter()


@router.post("/", response_model=Dataset, status_code=201)
async def create_dataset(
    dataset: DatasetCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session),
) -> Dataset:
    """Create a new dataset."""
    return cast(Dataset, dataset_service.create_dataset(db, dataset))


@router.get("/", response_model=List[Dataset])
async def list_datasets(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    omics_type: Optional[OmicsType] = Query(None, description="Filter by omics type"),
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    experiment_id: Optional[UUID] = Query(None, description="Filter by experiment ID"),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session),
) -> List[Dataset]:
    """List all datasets."""
    omics_type_str = omics_type.value if omics_type else None
    return cast(
        List[Dataset],
        dataset_service.get_datasets(
            db,
            skip=skip,
            limit=limit,
            name_filter=name,
            omics_type=omics_type_str,
            program_id=program_id,
            experiment_id=experiment_id,
        ),
    )


@router.get("/{dataset_id}", response_model=Dataset)
async def get_dataset(
    dataset_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
) -> Dataset:
    """Get a dataset by ID."""
    dataset = dataset_service.get_dataset(db, dataset_id)
    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")
    return cast(Dataset, dataset)


@router.post("/{dataset_id}/annotations", summary="Add annotation to dataset")
async def add_dataset_annotation(
    dataset_id: UUID,
    annotation: AnnotationCreate,
    db: AsyncSession = Depends(get_async_database_session),
) -> Dict[str, Any]:
    """Add a note/annotation to a dataset."""
    dataset = dataset_service.get_dataset(db, dataset_id)
    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    note = Note(
        entity_type="dataset",
        entity_id=cast(Any, ensure_uuid(dataset_id)),
        annotation_type=annotation.annotation_type,
        content=annotation.text,
    )
    db.add(note)
    await db.commit()
    await db.refresh(note)

    created_at_val = getattr(note, "created_at", None)
    return {
        "id": str(note.id),
        "entity_type": note.entity_type,
        "entity_id": str(note.entity_id),
        "text": note.content,
        "annotation_type": note.annotation_type,
        "created_at": created_at_val.isoformat() if isinstance(created_at_val, datetime) else None,
    }


@router.get("/{dataset_id}/notebook")
async def download_dataset_notebook(
    dataset_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
) -> Response:
    """
    Download an nbformat v4 notebook (JSON) for exploring a dataset.
    """
    dataset = dataset_service.get_dataset(db, dataset_id)
    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    nb = generate_dataset_notebook(str(dataset_id))

    dataset_name = getattr(dataset, "name", None) or f"dataset_{str(dataset_id)[:8]}"
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", dataset_name).strip("_") or f"dataset_{str(dataset_id)[:8]}"

    headers = {"Content-Disposition": f'attachment; filename="{safe_name}.ipynb"'}
    return Response(content=json.dumps(nb), media_type="application/json", headers=headers)


@router.patch("/{dataset_id}", response_model=Dataset)
async def update_dataset(
    dataset_id: UUID,
    dataset: DatasetUpdate,
    db: AsyncSession = Depends(get_async_database_session),
) -> Dataset:
    """
    Update a dataset using optimistic locking.

    Clients must provide the current `version` of the dataset. If the stored
    version differs, the request fails with HTTP 409 to prevent lost updates.
    """
    db_dataset = dataset_service.get_dataset(db, dataset_id)
    if not db_dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    update_data = dataset.model_dump(exclude_unset=True)
    expected_version = cast(int, update_data.pop("version"))
    program_ids = update_data.pop("program_ids", None)
    experiment_ids = update_data.pop("experiment_ids", None)

    # Handle omics_type enum conversion
    if "omics_type" in update_data and update_data["omics_type"]:
        update_data["omics_type"] = update_data["omics_type"].value

    # Update program relationships if provided
    if program_ids is not None:
        from amprenta_rag.api.services.programs import get_program

        db_dataset.programs.clear()
        for program_id in program_ids:
            program = get_program(db, program_id)
            if program:
                db_dataset.programs.append(program)

    # Update experiment relationships if provided
    if experiment_ids is not None:
        from amprenta_rag.api.services.experiments import get_experiment

        db_dataset.experiments.clear()
        for experiment_id in experiment_ids:
            experiment = get_experiment(db, experiment_id)
            if experiment:
                db_dataset.experiments.append(experiment)

    try:
        updated = update_with_lock(db_dataset, cast(Dict[str, Any], update_data), expected_version, db)
    except ConflictError:
        raise HTTPException(status_code=409, detail="Version conflict")

    db.refresh(updated)
    return cast(Dataset, updated)


@router.delete("/{dataset_id}", status_code=204)
async def delete_dataset(
    dataset_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
) -> None:
    """Delete a dataset."""
    success = dataset_service.delete_dataset(db, dataset_id)
    if not success:
        raise HTTPException(status_code=404, detail="Dataset not found")


@router.post("/find", response_model=DatasetFinderResponse)
async def find_datasets(
    request: DatasetFinderRequest,
) -> DatasetFinderResponse:
    """
    Find datasets across repositories using natural language query.
    
    Uses AI to extract search terms from natural language and searches
    across GEO, ArrayExpress, and Metabolomics Workbench databases.
    """
    try:
        from amprenta_rag.query.dataset_finder import find_datasets_by_nl
        
        result = find_datasets_by_nl(
            query=request.query,
            repositories=request.repositories,
            max_results=request.max_results,
        )
        
        # Convert to response schema
        return DatasetFinderResponse(
            query=result.query,
            extracted_terms=result.extracted_terms,
            results=[
                {
                    "accession": r.accession,
                    "title": r.title,
                    "description": r.description,
                    "source": r.source,
                    "species": r.species,
                    "tissue": r.tissue,
                    "disease": r.disease,
                    "assay_type": r.assay_type,
                    "sample_count": r.sample_count,
                    "url": r.url,
                    "score": r.score,
                }
                for r in result.results
            ],
            total_found=result.total_found,
            sources_searched=result.sources_searched,
            sources_failed=result.sources_failed,
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Dataset search failed: {str(e)}"
        )


@router.post("/{dataset_id}/enrich", response_model=MetadataEnrichmentResponse)
async def enrich_dataset_metadata(
    dataset_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
) -> MetadataEnrichmentResponse:
    """
    Enrich dataset metadata using LLM semantic extraction.
    
    Extracts structured metadata from dataset description/abstract using AI
    and merges it with existing metadata. LLM-extracted values take precedence.
    """
    try:
        from amprenta_rag.services.metadata_enrichment import enrich_dataset_metadata as enrich_func
        
        # Check if dataset exists first
        result = await db.execute(
            select(DatasetModel).filter(DatasetModel.id == dataset_id)
        )
        dataset = result.scalars().first()
        if not dataset:
            raise HTTPException(status_code=404, detail="Dataset not found")
        
        # Perform enrichment
        result = enrich_func(dataset_id)
        
        return MetadataEnrichmentResponse(
            dataset_id=result.dataset_id,
            success=result.success,
            enriched_fields=result.enriched_fields,
            extracted_metadata=result.extracted_metadata,
            error_message=result.error_message,
            processing_time_seconds=result.processing_time_seconds,
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Metadata enrichment failed: {str(e)}"
        )


@router.get("/{dataset_id}/enrichment-status", response_model=EnrichmentStatusResponse)
async def get_dataset_enrichment_status(
    dataset_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
) -> EnrichmentStatusResponse:
    """
    Get enrichment status for a dataset.
    
    Returns information about whether the dataset has been enriched,
    when it was enriched, and what fields are available.
    """
    try:
        from amprenta_rag.services.metadata_enrichment import get_enrichment_status
        
        # Check if dataset exists first
        result = await db.execute(
            select(DatasetModel).filter(DatasetModel.id == dataset_id)
        )
        dataset = result.scalars().first()
        if not dataset:
            raise HTTPException(status_code=404, detail="Dataset not found")
        
        status = get_enrichment_status(dataset_id)
        
        if "error" in status:
            return EnrichmentStatusResponse(
                enriched=False,
                error=status["error"]
            )
        
        return EnrichmentStatusResponse(
            enriched=status["enriched"],
            enrichment_timestamp=status["enrichment_timestamp"],
            enrichment_version=status["enrichment_version"],
            available_fields=status["available_fields"],
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to get enrichment status: {str(e)}"
        )


# --- Expression Overlay Schemas ---

class ExpressionOverlayRequest(BaseModel):
    """Request for expression overlay data."""
    gene_symbols: List[str]
    colormap: str = "diverging"


class ExpressionOverlayResponse(BaseModel):
    """Response for expression overlay data."""
    node_colors: Dict[str, str]
    expression_values: Dict[str, float]
    genes_found: int
    genes_missing: int


# --- Expression Overlay Endpoints ---

@router.post("/{dataset_id}/expression-overlay", response_model=ExpressionOverlayResponse)
async def get_dataset_expression_overlay(
    dataset_id: UUID,
    request: ExpressionOverlayRequest,
) -> ExpressionOverlayResponse:
    """Get expression values and colors for network overlay."""
    from amprenta_rag.services.expression_overlay import get_expression_overlay
    
    colors, expression = get_expression_overlay(
        request.gene_symbols,
        dataset_id,
        request.colormap,
    )
    
    return ExpressionOverlayResponse(
        node_colors=colors,
        expression_values=expression,
        genes_found=len(expression),
        genes_missing=len(request.gene_symbols) - len(expression),
    )

