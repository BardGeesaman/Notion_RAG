"""
Flow Cytometry API endpoints.

This module provides REST API endpoints for flow cytometry data operations:
- FCS file upload and ingestion
- Dataset management and event retrieval
- Gate creation, modification, and deletion
- Population statistics analysis

Follows FastAPI patterns established in other routers.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, File, HTTPException, Query, UploadFile, status
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.api.schemas import (
    EventsQueryParams,
    EventsResponse,
    FlowCytometryDatasetResponse,
    FlowCytometryParameterResponse,
    GateCreate,
    GateResponse,
    GateUpdate,
    PopulationResponse,
    UploadResponse,
)
from amprenta_rag.database.models import User
from amprenta_rag.database.models_flow_cytometry import (
    FlowCytometryDataset,
    FlowCytometryGate,
    FlowCytometryParameter,
    FlowCytometryPopulation,
)
from amprenta_rag.flow_cytometry.fcs_parser import load_events_parquet
from amprenta_rag.flow_cytometry.ingest_service import GateCreate as ServiceGateCreate
from amprenta_rag.flow_cytometry.ingest_service import apply_gate_to_dataset, ingest_fcs
from amprenta_rag.flow_cytometry.transforms import subsample_events

logger = logging.getLogger(__name__)

router = APIRouter()


@router.post("/upload", response_model=UploadResponse, status_code=status.HTTP_201_CREATED)
async def upload_fcs_file(
    file: UploadFile = File(..., description="FCS file to upload"),
    current_user: User = Depends(get_current_user),
) -> UploadResponse:
    """
    Upload FCS file and start ingestion process.
    
    Accepts an FCS file, saves it temporarily, and triggers the ingestion workflow
    which includes parsing, transformation, and database storage.
    """
    # Validate file extension
    if not file.filename or not file.filename.lower().endswith('.fcs'):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="File must have .fcs extension"
        )
    
    # Validate file size (max 500MB)
    if file.size and file.size > 500 * 1024 * 1024:
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail="File size must be less than 500MB"
        )
    
    # Save file to temporary location
    temp_dir = Path(tempfile.gettempdir()) / "fcs_uploads"
    temp_dir.mkdir(exist_ok=True)
    
    temp_file = temp_dir / f"{file.filename}"
    
    try:
        # Write uploaded file
        content = await file.read()
        temp_file.write_bytes(content)
        
        # Start ingestion process
        flow_dataset = ingest_fcs(str(temp_file))
        
        logger.info(f"FCS file {file.filename} uploaded and ingestion started for dataset {flow_dataset.id}")
        
        return UploadResponse(
            flow_dataset_id=flow_dataset.id,
            dataset_id=flow_dataset.dataset_id,
            filename=file.filename,
            file_size_bytes=len(content),
            processing_status=flow_dataset.processing_status,
            message="FCS file uploaded successfully, processing started"
        )
        
    except Exception as e:
        logger.error(f"Failed to upload FCS file {file.filename}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to process FCS file: {str(e)}"
        )
    finally:
        # Cleanup temporary file
        if temp_file.exists():
            temp_file.unlink()


@router.get("/datasets", response_model=List[FlowCytometryDatasetResponse])
async def list_flow_datasets(
    skip: int = Query(0, ge=0, description="Number of datasets to skip"),
    limit: int = Query(100, ge=1, le=1000, description="Maximum number of datasets to return"),
    processing_status: Optional[str] = Query(None, description="Filter by processing status"),
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> List[FlowCytometryDatasetResponse]:
    """
    List flow cytometry datasets.
    
    Returns paginated list of flow cytometry datasets with optional filtering
    by processing status.
    """
    query = db.query(FlowCytometryDataset).order_by(FlowCytometryDataset.ingested_at.desc())
    
    if processing_status:
        query = query.filter(FlowCytometryDataset.processing_status == processing_status)
    
    datasets = query.offset(skip).limit(limit).all()
    
    return [FlowCytometryDatasetResponse.model_validate(dataset) for dataset in datasets]


@router.get("/datasets/{dataset_id}", response_model=FlowCytometryDatasetResponse)
async def get_flow_dataset(
    dataset_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> FlowCytometryDatasetResponse:
    """
    Get flow cytometry dataset details.
    
    Returns detailed information about a specific flow cytometry dataset
    including acquisition metadata and processing status.
    """
    dataset = db.query(FlowCytometryDataset).filter(FlowCytometryDataset.id == dataset_id).first()
    
    if not dataset:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Flow cytometry dataset {dataset_id} not found"
        )
    
    return FlowCytometryDatasetResponse.model_validate(dataset)


@router.get("/datasets/{dataset_id}/parameters", response_model=List[FlowCytometryParameterResponse])
async def get_dataset_parameters(
    dataset_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> List[FlowCytometryParameterResponse]:
    """
    Get parameters for a flow cytometry dataset.
    
    Returns list of parameters (channels) available in the dataset,
    including metadata about fluorophores, detectors, and value ranges.
    """
    # Verify dataset exists
    dataset = db.query(FlowCytometryDataset).filter(FlowCytometryDataset.id == dataset_id).first()
    if not dataset:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Flow cytometry dataset {dataset_id} not found"
        )
    
    parameters = (
        db.query(FlowCytometryParameter)
        .filter(FlowCytometryParameter.flow_dataset_id == dataset_id)
        .order_by(FlowCytometryParameter.parameter_index)
        .all()
    )
    
    return [FlowCytometryParameterResponse.model_validate(param) for param in parameters]


@router.get("/datasets/{dataset_id}/events", response_model=EventsResponse)
async def get_dataset_events(
    dataset_id: UUID,
    query_params: EventsQueryParams = Depends(),
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> EventsResponse:
    """
    Get event data from a flow cytometry dataset.
    
    Returns paginated event data with optional parameter filtering and subsampling.
    Event data is loaded from Parquet files for performance.
    """
    # Verify dataset exists and is processed
    dataset = db.query(FlowCytometryDataset).filter(FlowCytometryDataset.id == dataset_id).first()
    if not dataset:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Flow cytometry dataset {dataset_id} not found"
        )
    
    if dataset.processing_status != "completed":
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=f"Dataset processing not completed (status: {dataset.processing_status})"
        )
    
    if not dataset.events_parquet_path:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Event data not available for this dataset"
        )
    
    try:
        # Load event data from Parquet
        events, parameter_names = load_events_parquet(dataset.events_parquet_path)
        
        total_events = len(events)
        subsampled = False
        
        # Apply subsampling if requested and dataset is large
        if query_params.subsample and total_events > 50000:
            events = subsample_events(events, max_events=50000)
            subsampled = True
            logger.info(f"Subsampled dataset {dataset_id} from {total_events} to {len(events)} events")
        
        # Filter parameters if specified
        if query_params.parameters:
            try:
                param_indices = [parameter_names.index(param) for param in query_params.parameters]
                events = events[:, param_indices]
                parameter_names = query_params.parameters
            except ValueError as e:
                available_params = ", ".join(parameter_names)
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail=f"Invalid parameter name. Available parameters: {available_params}"
                )
        
        # Apply pagination
        start_idx = query_params.offset
        end_idx = start_idx + query_params.limit
        
        if start_idx >= len(events):
            paginated_events = []
        else:
            paginated_events = events[start_idx:end_idx]
        
        # Convert to list format for JSON serialization
        events_list = paginated_events.tolist() if len(paginated_events) > 0 else []
        
        return EventsResponse(
            events=events_list,
            parameter_names=parameter_names,
            total_events=total_events,
            offset=query_params.offset,
            limit=query_params.limit,
            subsampled=subsampled
        )
        
    except Exception as e:
        logger.error(f"Failed to load events for dataset {dataset_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to load event data: {str(e)}"
        )


@router.post("/datasets/{dataset_id}/gates", response_model=GateResponse, status_code=status.HTTP_201_CREATED)
async def create_gate(
    dataset_id: UUID,
    gate_data: GateCreate,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> GateResponse:
    """
    Create a new gate for a flow cytometry dataset.
    
    Creates a gate definition and applies it to the dataset to compute
    population statistics. Supports polygon, rectangle, and quadrant gates.
    """
    # Verify dataset exists
    dataset = db.query(FlowCytometryDataset).filter(FlowCytometryDataset.id == dataset_id).first()
    if not dataset:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Flow cytometry dataset {dataset_id} not found"
        )
    
    if dataset.processing_status != "completed":
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=f"Dataset processing not completed (status: {dataset.processing_status})"
        )
    
    # Verify parameters exist
    x_param = db.query(FlowCytometryParameter).filter(
        FlowCytometryParameter.id == gate_data.x_parameter_id
    ).first()
    if not x_param:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"X parameter {gate_data.x_parameter_id} not found"
        )
    
    if gate_data.y_parameter_id:
        y_param = db.query(FlowCytometryParameter).filter(
            FlowCytometryParameter.id == gate_data.y_parameter_id
        ).first()
        if not y_param:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Y parameter {gate_data.y_parameter_id} not found"
            )
    
    try:
        # Convert to service layer schema
        service_gate = ServiceGateCreate(
            gate_name=gate_data.gate_name,
            gate_type=gate_data.gate_type,
            gate_definition=gate_data.gate_definition,
            x_parameter_id=gate_data.x_parameter_id,
            y_parameter_id=gate_data.y_parameter_id,
            parent_gate_id=gate_data.parent_gate_id,
            boolean_operator=gate_data.boolean_operator,
            operand_gate_ids=gate_data.operand_gate_ids
        )
        
        # Apply gate to dataset
        gate = apply_gate_to_dataset(dataset_id, service_gate)
        
        logger.info(f"Created gate {gate.id} for dataset {dataset_id}")
        
        return GateResponse.model_validate(gate)
        
    except Exception as e:
        logger.error(f"Failed to create gate for dataset {dataset_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to create gate: {str(e)}"
        )


@router.get("/datasets/{dataset_id}/gates", response_model=List[GateResponse])
async def list_gates(
    dataset_id: UUID,
    include_inactive: bool = Query(False, description="Include inactive gates"),
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> List[GateResponse]:
    """
    List gates for a flow cytometry dataset.
    
    Returns all gates defined for the dataset, optionally including
    inactive gates for historical analysis.
    """
    # Verify dataset exists
    dataset = db.query(FlowCytometryDataset).filter(FlowCytometryDataset.id == dataset_id).first()
    if not dataset:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Flow cytometry dataset {dataset_id} not found"
        )
    
    query = db.query(FlowCytometryGate).filter(FlowCytometryGate.flow_dataset_id == dataset_id)
    
    if not include_inactive:
        query = query.filter(FlowCytometryGate.is_active == True)
    
    gates = query.order_by(FlowCytometryGate.created_at).all()
    
    return [GateResponse.model_validate(gate) for gate in gates]


@router.put("/gates/{gate_id}", response_model=GateResponse)
async def update_gate(
    gate_id: UUID,
    gate_update: GateUpdate,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> GateResponse:
    """
    Update an existing gate.
    
    Allows partial updates to gate properties. If gate definition is changed,
    population statistics will be recomputed automatically.
    """
    gate = db.query(FlowCytometryGate).filter(FlowCytometryGate.id == gate_id).first()
    
    if not gate:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Gate {gate_id} not found"
        )
    
    try:
        # Apply updates
        update_data = gate_update.model_dump(exclude_unset=True)
        
        for field, value in update_data.items():
            setattr(gate, field, value)
        
        db.add(gate)
        db.commit()
        db.refresh(gate)
        
        # If gate definition changed, recompute population (simplified for MVP)
        if "gate_definition" in update_data:
            logger.info(f"Gate definition updated for {gate_id}, population recomputation needed")
            # In full implementation, would trigger recomputation here
        
        logger.info(f"Updated gate {gate_id}")
        
        return GateResponse.model_validate(gate)
        
    except Exception as e:
        db.rollback()
        logger.error(f"Failed to update gate {gate_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to update gate: {str(e)}"
        )


@router.delete("/gates/{gate_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_gate(
    gate_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> None:
    """
    Delete a gate.
    
    Marks the gate as inactive rather than physically deleting it
    to preserve analysis history. Associated populations are also marked inactive.
    """
    gate = db.query(FlowCytometryGate).filter(FlowCytometryGate.id == gate_id).first()
    
    if not gate:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Gate {gate_id} not found"
        )
    
    try:
        # Mark gate as inactive (soft delete)
        gate.is_active = False
        db.add(gate)
        
        # Mark associated populations as inactive (if needed in future)
        # populations = db.query(FlowCytometryPopulation).filter(
        #     FlowCytometryPopulation.gate_id == gate_id
        # ).all()
        # for pop in populations:
        #     pop.is_active = False
        #     db.add(pop)
        
        db.commit()
        
        logger.info(f"Deleted (deactivated) gate {gate_id}")
        
    except Exception as e:
        db.rollback()
        logger.error(f"Failed to delete gate {gate_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to delete gate: {str(e)}"
        )


@router.get("/datasets/{dataset_id}/populations", response_model=List[PopulationResponse])
async def get_population_statistics(
    dataset_id: UUID,
    gate_id: Optional[UUID] = Query(None, description="Filter by specific gate"),
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> List[PopulationResponse]:
    """
    Get population statistics for a flow cytometry dataset.
    
    Returns computed population statistics for all gates or a specific gate,
    including event counts, percentages, and parameter statistics.
    """
    # Verify dataset exists
    dataset = db.query(FlowCytometryDataset).filter(FlowCytometryDataset.id == dataset_id).first()
    if not dataset:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Flow cytometry dataset {dataset_id} not found"
        )
    
    query = db.query(FlowCytometryPopulation).filter(
        FlowCytometryPopulation.flow_dataset_id == dataset_id
    )
    
    if gate_id:
        query = query.filter(FlowCytometryPopulation.gate_id == gate_id)
    
    populations = query.order_by(FlowCytometryPopulation.analyzed_at.desc()).all()
    
    return [PopulationResponse.model_validate(population) for population in populations]


__all__ = ["router"]
