"""
Biophysical Assay API endpoints.

This module provides REST API endpoints for biophysical assay operations:
- SPR, MST, DSC file upload and ingestion
- Experiment management and data retrieval
- Kinetic/affinity/thermal analysis refitting
- Cross-assay comparison and profiling

Follows FastAPI patterns established in other routers.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, File, HTTPException, Query, UploadFile, status
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.api.async_dependencies import get_async_database_session
from amprenta_rag.api.schemas import (
    BiophysicalCompareResponse,
    BiophysicalUploadResponse,
    DSCExperimentResponse,
    MSTExperimentResponse,
    RefitRequest,
    SPRExperimentResponse,
)
from amprenta_rag.biophysical.ingest_service import (
    get_compound_biophysical_profile,
    ingest_dsc_file,
    ingest_mst_file,
    ingest_spr_file,
    reprocess_experiment,
)
from amprenta_rag.database.models import User
from amprenta_rag.database.models_biophysical import (
    DSCExperiment,
    MSTExperiment,
    SPRExperiment,
)

logger = logging.getLogger(__name__)

router = APIRouter()


# ============================================================================
# SPR Endpoints
# ============================================================================


@router.post("/spr/upload", response_model=BiophysicalUploadResponse, status_code=status.HTTP_201_CREATED)
async def upload_spr_file(
    file: UploadFile = File(...),
    experiment_id: Optional[UUID] = None,
    compound_id: Optional[UUID] = None,
    target_id: Optional[UUID] = None,
    target_name: Optional[str] = None,
    db: AsyncSession = Depends(get_async_database_session),
    current_user: User = Depends(get_current_user),
):
    """Upload and process SPR data file (CSV, TXT)."""
    
    if not file.filename:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Filename is required"
        )
    
    # Validate file extension
    allowed_extensions = ['.csv', '.txt']
    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in allowed_extensions:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid file type. Allowed: {', '.join(allowed_extensions)}"
        )
    
    try:
        # Save uploaded file to temporary location
        with tempfile.NamedTemporaryFile(delete=False, suffix=file_ext) as tmp:
            content = await file.read()
            tmp.write(content)
            tmp_path = tmp.name
        
        # Call ingest service
        experiment = ingest_spr_file(
            file_path=tmp_path,
            experiment_id=experiment_id,
            compound_id=compound_id,
            target_id=target_id,
            target_name=target_name,
            user_id=current_user.id
        )
        
        # Clean up temp file
        Path(tmp_path).unlink()
        
        return BiophysicalUploadResponse(
            experiment_id=experiment.id,
            filename=file.filename,
            file_size_bytes=len(content),
            processing_status=experiment.processing_status,
            message=f"SPR experiment {experiment.experiment_name or experiment.id} created successfully"
        )
        
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid SPR file: {str(e)}"
        )
    except Exception as e:
        logger.exception(f"Failed to process SPR file {file.filename}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to process SPR file: {str(e)}"
        )


@router.get("/spr", response_model=List[SPRExperimentResponse])
async def list_spr_experiments(
    compound_id: Optional[UUID] = Query(None, description="Filter by compound ID"),
    target_id: Optional[UUID] = Query(None, description="Filter by target ID"),
    skip: int = Query(0, ge=0, description="Number of experiments to skip"),
    limit: int = Query(100, ge=1, le=1000, description="Maximum number of experiments to return"),
    db: AsyncSession = Depends(get_async_database_session),
):
    """List SPR experiments with optional filtering."""
    
    try:
        query = select(SPRExperiment)
        
        # Apply filters
        if compound_id:
            query = query.filter(SPRExperiment.compound_id == compound_id)
        if target_id:
            query = query.filter(SPRExperiment.target_id == target_id)
        
        # Apply pagination
        paginated_query = query.order_by(SPRExperiment.created_at.desc()).offset(skip).limit(limit)
        result = await db.execute(paginated_query)
        experiments = result.scalars().all()
        
        # Detach from session to prevent DetachedInstanceError
        for exp in experiments:
            db.expunge(exp)
        
        return experiments
        
    except Exception as e:
        logger.exception("Failed to list SPR experiments")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve SPR experiments: {str(e)}"
        )


@router.get("/spr/{spr_id}", response_model=SPRExperimentResponse)
async def get_spr_experiment(
    spr_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
):
    """Get SPR experiment details by ID."""
    
    try:
        result = await db.execute(
            select(SPRExperiment).filter(SPRExperiment.id == spr_id)
        )
        experiment = result.scalars().first()
        
        if not experiment:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"SPR experiment {spr_id} not found"
            )
        
        # Detach from session
        db.expunge(experiment)
        return experiment
        
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Failed to get SPR experiment {spr_id}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve SPR experiment: {str(e)}"
        )


@router.post("/spr/{spr_id}/fit", response_model=SPRExperimentResponse)
async def refit_spr_kinetics(
    spr_id: UUID,
    request: RefitRequest,
    db: AsyncSession = Depends(get_async_database_session),
):
    """Refit SPR kinetics with specified model."""
    
    try:
        # Check if experiment exists
        result = await db.execute(
            select(SPRExperiment).filter(SPRExperiment.id == spr_id)
        )
        experiment = result.scalars().first()
        
        if not experiment:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"SPR experiment {spr_id} not found"
            )
        
        # Call reprocess service
        updated_experiment = reprocess_experiment(
            experiment_id=spr_id,
            assay_type="spr",
            model=request.model
        )
        
        if not updated_experiment:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Failed to reprocess experiment"
            )
        
        # Detach from session
        db.expunge(updated_experiment)
        return updated_experiment
        
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Failed to refit SPR experiment {spr_id}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to refit SPR kinetics: {str(e)}"
        )


# ============================================================================
# MST Endpoints
# ============================================================================


@router.post("/mst/upload", response_model=BiophysicalUploadResponse, status_code=status.HTTP_201_CREATED)
async def upload_mst_file(
    file: UploadFile = File(...),
    experiment_id: Optional[UUID] = None,
    compound_id: Optional[UUID] = None,
    target_id: Optional[UUID] = None,
    target_name: Optional[str] = None,
    db: AsyncSession = Depends(get_async_database_session),
    current_user: User = Depends(get_current_user),
):
    """Upload and process MST data file (CSV, XLSX)."""
    
    if not file.filename:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Filename is required"
        )
    
    # Validate file extension
    allowed_extensions = ['.csv', '.xlsx']
    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in allowed_extensions:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid file type. Allowed: {', '.join(allowed_extensions)}"
        )
    
    try:
        # Save uploaded file to temporary location
        with tempfile.NamedTemporaryFile(delete=False, suffix=file_ext) as tmp:
            content = await file.read()
            tmp.write(content)
            tmp_path = tmp.name
        
        # Call ingest service
        experiment = ingest_mst_file(
            file_path=tmp_path,
            experiment_id=experiment_id,
            compound_id=compound_id,
            target_id=target_id,
            target_name=target_name,
            user_id=current_user.id
        )
        
        # Clean up temp file
        Path(tmp_path).unlink()
        
        return BiophysicalUploadResponse(
            experiment_id=experiment.id,
            filename=file.filename,
            file_size_bytes=len(content),
            processing_status=experiment.processing_status,
            message=f"MST experiment {experiment.experiment_name or experiment.id} created successfully"
        )
        
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid MST file: {str(e)}"
        )
    except Exception as e:
        logger.exception(f"Failed to process MST file {file.filename}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to process MST file: {str(e)}"
        )


@router.get("/mst", response_model=List[MSTExperimentResponse])
async def list_mst_experiments(
    compound_id: Optional[UUID] = Query(None, description="Filter by compound ID"),
    target_id: Optional[UUID] = Query(None, description="Filter by target ID"),
    skip: int = Query(0, ge=0, description="Number of experiments to skip"),
    limit: int = Query(100, ge=1, le=1000, description="Maximum number of experiments to return"),
    db: AsyncSession = Depends(get_async_database_session),
):
    """List MST experiments with optional filtering."""
    
    try:
        query = select(MSTExperiment)
        
        # Apply filters
        if compound_id:
            query = query.filter(MSTExperiment.compound_id == compound_id)
        if target_id:
            query = query.filter(MSTExperiment.target_id == target_id)
        
        # Apply pagination
        paginated_query = query.order_by(MSTExperiment.created_at.desc()).offset(skip).limit(limit)
        result = await db.execute(paginated_query)
        experiments = result.scalars().all()
        
        # Detach from session to prevent DetachedInstanceError
        for exp in experiments:
            db.expunge(exp)
        
        return experiments
        
    except Exception as e:
        logger.exception("Failed to list MST experiments")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve MST experiments: {str(e)}"
        )


@router.get("/mst/{mst_id}", response_model=MSTExperimentResponse)
async def get_mst_experiment(
    mst_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
):
    """Get MST experiment details by ID."""
    
    try:
        result = await db.execute(
            select(MSTExperiment).filter(MSTExperiment.id == mst_id)
        )
        experiment = result.scalars().first()
        
        if not experiment:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"MST experiment {mst_id} not found"
            )
        
        # Detach from session
        db.expunge(experiment)
        return experiment
        
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Failed to get MST experiment {mst_id}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve MST experiment: {str(e)}"
        )


@router.post("/mst/{mst_id}/fit", response_model=MSTExperimentResponse)
async def refit_mst_affinity(
    mst_id: UUID,
    request: RefitRequest,
    db: AsyncSession = Depends(get_async_database_session),
):
    """Refit MST affinity with specified model."""
    
    try:
        # Check if experiment exists
        result = await db.execute(
            select(MSTExperiment).filter(MSTExperiment.id == mst_id)
        )
        experiment = result.scalars().first()
        
        if not experiment:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"MST experiment {mst_id} not found"
            )
        
        # Call reprocess service
        updated_experiment = reprocess_experiment(
            experiment_id=mst_id,
            assay_type="mst",
            model=request.model
        )
        
        if not updated_experiment:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Failed to reprocess experiment"
            )
        
        # Detach from session
        db.expunge(updated_experiment)
        return updated_experiment
        
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Failed to refit MST experiment {mst_id}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to refit MST affinity: {str(e)}"
        )


# ============================================================================
# DSC Endpoints
# ============================================================================


@router.post("/dsc/upload", response_model=BiophysicalUploadResponse, status_code=status.HTTP_201_CREATED)
async def upload_dsc_file(
    file: UploadFile = File(...),
    experiment_id: Optional[UUID] = None,
    compound_id: Optional[UUID] = None,
    protein_id: Optional[UUID] = None,
    protein_name: Optional[str] = None,
    db: AsyncSession = Depends(get_async_database_session),
    current_user: User = Depends(get_current_user),
):
    """Upload and process DSC data file (CSV, TXT)."""
    
    if not file.filename:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Filename is required"
        )
    
    # Validate file extension
    allowed_extensions = ['.csv', '.txt']
    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in allowed_extensions:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid file type. Allowed: {', '.join(allowed_extensions)}"
        )
    
    try:
        # Save uploaded file to temporary location
        with tempfile.NamedTemporaryFile(delete=False, suffix=file_ext) as tmp:
            content = await file.read()
            tmp.write(content)
            tmp_path = tmp.name
        
        # Call ingest service (map protein_id to target_id for DSC)
        experiment = ingest_dsc_file(
            file_path=tmp_path,
            experiment_id=experiment_id,
            compound_id=compound_id,
            target_id=protein_id,
            target_name=protein_name,
            user_id=current_user.id
        )
        
        # Clean up temp file
        Path(tmp_path).unlink()
        
        return BiophysicalUploadResponse(
            experiment_id=experiment.id,
            filename=file.filename,
            file_size_bytes=len(content),
            processing_status=experiment.processing_status,
            message=f"DSC experiment {experiment.experiment_name or experiment.id} created successfully"
        )
        
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid DSC file: {str(e)}"
        )
    except Exception as e:
        logger.exception(f"Failed to process DSC file {file.filename}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to process DSC file: {str(e)}"
        )


@router.get("/dsc", response_model=List[DSCExperimentResponse])
async def list_dsc_experiments(
    compound_id: Optional[UUID] = Query(None, description="Filter by compound ID"),
    protein_id: Optional[UUID] = Query(None, description="Filter by protein ID"),
    skip: int = Query(0, ge=0, description="Number of experiments to skip"),
    limit: int = Query(100, ge=1, le=1000, description="Maximum number of experiments to return"),
    db: AsyncSession = Depends(get_async_database_session),
):
    """List DSC experiments with optional filtering."""
    
    try:
        query = select(DSCExperiment)
        
        # Apply filters (map protein_id to target_id for DSC)
        if compound_id:
            query = query.filter(DSCExperiment.compound_id == compound_id)
        if protein_id:
            query = query.filter(DSCExperiment.target_id == protein_id)
        
        # Apply pagination
        paginated_query = query.order_by(DSCExperiment.created_at.desc()).offset(skip).limit(limit)
        result = await db.execute(paginated_query)
        experiments = result.scalars().all()
        
        # Detach from session to prevent DetachedInstanceError
        for exp in experiments:
            db.expunge(exp)
        
        return experiments
        
    except Exception as e:
        logger.exception("Failed to list DSC experiments")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve DSC experiments: {str(e)}"
        )


@router.get("/dsc/{dsc_id}", response_model=DSCExperimentResponse)
async def get_dsc_experiment(
    dsc_id: UUID,
    db: AsyncSession = Depends(get_async_database_session),
):
    """Get DSC experiment details by ID."""
    
    try:
        result = await db.execute(
            select(DSCExperiment).filter(DSCExperiment.id == dsc_id)
        )
        experiment = result.scalars().first()
        
        if not experiment:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"DSC experiment {dsc_id} not found"
            )
        
        # Detach from session
        db.expunge(experiment)
        return experiment
        
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Failed to get DSC experiment {dsc_id}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve DSC experiment: {str(e)}"
        )


@router.post("/dsc/{dsc_id}/fit", response_model=DSCExperimentResponse)
async def refit_dsc_thermal(
    dsc_id: UUID,
    request: RefitRequest,
    db: AsyncSession = Depends(get_async_database_session),
):
    """Refit DSC thermal analysis with specified model."""
    
    try:
        # Check if experiment exists
        result = await db.execute(
            select(DSCExperiment).filter(DSCExperiment.id == dsc_id)
        )
        experiment = result.scalars().first()
        
        if not experiment:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"DSC experiment {dsc_id} not found"
            )
        
        # Call reprocess service
        updated_experiment = reprocess_experiment(
            experiment_id=dsc_id,
            assay_type="dsc",
            model=request.model
        )
        
        if not updated_experiment:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Failed to reprocess experiment"
            )
        
        # Detach from session
        db.expunge(updated_experiment)
        return updated_experiment
        
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Failed to refit DSC experiment {dsc_id}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to refit DSC thermal analysis: {str(e)}"
        )


# ============================================================================
# Cross-Assay Endpoint
# ============================================================================


@router.get("/compare", response_model=BiophysicalCompareResponse)
async def compare_biophysical_results(
    compound_id: UUID = Query(..., description="Compound ID to compare across assays"),
    db: AsyncSession = Depends(get_async_database_session),
):
    """Compare biophysical results across SPR, MST, and DSC assays for a compound."""
    
    try:
        # Get compound biophysical profile
        profile = get_compound_biophysical_profile(compound_id)
        
        if not profile:
            # Return empty results if no data found
            return BiophysicalCompareResponse(
                compound_id=compound_id,
                spr_results=[],
                mst_results=[],
                dsc_results=[]
            )
        
        # Extract results by assay type
        spr_results = [exp for exp in profile.get('spr_experiments', [])]
        mst_results = [exp for exp in profile.get('mst_experiments', [])]
        dsc_results = [exp for exp in profile.get('dsc_experiments', [])]
        
        # Detach all objects from session
        for exp in spr_results:
            db.expunge(exp)
        for exp in mst_results:
            db.expunge(exp)
        for exp in dsc_results:
            db.expunge(exp)
        
        return BiophysicalCompareResponse(
            compound_id=compound_id,
            spr_results=spr_results,
            mst_results=mst_results,
            dsc_results=dsc_results
        )
        
    except Exception as e:
        logger.exception(f"Failed to compare biophysical results for compound {compound_id}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve biophysical comparison: {str(e)}"
        )
