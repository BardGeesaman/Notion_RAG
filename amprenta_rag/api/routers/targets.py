"""Target management API endpoints."""

from uuid import UUID
from typing import List, Optional
from fastapi import APIRouter, HTTPException, Query, Depends, status
from sqlalchemy.orm import Session

from amprenta_rag.api import schemas
from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.database.base import get_db
from amprenta_rag.models.auth import User
from amprenta_rag.services import target_service
import logging

logger = logging.getLogger(__name__)

router = APIRouter()


@router.get(
    "",
    summary="List targets",
    description="Retrieve a list of biological targets with optional filtering.",
    response_model=List[schemas.TargetResponse]
)
def list_targets(
    target_class: Optional[str] = Query(None, description="Filter by target class (kinase, gpcr, etc.)"),
    target_family: Optional[str] = Query(None, description="Filter by target family"),
    lifecycle_status: str = Query("active", description="Filter by lifecycle status"),
    skip: int = Query(0, ge=0, description="Number of records to skip"),
    limit: int = Query(100, ge=1, le=500, description="Maximum number of records to return"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """List targets with optional filters.
    
    Returns a paginated list of biological targets. Can be filtered by
    target class, family, and lifecycle status.
    """
    try:
        targets = target_service.list_targets(
            target_class=target_class,
            target_family=target_family,
            lifecycle_status=lifecycle_status,
            skip=skip,
            limit=limit,
            db=db
        )
        return targets
    except Exception as e:
        logger.error(f"Error listing targets: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to retrieve targets"
        )


@router.get(
    "/{target_id}",
    summary="Get target by ID",
    description="Retrieve detailed information about a specific target.",
    response_model=schemas.TargetResponse
)
def get_target(
    target_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Get target by ID.
    
    Returns detailed information about a biological target including
    druggability scores, validation status, and external identifiers.
    """
    target = target_service.get_target(target_id, db=db)
    if not target:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Target {target_id} not found"
        )
    return target


@router.post(
    "",
    summary="Create new target",
    description="Create a new biological target for drug discovery.",
    response_model=schemas.TargetResponse,
    status_code=status.HTTP_201_CREATED
)
def create_target(
    target_data: schemas.TargetCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Create a new target.
    
    Creates a new biological target with the provided information.
    Target names must be unique across the system.
    """
    try:
        # Check if target name already exists
        existing = target_service.get_target_by_name(target_data.name, db=db)
        if existing:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Target with name '{target_data.name}' already exists"
            )
        
        target = target_service.create_target(
            name=target_data.name,
            gene_symbol=target_data.gene_symbol,
            uniprot_id=target_data.uniprot_id,
            target_class=target_data.target_class,
            target_family=target_data.target_family,
            description=target_data.description,
            db=db
        )
        return target
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error creating target: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to create target"
        )


@router.put(
    "/{target_id}",
    summary="Update target",
    description="Update an existing target's properties.",
    response_model=schemas.TargetResponse
)
def update_target(
    target_id: UUID,
    update_data: schemas.TargetUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Update target fields.
    
    Updates the specified fields of an existing target. Only provided
    fields will be updated; others remain unchanged.
    """
    try:
        # Convert Pydantic model to dict, excluding unset values
        updates = update_data.model_dump(exclude_unset=True)
        if not updates:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="No fields provided for update"
            )
        
        target = target_service.update_target(target_id, updates, db=db)
        if not target:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Target {target_id} not found"
            )
        return target
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error updating target {target_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to update target"
        )


@router.get(
    "/{target_id}/assays",
    summary="Get target assays",
    description="Retrieve experiments and assays linked to a target.",
    response_model=List[schemas.TargetAssayResponse]
)
def get_target_assays(
    target_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Get experiments/assays linked to target.
    
    Returns a list of experiments and assays that have compounds
    tested against this target.
    """
    try:
        # Verify target exists
        target = target_service.get_target(target_id, db=db)
        if not target:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Target {target_id} not found"
            )
        
        assays = target_service.get_target_assays(target_id, db=db)
        return assays
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error retrieving assays for target {target_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to retrieve target assays"
        )


@router.get(
    "/{target_id}/compounds",
    summary="Get target compounds",
    description="Retrieve compounds with activity data against a target.",
    response_model=List[schemas.CompoundActivityResponse]
)
def get_target_compounds(
    target_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Get compounds with activity data for target.
    
    Returns compounds that have been tested against this target
    along with their activity measurements (IC50, Ki, etc.).
    """
    try:
        # Verify target exists
        target = target_service.get_target(target_id, db=db)
        if not target:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Target {target_id} not found"
            )
        
        compounds = target_service.get_target_compounds(target_id, db=db)
        return compounds
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error retrieving compounds for target {target_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to retrieve target compounds"
        )


@router.get(
    "/{target_id}/landscape",
    summary="Get competitive landscape",
    description="Retrieve competitive landscape information for a target.",
    response_model=dict
)
def get_target_landscape(
    target_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Get competitive landscape for target.
    
    Returns competitive landscape analysis including published compounds,
    clinical trials, and market information. Currently a stub implementation.
    """
    # Verify target exists
    target = target_service.get_target(target_id, db=db)
    if not target:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Target {target_id} not found"
        )
    
    # Stub implementation - to be enhanced with external data sources
    return {
        "target_id": target_id,
        "target_name": target.name,
        "published_compounds": 0,
        "clinical_trials": 0,
        "approved_drugs": 0,
        "patent_landscape": {},
        "key_players": [],
        "market_size": None,
        "note": "Competitive landscape analysis coming soon"
    }


@router.post(
    "/{target_id}/druggability",
    summary="Calculate druggability",
    description="Calculate druggability score for a target based on various factors.",
    response_model=schemas.DruggabilityResponse
)
def calculate_target_druggability(
    target_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Calculate druggability score for target.
    
    Calculates a druggability score (0-1) based on factors such as:
    - UniProt annotation availability
    - Known binding pockets
    - Validation status
    - Available activity data
    """
    try:
        result = target_service.calculate_druggability(target_id, db=db)
        
        if "error" in result:
            if "not found" in result["error"]:
                raise HTTPException(
                    status_code=status.HTTP_404_NOT_FOUND,
                    detail=f"Target {target_id} not found"
                )
            else:
                raise HTTPException(
                    status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail=result["error"]
                )
        
        return result
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error calculating druggability for target {target_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to calculate druggability"
        )


@router.get(
    "/search/{query}",
    summary="Search targets",
    description="Search targets by name, gene symbol, or description.",
    response_model=List[schemas.TargetResponse]
)
def search_targets(
    query: str,
    limit: int = Query(20, ge=1, le=100, description="Maximum number of results"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Search targets by name, gene symbol, or description.
    
    Performs a text search across target names, gene symbols, and descriptions
    to find matching targets.
    """
    try:
        targets = target_service.search_targets(query, limit=limit, db=db)
        return targets
    except Exception as e:
        logger.error(f"Error searching targets with query '{query}': {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to search targets"
        )
