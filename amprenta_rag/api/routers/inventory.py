"""Compound inventory API endpoints."""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, status
from sqlalchemy.orm import Session

from amprenta_rag.api import schemas
from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.models.auth import User
from amprenta_rag.services import inventory as service
import logging

logger = logging.getLogger(__name__)

router = APIRouter()


# ============================================================================
# COMPOUND SAMPLE ENDPOINTS (8)
# ============================================================================

@router.post(
    "/samples", 
    response_model=schemas.CompoundSampleResponse, 
    status_code=status.HTTP_201_CREATED,
    summary="Create compound sample",
    description="Register a new compound sample in inventory."
)
def create_compound_sample(
    data: schemas.CompoundSampleCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Register a new compound sample.
    
    Creates a new compound sample with barcode generation, storage tracking,
    and concentration/expiry management for inventory control.
    """
    try:
        sample = service.create_compound_sample(
            db=db,
            compound_id=data.compound_id,
            quantity=data.quantity,
            quantity_unit=data.quantity_unit,
            concentration=data.concentration,
            concentration_unit=data.concentration_unit,
            solvent=data.solvent,
            format=data.format,
            batch_lot=data.batch_lot,
            expiry_date=data.expiry_date,
            storage_location_id=data.storage_location_id,
            position=data.position,
            plate_id=data.plate_id,
            well_position=data.well_position,
            barcode=data.barcode,
            created_by_id=current_user.id,
            notes=data.notes,
        )
        return schemas.CompoundSampleResponse.model_validate(sample)
    except ValueError as e:
        logger.error(f"Failed to create compound sample: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.get(
    "/samples", 
    response_model=List[schemas.CompoundSampleResponse],
    summary="List compound samples",
    description="List compound samples with optional filtering."
)
def list_compound_samples(
    compound_id: Optional[UUID] = Query(None, description="Filter by compound"),
    status: Optional[str] = Query(None, description="Filter by sample status"),
    storage_location_id: Optional[UUID] = Query(None, description="Filter by storage location"),
    plate_id: Optional[UUID] = Query(None, description="Filter by plate"),
    limit: int = Query(200, ge=1, le=1000, description="Maximum results"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """List compound samples with optional filters.
    
    Returns paginated list of compound samples with filtering options
    for compound, status, storage location, and plate.
    """
    samples = service.list_compound_samples(
        db=db,
        compound_id=compound_id,
        status=status,
        storage_location_id=storage_location_id,
        plate_id=plate_id,
        limit=limit,
    )
    return [schemas.CompoundSampleResponse.model_validate(s) for s in samples]


@router.get(
    "/samples/{sample_id}", 
    response_model=schemas.CompoundSampleResponse,
    summary="Get compound sample",
    description="Get detailed information about a specific compound sample."
)
def get_compound_sample(
    sample_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get compound sample by ID.
    
    Returns detailed information about a compound sample including
    storage location, concentration, and expiry information.
    """
    sample = service.get_compound_sample(db, sample_id)
    if not sample:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"Sample {sample_id} not found"
        )
    return schemas.CompoundSampleResponse.model_validate(sample)


@router.get(
    "/samples/barcode/{barcode}", 
    response_model=schemas.CompoundSampleResponse,
    summary="Lookup sample by barcode",
    description="Find compound sample using barcode scan."
)
def get_sample_by_barcode(
    barcode: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Lookup compound sample by barcode.
    
    Enables barcode scanning for quick sample identification
    and retrieval in laboratory workflows.
    """
    sample = service.get_compound_sample_by_barcode(db, barcode)
    if not sample:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"Sample with barcode '{barcode}' not found"
        )
    return schemas.CompoundSampleResponse.model_validate(sample)


@router.put(
    "/samples/{sample_id}", 
    response_model=schemas.CompoundSampleResponse,
    summary="Update compound sample",
    description="Update sample properties like quantity and storage location."
)
def update_compound_sample(
    sample_id: UUID,
    data: schemas.CompoundSampleUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Update compound sample.
    
    Updates sample properties with automatic depletion logic
    when quantity reaches zero.
    """
    try:
        sample = service.update_compound_sample(
            db=db,
            sample_id=sample_id,
            quantity=data.quantity,
            status=data.status,
            storage_location_id=data.storage_location_id,
            position=data.position,
            notes=data.notes,
        )
        return schemas.CompoundSampleResponse.model_validate(sample)
    except ValueError as e:
        logger.error(f"Failed to update sample {sample_id}: {e}")
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))


@router.post(
    "/samples/{sample_id}/transfer", 
    response_model=dict,
    summary="Transfer sample",
    description="Transfer sample to new storage location with audit trail."
)
def transfer_compound_sample(
    sample_id: UUID,
    data: schemas.CompoundSampleTransfer,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Transfer sample to new storage location.
    
    Creates audit trail record using SampleTransfer model
    and updates sample location.
    """
    try:
        transfer = service.transfer_compound_sample(
            db=db,
            sample_id=sample_id,
            to_location_id=data.to_location_id,
            transferred_by_id=current_user.id,
            new_position=data.new_position,
            notes=data.notes,
        )
        return {"transfer_id": str(transfer.id), "message": "Transfer recorded"}
    except ValueError as e:
        logger.error(f"Failed to transfer sample {sample_id}: {e}")
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))


@router.post(
    "/samples/{sample_id}/dispense", 
    response_model=schemas.CompoundSampleResponse,
    summary="Dispense from sample",
    description="Record dispensing from sample and deduct quantity."
)
def dispense_from_sample(
    sample_id: UUID,
    data: schemas.CompoundSampleDispense,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Record dispensing from a sample (deduct quantity).
    
    Validates sufficient quantity and automatically marks
    sample as depleted when quantity reaches zero.
    """
    sample = service.get_compound_sample(db, sample_id)
    if not sample:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail="Sample not found"
        )
    
    if sample.quantity < data.quantity:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Insufficient quantity: {sample.quantity} available, {data.quantity} requested"
        )
    
    new_quantity = sample.quantity - data.quantity
    updated = service.update_compound_sample(db, sample_id, quantity=new_quantity, notes=data.notes)
    return schemas.CompoundSampleResponse.model_validate(updated)


@router.delete(
    "/samples/{sample_id}", 
    status_code=status.HTTP_204_NO_CONTENT,
    summary="Archive sample",
    description="Soft delete sample by setting status to archived."
)
def delete_compound_sample(
    sample_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Soft delete sample (set status=archived).
    
    Archives sample rather than hard deletion for audit trail.
    Sample can be recovered by updating status back to active.
    """
    try:
        service.update_compound_sample(db, sample_id, status="archived")
    except ValueError as e:
        logger.error(f"Failed to archive sample {sample_id}: {e}")
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))


# ============================================================================
# COMPOUND PLATE ENDPOINTS (4)
# ============================================================================

@router.post(
    "/plates", 
    response_model=schemas.CompoundPlateResponse, 
    status_code=status.HTTP_201_CREATED,
    summary="Create compound plate",
    description="Create a new compound storage plate."
)
def create_compound_plate(
    data: schemas.CompoundPlateCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Create a new compound plate.
    
    Creates compound storage plate with barcode validation
    and storage location tracking.
    """
    try:
        plate = service.create_compound_plate(
            db=db,
            barcode=data.barcode,
            plate_format=data.plate_format,
            plate_type=data.plate_type,
            storage_location_id=data.storage_location_id,
            created_by_id=current_user.id,
        )
        return schemas.CompoundPlateResponse.model_validate(plate)
    except ValueError as e:
        logger.error(f"Failed to create plate: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.get(
    "/plates", 
    response_model=List[schemas.CompoundPlateResponse],
    summary="List compound plates",
    description="List compound plates with optional filtering."
)
def list_plates(
    status: Optional[str] = Query(None, description="Filter by plate status"),
    storage_location_id: Optional[UUID] = Query(None, description="Filter by storage location"),
    limit: int = Query(100, ge=1, le=500, description="Maximum results"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """List compound plates with optional filters.
    
    Returns paginated list of compound plates with filtering
    options for status and storage location.
    """
    plates = service.list_plates(db, status=status, storage_location_id=storage_location_id, limit=limit)
    return [schemas.CompoundPlateResponse.model_validate(p) for p in plates]


@router.get(
    "/plates/{plate_id}", 
    response_model=dict,
    summary="Get plate with contents",
    description="Get plate information including all well contents."
)
def get_plate(
    plate_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get plate with well contents.
    
    Returns plate metadata along with all samples
    stored in the plate wells.
    """
    plate = service.get_plate(db, plate_id)
    if not plate:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail="Plate not found"
        )
    
    contents = service.get_plate_contents(db, plate_id)
    return {
        "plate": schemas.CompoundPlateResponse.model_validate(plate),
        "wells": [schemas.CompoundSampleResponse.model_validate(s) for s in contents],
    }


@router.get(
    "/plates/barcode/{barcode}", 
    response_model=schemas.CompoundPlateResponse,
    summary="Lookup plate by barcode",
    description="Find compound plate using barcode scan."
)
def get_plate_by_barcode(
    barcode: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Lookup plate by barcode.
    
    Enables barcode scanning for quick plate identification
    and retrieval in laboratory workflows.
    """
    plate = service.get_plate_by_barcode(db, barcode)
    if not plate:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"Plate with barcode '{barcode}' not found"
        )
    return schemas.CompoundPlateResponse.model_validate(plate)


# ============================================================================
# REQUEST ENDPOINTS (5)
# ============================================================================

@router.post(
    "/requests", 
    response_model=schemas.CompoundRequestResponse, 
    status_code=status.HTTP_201_CREATED,
    summary="Create compound request",
    description="Create a new request for compound samples."
)
def create_request(
    data: schemas.CompoundRequestCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Create a compound request.
    
    Creates request for specific sample or any sample of a compound.
    Supports flexible targeting and priority-based processing.
    """
    try:
        request = service.create_request(
            db=db,
            requester_id=current_user.id,
            sample_id=data.sample_id,
            compound_id=data.compound_id,
            requested_quantity=data.requested_quantity,
            quantity_unit=data.quantity_unit,
            purpose=data.purpose,
            priority=data.priority,
            notes=data.notes,
        )
        return schemas.CompoundRequestResponse.model_validate(request)
    except ValueError as e:
        logger.error(f"Failed to create request: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.get(
    "/requests", 
    response_model=List[schemas.CompoundRequestResponse],
    summary="List compound requests",
    description="List compound requests with optional filtering."
)
def list_requests(
    status: Optional[str] = Query(None, description="Filter by request status"),
    requester_id: Optional[UUID] = Query(None, description="Filter by requester"),
    pending_only: bool = Query(False, description="Show only pending requests"),
    limit: int = Query(100, ge=1, le=500, description="Maximum results"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """List requests with optional filters.
    
    Returns paginated list of compound requests with filtering
    options for status, requester, and pending-only view.
    """
    if pending_only:
        requests = service.list_pending_requests(db, limit=limit)
    else:
        requests = service.list_requests(db, status=status, requester_id=requester_id, limit=limit)
    return [schemas.CompoundRequestResponse.model_validate(r) for r in requests]


@router.get(
    "/requests/{request_id}", 
    response_model=schemas.CompoundRequestResponse,
    summary="Get compound request",
    description="Get detailed information about a specific request."
)
def get_request(
    request_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get request by ID.
    
    Returns detailed request information including approval
    and fulfillment status with timestamps.
    """
    request = service.get_request(db, request_id)
    if not request:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail="Request not found"
        )
    return schemas.CompoundRequestResponse.model_validate(request)


@router.put(
    "/requests/{request_id}/approve", 
    response_model=schemas.CompoundRequestResponse,
    summary="Approve request",
    description="Approve a compound request for fulfillment."
)
def approve_request(
    request_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Approve a compound request.
    
    Transitions request from 'requested' to 'approved' status
    with approval timestamp and user tracking.
    """
    try:
        request = service.approve_request(db, request_id, approved_by_id=current_user.id)
        return schemas.CompoundRequestResponse.model_validate(request)
    except ValueError as e:
        logger.error(f"Failed to approve request {request_id}: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.put(
    "/requests/{request_id}/fulfill", 
    response_model=schemas.CompoundRequestResponse,
    summary="Fulfill request",
    description="Fulfill request and deduct quantity from sample."
)
def fulfill_request(
    request_id: UUID,
    data: schemas.CompoundRequestFulfill,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Fulfill a compound request.
    
    Validates sufficient quantity, deducts from sample,
    and marks request as fulfilled with audit trail.
    """
    try:
        request = service.fulfill_request(
            db=db,
            request_id=request_id,
            fulfilled_by_id=current_user.id,
            sample_id=data.sample_id,
        )
        return schemas.CompoundRequestResponse.model_validate(request)
    except ValueError as e:
        logger.error(f"Failed to fulfill request {request_id}: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.put(
    "/requests/{request_id}/reject", 
    response_model=schemas.CompoundRequestResponse,
    summary="Reject request",
    description="Reject request with reason for audit trail."
)
def reject_request(
    request_id: UUID,
    data: schemas.CompoundRequestReject,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Reject a compound request.
    
    Rejects request with mandatory reason for audit trail
    and compliance tracking.
    """
    try:
        request = service.reject_request(
            db=db,
            request_id=request_id,
            rejected_by_id=current_user.id,
            rejection_reason=data.rejection_reason,
        )
        return schemas.CompoundRequestResponse.model_validate(request)
    except ValueError as e:
        logger.error(f"Failed to reject request {request_id}: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


# ============================================================================
# BARCODE LOOKUP (1)
# ============================================================================

@router.get(
    "/barcode/{barcode}", 
    response_model=schemas.BarcodeLookupResponse,
    summary="Unified barcode lookup",
    description="Lookup sample or plate by barcode scan."
)
def lookup_barcode(
    barcode: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Unified barcode lookup (sample or plate).
    
    Searches both samples and plates for barcode match,
    returning entity type and details for laboratory workflows.
    """
    result = service.lookup_barcode(db, barcode)
    
    response = {"type": result["type"], "sample": None, "plate": None}
    if result["type"] == "sample":
        response["sample"] = schemas.CompoundSampleResponse.model_validate(result["entity"])
    elif result["type"] == "plate":
        response["plate"] = schemas.CompoundPlateResponse.model_validate(result["entity"])
    
    return response


# ============================================================================
# ALERTS (3)
# ============================================================================

@router.get(
    "/alerts/low-stock", 
    response_model=List[schemas.CompoundSampleResponse],
    summary="Low stock alerts",
    description="Get samples below quantity threshold."
)
def get_low_stock_alerts(
    threshold: float = Query(100.0, ge=0, description="Quantity threshold"),
    limit: int = Query(100, ge=1, le=500, description="Maximum results"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get samples below quantity threshold.
    
    Returns samples that are running low on stock
    for proactive reordering and inventory management.
    """
    samples = service.get_low_stock_samples(db, threshold=threshold, limit=limit)
    return [schemas.CompoundSampleResponse.model_validate(s) for s in samples]


@router.get(
    "/alerts/expiring", 
    response_model=List[schemas.CompoundSampleResponse],
    summary="Expiring samples",
    description="Get samples expiring within specified days."
)
def get_expiring_alerts(
    days: int = Query(30, ge=1, le=365, description="Days until expiry"),
    limit: int = Query(100, ge=1, le=500, description="Maximum results"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get samples expiring within N days.
    
    Returns samples approaching expiry for proactive
    management and usage prioritization.
    """
    samples = service.get_expiring_samples(db, days=days, limit=limit)
    return [schemas.CompoundSampleResponse.model_validate(s) for s in samples]


@router.get(
    "/alerts/expired", 
    response_model=List[schemas.CompoundSampleResponse],
    summary="Expired samples",
    description="Get samples past expiry date for cleanup."
)
def get_expired_alerts(
    limit: int = Query(100, ge=1, le=500, description="Maximum results"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get samples past expiry date.
    
    Returns expired samples for cleanup, disposal,
    or status updates in inventory management.
    """
    samples = service.get_expired_samples(db, limit=limit)
    return [schemas.CompoundSampleResponse.model_validate(s) for s in samples]
