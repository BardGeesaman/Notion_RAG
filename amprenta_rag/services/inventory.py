"""Compound inventory service layer."""

from __future__ import annotations

import logging
import uuid
from datetime import date, datetime, timedelta
from typing import List, Optional
from uuid import UUID

from sqlalchemy import and_, or_
from sqlalchemy.orm import Session

from amprenta_rag.database.models import CompoundPlate, CompoundRequest, Sample, SampleTransfer, User
from amprenta_rag.models.chemistry import Compound

logger = logging.getLogger(__name__)


# ============================================================================
# BARCODE GENERATION
# ============================================================================

def generate_barcode(prefix: str = "COMP") -> str:
    """Generate unique barcode in format PREFIX-YYYYMMDD-XXXXXX."""
    timestamp = datetime.utcnow().strftime("%Y%m%d")
    unique_id = uuid.uuid4().hex[:6].upper()
    return f"{prefix}-{timestamp}-{unique_id}"


# ============================================================================
# COMPOUND SAMPLE CRUD (uses extended Sample model)
# ============================================================================

def create_compound_sample(
    db: Session,
    compound_id: UUID,
    quantity: float,
    quantity_unit: str = "µL",
    concentration: Optional[float] = None,
    concentration_unit: Optional[str] = None,
    solvent: Optional[str] = None,
    format: Optional[str] = None,
    batch_lot: Optional[str] = None,
    expiry_date: Optional[date] = None,
    storage_location_id: Optional[UUID] = None,
    position: Optional[str] = None,
    plate_id: Optional[UUID] = None,
    well_position: Optional[str] = None,
    barcode: Optional[str] = None,
    created_by_id: Optional[UUID] = None,
    notes: Optional[str] = None,
) -> Sample:
    """
    Register a new compound sample (tube/vial).
    
    P1 FIX: Validates barcode uniqueness before creation.
    """
    # Generate barcode if not provided
    if not barcode:
        barcode = generate_barcode()
    
    # P1 FIX: Check barcode uniqueness
    existing = db.query(Sample).filter(Sample.barcode == barcode).first()
    if existing:
        raise ValueError(f"Barcode {barcode} already in use")
    
    # Verify compound exists
    compound = db.query(Compound).filter(Compound.id == compound_id).first()
    if not compound:
        raise ValueError(f"Compound {compound_id} not found")
    
    sample = Sample(
        name=f"{compound.compound_id} - {barcode}",
        sample_type="compound",
        barcode=barcode,
        compound_id=compound_id,
        quantity=quantity,
        unit=quantity_unit,
        concentration=concentration,
        concentration_unit=concentration_unit,
        solvent=solvent,
        format=format,
        batch_lot=batch_lot,
        expiry_date=expiry_date,
        storage_location_id=storage_location_id,
        position=position,
        plate_id=plate_id,
        well_position=well_position,
        status="available",
        created_by_id=created_by_id,
        notes=notes,
    )
    db.add(sample)
    db.commit()
    db.refresh(sample)
    logger.info(f"Created compound sample {sample.id} with barcode {barcode}")
    return sample


def get_compound_sample(db: Session, sample_id: UUID) -> Optional[Sample]:
    """Get compound sample by ID."""
    return db.query(Sample).filter(
        Sample.id == sample_id,
        Sample.sample_type == "compound"
    ).first()


def get_compound_sample_by_barcode(db: Session, barcode: str) -> Optional[Sample]:
    """Get compound sample by barcode."""
    return db.query(Sample).filter(Sample.barcode == barcode).first()


def list_compound_samples(
    db: Session,
    compound_id: Optional[UUID] = None,
    status: Optional[str] = None,
    storage_location_id: Optional[UUID] = None,
    plate_id: Optional[UUID] = None,
    limit: int = 200,
) -> List[Sample]:
    """List compound samples with optional filters."""
    query = db.query(Sample).filter(Sample.sample_type == "compound")
    
    if compound_id:
        query = query.filter(Sample.compound_id == compound_id)
    if status:
        query = query.filter(Sample.status == status)
    if storage_location_id:
        query = query.filter(Sample.storage_location_id == storage_location_id)
    if plate_id:
        query = query.filter(Sample.plate_id == plate_id)
    
    return query.order_by(Sample.created_at.desc()).limit(limit).all()


def update_compound_sample(
    db: Session,
    sample_id: UUID,
    quantity: Optional[float] = None,
    status: Optional[str] = None,
    storage_location_id: Optional[UUID] = None,
    position: Optional[str] = None,
    notes: Optional[str] = None,
) -> Sample:
    """Update compound sample fields."""
    sample = get_compound_sample(db, sample_id)
    if not sample:
        raise ValueError(f"Sample {sample_id} not found")
    
    if quantity is not None:
        sample.quantity = quantity
        if quantity == 0:
            sample.status = "depleted"
    if status is not None:
        sample.status = status
    if storage_location_id is not None:
        sample.storage_location_id = storage_location_id
    if position is not None:
        sample.position = position
    if notes is not None:
        sample.notes = notes
    
    db.commit()
    db.refresh(sample)
    return sample


def transfer_compound_sample(
    db: Session,
    sample_id: UUID,
    to_location_id: UUID,
    transferred_by_id: UUID,
    new_position: Optional[str] = None,
    notes: Optional[str] = None,
) -> SampleTransfer:
    """
    Transfer sample to new location and create audit record.
    
    P1 FIX: Uses existing SampleTransfer model for audit trail.
    """
    sample = get_compound_sample(db, sample_id)
    if not sample:
        raise ValueError(f"Sample {sample_id} not found")
    
    from_location_id = sample.storage_location_id
    
    # Create transfer record
    transfer = SampleTransfer(
        sample_id=sample_id,
        from_location_id=from_location_id,
        to_location_id=to_location_id,
        transferred_by_id=transferred_by_id,
        transferred_at=datetime.utcnow(),
        notes=notes,
    )
    db.add(transfer)
    
    # Update sample location
    sample.storage_location_id = to_location_id
    if new_position:
        sample.position = new_position
    
    db.commit()
    db.refresh(transfer)
    logger.info(f"Transferred sample {sample_id} from {from_location_id} to {to_location_id}")
    return transfer


def deplete_compound_sample(db: Session, sample_id: UUID) -> Sample:
    """Mark sample as depleted."""
    return update_compound_sample(db, sample_id, quantity=0, status="depleted")


# ============================================================================
# COMPOUND PLATE
# ============================================================================

def create_compound_plate(
    db: Session,
    barcode: str,
    plate_format: str,
    plate_type: str,
    storage_location_id: Optional[UUID] = None,
    created_by_id: Optional[UUID] = None,
) -> CompoundPlate:
    """Create a new compound plate."""
    # Check barcode uniqueness
    existing = db.query(CompoundPlate).filter(CompoundPlate.barcode == barcode).first()
    if existing:
        raise ValueError(f"Plate barcode {barcode} already in use")
    
    plate = CompoundPlate(
        barcode=barcode,
        plate_format=plate_format,
        plate_type=plate_type,
        storage_location_id=storage_location_id,
        status="active",
        created_by_id=created_by_id,
    )
    db.add(plate)
    db.commit()
    db.refresh(plate)
    logger.info(f"Created compound plate {plate.id} with barcode {barcode}")
    return plate


def get_plate(db: Session, plate_id: UUID) -> Optional[CompoundPlate]:
    """Get plate by ID."""
    return db.query(CompoundPlate).filter(CompoundPlate.id == plate_id).first()


def get_plate_by_barcode(db: Session, barcode: str) -> Optional[CompoundPlate]:
    """Get plate by barcode."""
    return db.query(CompoundPlate).filter(CompoundPlate.barcode == barcode).first()


def list_plates(
    db: Session,
    status: Optional[str] = None,
    storage_location_id: Optional[UUID] = None,
    limit: int = 100,
) -> List[CompoundPlate]:
    """List compound plates with optional filters."""
    query = db.query(CompoundPlate)
    if status:
        query = query.filter(CompoundPlate.status == status)
    if storage_location_id:
        query = query.filter(CompoundPlate.storage_location_id == storage_location_id)
    return query.order_by(CompoundPlate.created_at.desc()).limit(limit).all()


def get_plate_contents(db: Session, plate_id: UUID) -> List[Sample]:
    """Get all samples in a plate's wells."""
    return db.query(Sample).filter(
        Sample.plate_id == plate_id,
        Sample.sample_type == "compound"
    ).order_by(Sample.well_position).all()


# ============================================================================
# REQUEST WORKFLOW
# ============================================================================

def create_request(
    db: Session,
    requester_id: UUID,
    requested_quantity: float,
    quantity_unit: str = "µL",
    sample_id: Optional[UUID] = None,
    compound_id: Optional[UUID] = None,
    purpose: Optional[str] = None,
    priority: str = "normal",
    notes: Optional[str] = None,
) -> CompoundRequest:
    """Create a compound request."""
    if not sample_id and not compound_id:
        raise ValueError("Either sample_id or compound_id must be provided")
    
    request = CompoundRequest(
        sample_id=sample_id,
        compound_id=compound_id,
        requester_id=requester_id,
        requested_quantity=requested_quantity,
        quantity_unit=quantity_unit,
        purpose=purpose,
        priority=priority,
        status="requested",
        requested_at=datetime.utcnow(),
        notes=notes,
    )
    db.add(request)
    db.commit()
    db.refresh(request)
    logger.info(f"Created compound request {request.id}")
    return request


def approve_request(
    db: Session,
    request_id: UUID,
    approved_by_id: UUID,
) -> CompoundRequest:
    """Approve a compound request."""
    request = db.query(CompoundRequest).filter(CompoundRequest.id == request_id).first()
    if not request:
        raise ValueError(f"Request {request_id} not found")
    if request.status != "requested":
        raise ValueError(f"Request is not in 'requested' status (current: {request.status})")
    
    request.status = "approved"
    request.approved_at = datetime.utcnow()
    request.approved_by_id = approved_by_id
    
    db.commit()
    db.refresh(request)
    logger.info(f"Approved request {request_id}")
    return request


def fulfill_request(
    db: Session,
    request_id: UUID,
    fulfilled_by_id: UUID,
    sample_id: Optional[UUID] = None,
) -> CompoundRequest:
    """
    Fulfill a compound request and deduct quantity.
    
    P1 FIX: Validates sufficient quantity before fulfillment.
    """
    request = db.query(CompoundRequest).filter(CompoundRequest.id == request_id).first()
    if not request:
        raise ValueError(f"Request {request_id} not found")
    if request.status not in ("requested", "approved"):
        raise ValueError(f"Request cannot be fulfilled (current: {request.status})")
    
    # Determine which sample to use
    target_sample_id = sample_id or request.sample_id
    if not target_sample_id:
        raise ValueError("No sample specified for fulfillment")
    
    sample = get_compound_sample(db, target_sample_id)
    if not sample:
        raise ValueError(f"Sample {target_sample_id} not found")
    
    # P1 FIX: Validate quantity
    if sample.quantity < request.requested_quantity:
        raise ValueError(
            f"Insufficient quantity: {sample.quantity} {sample.unit} available, "
            f"{request.requested_quantity} {request.quantity_unit} requested"
        )
    
    # Deduct quantity
    sample.quantity -= request.requested_quantity
    if sample.quantity == 0:
        sample.status = "depleted"
    
    # Mark fulfilled
    request.status = "fulfilled"
    request.fulfilled_by_id = fulfilled_by_id
    request.fulfilled_at = datetime.utcnow()
    if not request.sample_id:
        request.sample_id = target_sample_id
    
    db.commit()
    db.refresh(request)
    logger.info(f"Fulfilled request {request_id} from sample {target_sample_id}")
    return request


def reject_request(
    db: Session,
    request_id: UUID,
    rejected_by_id: UUID,
    rejection_reason: str,
) -> CompoundRequest:
    """Reject a compound request with reason."""
    request = db.query(CompoundRequest).filter(CompoundRequest.id == request_id).first()
    if not request:
        raise ValueError(f"Request {request_id} not found")
    
    request.status = "rejected"
    request.rejection_reason = rejection_reason
    request.approved_by_id = rejected_by_id  # Reuse for rejection tracking
    request.approved_at = datetime.utcnow()
    
    db.commit()
    db.refresh(request)
    logger.info(f"Rejected request {request_id}: {rejection_reason}")
    return request


def list_pending_requests(db: Session, limit: int = 100) -> List[CompoundRequest]:
    """List pending requests (requested or approved, not yet fulfilled)."""
    return db.query(CompoundRequest).filter(
        CompoundRequest.status.in_(["requested", "approved"])
    ).order_by(
        CompoundRequest.priority.desc(),
        CompoundRequest.requested_at.asc()
    ).limit(limit).all()


def get_request(db: Session, request_id: UUID) -> Optional[CompoundRequest]:
    """Get request by ID."""
    return db.query(CompoundRequest).filter(CompoundRequest.id == request_id).first()


def list_requests(
    db: Session,
    status: Optional[str] = None,
    requester_id: Optional[UUID] = None,
    limit: int = 100,
) -> List[CompoundRequest]:
    """List requests with optional filters."""
    query = db.query(CompoundRequest)
    if status:
        query = query.filter(CompoundRequest.status == status)
    if requester_id:
        query = query.filter(CompoundRequest.requester_id == requester_id)
    return query.order_by(CompoundRequest.requested_at.desc()).limit(limit).all()


# ============================================================================
# COMPOUND PLATE
# ============================================================================

def create_compound_plate(
    db: Session,
    barcode: str,
    plate_format: str,
    plate_type: str,
    storage_location_id: Optional[UUID] = None,
    created_by_id: Optional[UUID] = None,
) -> CompoundPlate:
    """Create a new compound plate."""
    # Check barcode uniqueness
    existing = db.query(CompoundPlate).filter(CompoundPlate.barcode == barcode).first()
    if existing:
        raise ValueError(f"Plate barcode {barcode} already in use")
    
    plate = CompoundPlate(
        barcode=barcode,
        plate_format=plate_format,
        plate_type=plate_type,
        storage_location_id=storage_location_id,
        status="active",
        created_by_id=created_by_id,
    )
    db.add(plate)
    db.commit()
    db.refresh(plate)
    logger.info(f"Created compound plate {plate.id} with barcode {barcode}")
    return plate


def get_plate(db: Session, plate_id: UUID) -> Optional[CompoundPlate]:
    """Get plate by ID."""
    return db.query(CompoundPlate).filter(CompoundPlate.id == plate_id).first()


def get_plate_by_barcode(db: Session, barcode: str) -> Optional[CompoundPlate]:
    """Get plate by barcode."""
    return db.query(CompoundPlate).filter(CompoundPlate.barcode == barcode).first()


def list_plates(
    db: Session,
    status: Optional[str] = None,
    storage_location_id: Optional[UUID] = None,
    limit: int = 100,
) -> List[CompoundPlate]:
    """List compound plates with optional filters."""
    query = db.query(CompoundPlate)
    if status:
        query = query.filter(CompoundPlate.status == status)
    if storage_location_id:
        query = query.filter(CompoundPlate.storage_location_id == storage_location_id)
    return query.order_by(CompoundPlate.created_at.desc()).limit(limit).all()


def get_plate_contents(db: Session, plate_id: UUID) -> List[Sample]:
    """Get all samples in a plate's wells."""
    return db.query(Sample).filter(
        Sample.plate_id == plate_id,
        Sample.sample_type == "compound"
    ).order_by(Sample.well_position).all()


# ============================================================================
# BARCODE LOOKUP
# ============================================================================

def lookup_barcode(db: Session, barcode: str) -> dict:
    """
    Unified barcode lookup - returns sample or plate.
    
    Returns dict with 'type' ('sample' or 'plate') and 'entity'.
    """
    sample = get_compound_sample_by_barcode(db, barcode)
    if sample:
        return {"type": "sample", "entity": sample}
    
    plate = get_plate_by_barcode(db, barcode)
    if plate:
        return {"type": "plate", "entity": plate}
    
    return {"type": None, "entity": None}


# ============================================================================
# ALERT QUERIES (for dashboard)
# ============================================================================

def get_low_stock_samples(
    db: Session,
    threshold: float = 100.0,
    limit: int = 100,
) -> List[Sample]:
    """Get samples below quantity threshold."""
    return db.query(Sample).filter(
        Sample.sample_type == "compound",
        Sample.status == "available",
        Sample.quantity < threshold,
        Sample.quantity > 0,
    ).order_by(Sample.quantity.asc()).limit(limit).all()


def get_expiring_samples(
    db: Session,
    days: int = 30,
    limit: int = 100,
) -> List[Sample]:
    """Get samples expiring within N days."""
    cutoff = date.today() + timedelta(days=days)
    return db.query(Sample).filter(
        Sample.sample_type == "compound",
        Sample.status == "available",
        Sample.expiry_date.isnot(None),
        Sample.expiry_date <= cutoff,
        Sample.expiry_date >= date.today(),
    ).order_by(Sample.expiry_date.asc()).limit(limit).all()


def get_expired_samples(db: Session, limit: int = 100) -> List[Sample]:
    """Get samples past expiry date."""
    return db.query(Sample).filter(
        Sample.sample_type == "compound",
        Sample.expiry_date.isnot(None),
        Sample.expiry_date < date.today(),
    ).order_by(Sample.expiry_date.asc()).limit(limit).all()
