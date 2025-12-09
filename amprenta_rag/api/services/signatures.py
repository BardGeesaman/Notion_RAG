"""
CRUD services for Signatures.
"""

from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import (
    Signature as SignatureModel,
    SignatureComponent as SignatureComponentModel,
    Feature as FeatureModel,
)
from amprenta_rag.api.schemas import SignatureCreate, SignatureUpdate
import uuid


def create_signature(db: Session, signature: SignatureCreate) -> SignatureModel:
    """Create a new signature."""
    db_signature = SignatureModel(
        id=uuid.uuid4(),
        name=signature.name,
        description=signature.description,
    )
    
    # Add components
    for comp in signature.components:
        db_component = SignatureComponentModel(
            id=uuid.uuid4(),
            signature_id=db_signature.id,
            feature_id=comp.feature_id,
            feature_name=comp.feature_name,
            feature_type=comp.feature_type.value,
            direction=comp.direction.value if comp.direction else None,
            weight=comp.weight or 1.0,
        )
        db_signature.components.append(db_component)
        
        # Link to feature if feature_id provided
        if comp.feature_id:
            feature = db.query(FeatureModel).filter(FeatureModel.id == comp.feature_id).first()
            if feature:
                db_signature.features.append(feature)
    
    # Auto-compute modalities from components
    modalities = list(set(comp.feature_type.value for comp in signature.components))
    db_signature.modalities = modalities
    
    # Add program relationships
    if signature.program_ids:
        from amprenta_rag.api.services.programs import get_program
        for program_id in signature.program_ids:
            program = get_program(db, program_id)
            if program:
                db_signature.programs.append(program)
    
    db.add(db_signature)
    db.commit()
    db.refresh(db_signature)
    return db_signature


def get_signature(db: Session, signature_id: UUID) -> Optional[SignatureModel]:
    """Get a signature by ID."""
    return db.query(SignatureModel).filter(SignatureModel.id == signature_id).first()


def get_signature_by_name(db: Session, name: str) -> Optional[SignatureModel]:
    """Get a signature by name."""
    return db.query(SignatureModel).filter(SignatureModel.name == name).first()


def get_signatures(
    db: Session,
    skip: int = 0,
    limit: int = 100,
    name_filter: Optional[str] = None,
    program_id: Optional[UUID] = None,
) -> List[SignatureModel]:
    """Get all signatures with optional filtering."""
    query = db.query(SignatureModel)
    
    if name_filter:
        query = query.filter(SignatureModel.name.ilike(f"%{name_filter}%"))
    
    if program_id:
        query = query.filter(SignatureModel.programs.any(id=program_id))
    
    return query.offset(skip).limit(limit).all()


def update_signature(
    db: Session,
    signature_id: UUID,
    signature: SignatureUpdate,
) -> Optional[SignatureModel]:
    """Update a signature."""
    db_signature = get_signature(db, signature_id)
    if not db_signature:
        return None
    
    update_data = signature.model_dump(exclude_unset=True)
    components = update_data.pop("components", None)
    program_ids = update_data.pop("program_ids", None)
    
    for field, value in update_data.items():
        setattr(db_signature, field, value)
    
    # Update components if provided
    if components is not None:
        # Clear existing components
        db_signature.components.clear()
        db_signature.features.clear()
        
        # Add new components
        modalities = []
        for comp in components:
            db_component = SignatureComponentModel(
                id=uuid.uuid4(),
                signature_id=db_signature.id,
                feature_id=comp.feature_id,
                feature_name=comp.feature_name,
                feature_type=comp.feature_type.value,
                direction=comp.direction.value if comp.direction else None,
                weight=comp.weight or 1.0,
            )
            db_signature.components.append(db_component)
            modalities.append(comp.feature_type.value)
            
            # Link to feature if feature_id provided
            if comp.feature_id:
                feature = db.query(FeatureModel).filter(FeatureModel.id == comp.feature_id).first()
                if feature:
                    db_signature.features.append(feature)
        
        # Update modalities
        db_signature.modalities = list(set(modalities))
    
    # Update program relationships if provided
    if program_ids is not None:
        from amprenta_rag.api.services.programs import get_program
        db_signature.programs.clear()
        for program_id in program_ids:
            program = get_program(db, program_id)
            if program:
                db_signature.programs.append(program)
    
    db.commit()
    db.refresh(db_signature)
    return db_signature


def delete_signature(db: Session, signature_id: UUID) -> bool:
    """Delete a signature."""
    db_signature = get_signature(db, signature_id)
    if not db_signature:
        return False
    
    db.delete(db_signature)
    db.commit()
    return True

