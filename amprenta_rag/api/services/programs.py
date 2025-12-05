"""
CRUD services for Programs.
"""

from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session
from sqlalchemy import select

from amprenta_rag.database.models import Program as ProgramModel
from amprenta_rag.api.schemas import ProgramCreate, ProgramUpdate
import uuid


def create_program(db: Session, program: ProgramCreate) -> ProgramModel:
    """Create a new program."""
    db_program = ProgramModel(
        id=uuid.uuid4(),
        name=program.name,
        description=program.description,
        disease=program.disease or [],
    )
    db.add(db_program)
    db.commit()
    db.refresh(db_program)
    return db_program


def get_program(db: Session, program_id: UUID) -> Optional[ProgramModel]:
    """Get a program by ID."""
    return db.query(ProgramModel).filter(ProgramModel.id == program_id).first()


def get_programs(
    db: Session,
    skip: int = 0,
    limit: int = 100,
    name_filter: Optional[str] = None,
) -> List[ProgramModel]:
    """Get all programs with optional filtering."""
    query = db.query(ProgramModel)
    
    if name_filter:
        query = query.filter(ProgramModel.name.ilike(f"%{name_filter}%"))
    
    return query.offset(skip).limit(limit).all()


def update_program(
    db: Session,
    program_id: UUID,
    program: ProgramUpdate,
) -> Optional[ProgramModel]:
    """Update a program."""
    db_program = get_program(db, program_id)
    if not db_program:
        return None
    
    update_data = program.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        # Handle None values for list fields (convert to empty list)
        if field in ['disease'] and value is None:
            value = []
        setattr(db_program, field, value)
    
    db.commit()
    db.refresh(db_program)
    return db_program


def delete_program(db: Session, program_id: UUID) -> bool:
    """Delete a program."""
    db_program = get_program(db, program_id)
    if not db_program:
        return False
    
    db.delete(db_program)
    db.commit()
    return True

