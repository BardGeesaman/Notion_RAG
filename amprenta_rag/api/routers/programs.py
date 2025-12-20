"""
API router for Programs.
"""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import Program, ProgramCreate, ProgramUpdate
from amprenta_rag.api.services import programs as program_service

router = APIRouter()


@router.post("/", response_model=Program, status_code=201)
async def create_program(
    program: ProgramCreate,
    db: Session = Depends(get_database_session),
):
    """Create a new program."""
    return program_service.create_program(db, program)


@router.get("/", response_model=List[Program])
async def list_programs(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    db: Session = Depends(get_database_session),
):
    """List all programs."""
    return program_service.get_programs(db, skip=skip, limit=limit, name_filter=name)


@router.get("/{program_id}", response_model=Program)
async def get_program(
    program_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Get a program by ID."""
    program = program_service.get_program(db, program_id)
    if not program:
        raise HTTPException(status_code=404, detail="Program not found")
    return program


@router.patch("/{program_id}", response_model=Program)
async def update_program(
    program_id: UUID,
    program: ProgramUpdate,
    db: Session = Depends(get_database_session),
):
    """Update a program."""
    updated = program_service.update_program(db, program_id, program)
    if not updated:
        raise HTTPException(status_code=404, detail="Program not found")
    return updated


@router.delete("/{program_id}", status_code=204)
async def delete_program(
    program_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Delete a program."""
    success = program_service.delete_program(db, program_id)
    if not success:
        raise HTTPException(status_code=404, detail="Program not found")

