"""
Program CRUD operations.

Split out of `amprenta_rag.database.crud` to keep domain operations smaller and
more maintainable. The facade module re-exports these functions for backwards
compatibility.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, cast
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Program
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_program(
    db: Session,
    name: str,
    description: Optional[str] = None,
    disease: Optional[List[str]] = None,
    notion_page_id: Optional[str] = None,
    external_ids: Optional[Dict[str, Any]] = None,
    commit: bool = True,
) -> Program:
    """
    Create a new program in the database.

    Args:
        db: Database session
        name: Program name
        description: Program description
        disease: List of disease terms
        notion_page_id: Notion page ID (for migration)
        external_ids: External system IDs
        commit: Whether to commit the transaction

    Returns:
        Created Program object
    """
    disease_list: List[str] = list(disease) if disease else []
    program = Program(
        name=name,
        description=description,
        disease=cast(Any, disease_list),
        notion_page_id=notion_page_id,
        external_ids=external_ids or {},
    )

    db.add(program)

    if commit:
        db.commit()
        db.refresh(program)
        logger.info(
            "[CRUD][PROGRAM] Created program: %s (ID: %s)",
            program.name,
            program.id,
        )

    return program


def get_program_by_id(db: Session, program_id: UUID) -> Optional[Program]:
    """Get program by UUID."""
    return db.query(Program).filter(Program.id == program_id).first()


def get_program_by_name(db: Session, name: str) -> Optional[Program]:
    """Get program by name."""
    return db.query(Program).filter(Program.name == name).first()


def get_program_by_notion_id(db: Session, notion_page_id: str) -> Optional[Program]:
    """Get program by Notion page ID."""
    return db.query(Program).filter(Program.notion_page_id == notion_page_id).first()


def get_or_create_program(
    db: Session,
    name: str,
    description: Optional[str] = None,
    disease: Optional[List[str]] = None,
    notion_page_id: Optional[str] = None,
) -> tuple[Program, bool]:
    """
    Get existing program or create new one.

    Returns:
        Tuple of (Program, created: bool)
    """
    # Try to find by notion_page_id first (if provided)
    if notion_page_id:
        program = get_program_by_notion_id(db, notion_page_id)
        if program:
            return program, False

    # Try to find by name
    program = get_program_by_name(db, name)
    if program:
        return program, False

    # Create new
    program = create_program(
        db,
        name=name,
        description=description,
        disease=disease,
        notion_page_id=notion_page_id,
    )
    return program, True


def list_programs(
    db: Session,
    limit: Optional[int] = None,
    offset: int = 0,
) -> List[Program]:
    """List all programs with optional pagination."""
    query = db.query(Program).order_by(Program.name)

    if limit:
        query = query.limit(limit).offset(offset)

    return query.all()


