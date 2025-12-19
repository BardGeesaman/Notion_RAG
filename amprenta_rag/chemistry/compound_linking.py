"""Compound linking utilities using PostgreSQL."""
from typing import List
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Compound, Program
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def link_compound_to_signature(compound_id: str, signature_id: str) -> bool:
    """Link a compound to a signature."""
    db_gen = get_db()
    db = next(db_gen)
    try:
        compound = db.query(Compound).filter(Compound.compound_id == compound_id).first()
        if not compound:
            return False

        external_ids = compound.external_ids or {}
        signatures = external_ids.get("signatures", [])
        if signature_id not in signatures:
            signatures.append(signature_id)
        external_ids["signatures"] = signatures
        compound.external_ids = external_ids

        db.commit()
        return True
    except Exception as e:
        logger.warning("[COMPOUND_LINK] Failed linking compound %s to signature %s: %r", compound_id, signature_id, e)
        db.rollback()
        return False
    finally:
        db_gen.close()


def link_compound_to_program(compound_id: str, program_id: str) -> bool:
    """Link a compound to a program."""
    db_gen = get_db()
    db = next(db_gen)
    try:
        compound = db.query(Compound).filter(Compound.compound_id == compound_id).first()
        program = db.query(Program).filter(Program.id == UUID(program_id)).first()

        if not compound or not program:
            return False

        if program not in compound.programs:
            compound.programs.append(program)

        db.commit()
        return True
    except Exception as e:
        logger.warning("[COMPOUND_LINK] Failed linking compound %s to program %s: %r", compound_id, program_id, e)
        db.rollback()
        return False
    finally:
        db_gen.close()


def get_compounds_for_signature(signature_id: str) -> List[str]:
    """Get compound IDs linked to a signature."""
    db_gen = get_db()
    db = next(db_gen)
    try:
        compounds = db.query(Compound).all()
        result = []
        for c in compounds:
            ext_ids = c.external_ids or {}
            if signature_id in ext_ids.get("signatures", []):
                result.append(c.compound_id)
        return result
    finally:
        db_gen.close()


def get_compounds_for_program(program_id: str) -> List[str]:
    """Get compound IDs linked to a program."""
    db_gen = get_db()
    db = next(db_gen)
    try:
        program = db.query(Program).filter(Program.id == UUID(program_id)).first()
        if not program:
            return []
        return [c.compound_id for c in program.compounds]
    finally:
        db_gen.close()
