"""Comparison utilities for comparing entities."""
from __future__ import annotations

from typing import Dict, Any
from uuid import UUID

from amprenta_rag.database.models import Experiment, Compound
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def compare_experiments(id1: str, id2: str, db) -> Dict[str, Any]:
    """
    Compare two experiments.

    Args:
        id1: UUID of first experiment
        id2: UUID of second experiment
        db: Database session

    Returns:
        Dict with item1, item2, and differences
    """
    exp1 = db.query(Experiment).filter(Experiment.id == (UUID(id1) if isinstance(id1, str) else id1)).first()
    exp2 = db.query(Experiment).filter(Experiment.id == (UUID(id2) if isinstance(id2, str) else id2)).first()

    if not exp1 or not exp2:
        raise ValueError("One or both experiments not found")

    item1 = {
        "id": str(exp1.id),
        "name": exp1.name,
        "description": exp1.description,
        "design_type": exp1.design_type,
    }

    item2 = {
        "id": str(exp2.id),
        "name": exp2.name,
        "description": exp2.description,
        "design_type": exp2.design_type,
    }

    # Find differences
    differences = []
    if exp1.name != exp2.name:
        differences.append("name")
    if exp1.description != exp2.description:
        differences.append("description")
    if exp1.design_type != exp2.design_type:
        differences.append("design_type")

    return {
        "item1": item1,
        "item2": item2,
        "differences": differences,
    }


def compare_compounds(id1: str, id2: str, db) -> Dict[str, Any]:
    """
    Compare two compounds.

    Args:
        id1: UUID of first compound
        id2: UUID of second compound
        db: Database session

    Returns:
        Dict with item1, item2, and differences
    """
    comp1 = db.query(Compound).filter(Compound.id == (UUID(id1) if isinstance(id1, str) else id1)).first()
    comp2 = db.query(Compound).filter(Compound.id == (UUID(id2) if isinstance(id2, str) else id2)).first()

    if not comp1 or not comp2:
        raise ValueError("One or both compounds not found")

    item1 = {
        "id": str(comp1.id),
        "compound_id": comp1.compound_id,
        "smiles": comp1.smiles,
        "molecular_weight": comp1.molecular_weight,
        "logp": comp1.logp,
    }

    item2 = {
        "id": str(comp2.id),
        "compound_id": comp2.compound_id,
        "smiles": comp2.smiles,
        "molecular_weight": comp2.molecular_weight,
        "logp": comp2.logp,
    }

    # Find differences
    differences = []
    if comp1.compound_id != comp2.compound_id:
        differences.append("compound_id")
    if comp1.smiles != comp2.smiles:
        differences.append("smiles")
    if comp1.molecular_weight != comp2.molecular_weight:
        differences.append("molecular_weight")
    if comp1.logp != comp2.logp:
        differences.append("logp")

    return {
        "item1": item1,
        "item2": item2,
        "differences": differences,
    }
