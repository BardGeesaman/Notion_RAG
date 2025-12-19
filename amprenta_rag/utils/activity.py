"""Activity service for retrieving recent activity and statistics."""
from __future__ import annotations

from typing import List, Dict, Any

from amprenta_rag.database.models import Experiment, DiscoveredStudy, Compound, Dataset, User
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_recent_experiments(db, limit: int = 5) -> List[Dict[str, Any]]:
    """
    Get recent experiments ordered by creation date.

    Args:
        db: Database session
        limit: Maximum number of experiments to return

    Returns:
        List of dicts with id, name, created_at, created_by
    """
    experiments = (
        db.query(Experiment)
        .order_by(Experiment.created_at.desc())
        .limit(limit)
        .all()
    )

    results = []
    for exp in experiments:
        results.append({
            "id": str(exp.id),
            "name": exp.name,
            "created_at": exp.created_at.isoformat() if exp.created_at else None,
            "created_by": exp.created_by.username if exp.created_by else None,
        })

    return results


def get_recent_discoveries(db, limit: int = 5) -> List[Dict[str, Any]]:
    """
    Get recent discovered studies ordered by discovery date.

    Args:
        db: Database session
        limit: Maximum number of discoveries to return

    Returns:
        List of dicts with discovery information
    """
    discoveries = (
        db.query(DiscoveredStudy)
        .order_by(DiscoveredStudy.discovered_at.desc())
        .limit(limit)
        .all()
    )

    results = []
    for disc in discoveries:
        results.append({
            "id": str(disc.id),
            "study_id": disc.study_id,
            "repository": disc.repository,
            "title": disc.title,
            "discovered_at": disc.discovered_at.isoformat() if disc.discovered_at else None,
            "status": disc.status,
        })

    return results


def get_recent_compounds(db, limit: int = 5) -> List[Dict[str, Any]]:
    """
    Get recent compounds ordered by creation date.

    Args:
        db: Database session
        limit: Maximum number of compounds to return

    Returns:
        List of dicts with compound information
    """
    compounds = (
        db.query(Compound)
        .order_by(Compound.created_at.desc())
        .limit(limit)
        .all()
    )

    results = []
    for comp in compounds:
        results.append({
            "id": str(comp.id),
            "compound_id": comp.compound_id,
            "smiles": comp.smiles[:50] + "..." if comp.smiles and len(comp.smiles) > 50 else comp.smiles,
            "created_at": comp.created_at.isoformat() if comp.created_at else None,
        })

    return results


def get_activity_stats(db) -> Dict[str, int]:
    """
    Get activity statistics (counts of entities).

    Args:
        db: Database session

    Returns:
        Dict with total_experiments, total_compounds, total_datasets, total_users
    """
    stats = {
        "total_experiments": db.query(Experiment).count(),
        "total_compounds": db.query(Compound).count(),
        "total_datasets": db.query(Dataset).count(),
        "total_users": db.query(User).filter(User.is_active).count(),
    }

    return stats
