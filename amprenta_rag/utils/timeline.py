"""Timeline utilities for unified activity timeline."""
from __future__ import annotations

from typing import List, Dict, Any, Optional
from datetime import datetime

from amprenta_rag.database.models import Experiment, Compound, Sample, DiscoveredStudy
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_timeline(limit: int = 50, user_id: Optional[str] = None, db = None) -> List[Dict[str, Any]]:
    """
    Get unified timeline of recent activities across entity types.
    
    Args:
        limit: Maximum number of items to return
        user_id: Optional user ID to filter by creator
        db: Database session
        
    Returns:
        List of dicts with type, id, name, timestamp, user, sorted by timestamp descending
    """
    timeline_items = []
    
    # Get recent experiments
    exp_query = db.query(Experiment).order_by(Experiment.created_at.desc())
    if user_id:
        from uuid import UUID
        exp_query = exp_query.filter(Experiment.created_by_id == (UUID(user_id) if isinstance(user_id, str) else user_id))
    experiments = exp_query.limit(limit).all()
    
    for exp in experiments:
        timeline_items.append({
            "type": "experiment",
            "id": str(exp.id),
            "name": exp.name,
            "timestamp": exp.created_at,
            "user": exp.created_by.username if exp.created_by else None,
        })
    
    # Get recent compounds
    comp_query = db.query(Compound).order_by(Compound.created_at.desc())
    compounds = comp_query.limit(limit).all()
    
    for comp in compounds:
        timeline_items.append({
            "type": "compound",
            "id": str(comp.id),
            "name": comp.compound_id,
            "timestamp": comp.created_at,
            "user": None,  # Compounds don't have created_by
        })
    
    # Get recent samples
    sample_query = db.query(Sample).order_by(Sample.created_at.desc())
    if user_id:
        from uuid import UUID
        sample_query = sample_query.filter(Sample.created_by_id == (UUID(user_id) if isinstance(user_id, str) else user_id))
    samples = sample_query.limit(limit).all()
    
    for sample in samples:
        timeline_items.append({
            "type": "sample",
            "id": str(sample.id),
            "name": sample.name,
            "timestamp": sample.created_at,
            "user": sample.created_by.username if sample.created_by else None,
        })
    
    # Get recent discoveries
    disc_query = db.query(DiscoveredStudy).order_by(DiscoveredStudy.discovered_at.desc())
    discoveries = disc_query.limit(limit).all()
    
    for disc in discoveries:
        timeline_items.append({
            "type": "discovery",
            "id": str(disc.id),
            "name": disc.study_id,
            "timestamp": disc.discovered_at,
            "user": None,  # Discoveries don't have created_by
        })
    
    # Sort by timestamp descending
    timeline_items.sort(key=lambda x: x["timestamp"] if x["timestamp"] else datetime.min, reverse=True)
    
    # Limit to requested number
    return timeline_items[:limit]
