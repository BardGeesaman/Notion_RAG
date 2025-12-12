"""Dashboard widget utilities for displaying metrics."""
from __future__ import annotations

from datetime import datetime, timedelta

from amprenta_rag.database.models import Experiment, Compound, Sample, DiscoveredStudy


def get_experiment_count(db) -> int:
    """
    Get total count of experiments.
    
    Args:
        db: Database session
        
    Returns:
        Total number of experiments
    """
    return db.query(Experiment).count()


def get_compound_count(db) -> int:
    """
    Get total count of compounds.
    
    Args:
        db: Database session
        
    Returns:
        Total number of compounds
    """
    return db.query(Compound).count()


def get_sample_count(db) -> int:
    """
    Get total count of samples.
    
    Args:
        db: Database session
        
    Returns:
        Total number of samples
    """
    return db.query(Sample).count()


def get_discovery_count(days: int = 7, db=None) -> int:
    """
    Get count of discoveries in the last N days.
    
    Args:
        days: Number of days to look back (default: 7)
        db: Database session
        
    Returns:
        Number of discoveries in the last N days
    """
    if db is None:
        return 0
    
    cutoff = datetime.utcnow() - timedelta(days=days)
    return db.query(DiscoveredStudy).filter(
        DiscoveredStudy.discovered_at >= cutoff
    ).count()
