"""
Portfolio service for compound collection aggregation.

Aggregates ADMET, SAR, alerts, and MOBO data for portfolio view.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Compound
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_portfolio_summary(db: Session) -> Dict[str, Any]:
    """
    Get portfolio overview summary.

    Args:
        db: Database session

    Returns:
        Dictionary with compound counts, scaffold count, date range
    """
    try:
        total_compounds = db.query(Compound).count()
        
        # Get date range
        first_compound = db.query(Compound).order_by(Compound.created_at.asc()).first()
        latest_compound = db.query(Compound).order_by(Compound.created_at.desc()).first()
        
        date_from = first_compound.created_at if first_compound and first_compound.created_at else None
        date_to = latest_compound.created_at if latest_compound and latest_compound.created_at else None
        
        # Scaffold count (placeholder - would need scaffold analysis)
        scaffold_count = 0
        
        return {
            "total_compounds": total_compounds,
            "scaffold_count": scaffold_count,
            "date_from": date_from.isoformat() if date_from else None,
            "date_to": date_to.isoformat() if date_to else None,
        }
    
    except Exception as e:
        logger.error("[PORTFOLIO] Failed to get summary: %r", e)
        return {
            "total_compounds": 0,
            "scaffold_count": 0,
            "date_from": None,
            "date_to": None,
        }


def get_admet_summary(db: Session) -> Dict[str, int]:
    """
    Get ADMET traffic light summary.

    Args:
        db: Database session

    Returns:
        Dictionary with green, yellow, red counts
    """
    try:
        # Placeholder - would need alert_checker integration
        return {
            "green": 0,
            "yellow": 0,
            "red": 0,
            "unknown": db.query(Compound).count(),
        }
    
    except Exception as e:
        logger.error("[PORTFOLIO] Failed to get ADMET summary: %r", e)
        return {"green": 0, "yellow": 0, "red": 0, "unknown": 0}


def get_sar_gaps(db: Session, min_compounds: int = 3) -> List[Dict[str, Any]]:
    """
    Get scaffolds with sparse SAR coverage.

    Args:
        db: Database session
        min_compounds: Minimum compounds per scaffold

    Returns:
        List of scaffolds with <min_compounds compounds
    """
    try:
        # Placeholder - would need scaffold clustering
        return []
    
    except Exception as e:
        logger.error("[PORTFOLIO] Failed to get SAR gaps: %r", e)
        return []


def get_recommendations(db: Session, limit: int = 5) -> List[Dict[str, Any]]:
    """
    Get top compound recommendations.

    Args:
        db: Database session
        limit: Maximum recommendations to return

    Returns:
        List of recommended compounds with scores
    """
    try:
        # Get recent compounds as placeholder
        compounds = db.query(Compound).order_by(Compound.created_at.desc()).limit(limit).all()
        
        results = []
        for comp in compounds:
            results.append({
                "compound_id": str(comp.compound_id) if comp.compound_id else str(comp.id),
                "smiles": comp.canonical_smiles or comp.smiles,
                "score": 0.0,  # Placeholder - would need MOBO integration
                "reason": "Recent addition",
            })
        
        return results
    
    except Exception as e:
        logger.error("[PORTFOLIO] Failed to get recommendations: %r", e)
        return []

