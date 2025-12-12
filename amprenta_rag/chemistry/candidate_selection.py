"""Candidate selection and TPP scoring service."""
from __future__ import annotations

from typing import Dict, Any, Optional
from uuid import UUID

from amprenta_rag.database.models import Compound, TargetProductProfile, CandidateNomination
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_traffic_light(value: float, min_val: float, max_val: float) -> str:
    """
    Calculate traffic light status for a value against bounds.
    
    Args:
        value: The value to check
        min_val: Minimum acceptable value
        max_val: Maximum acceptable value
        
    Returns:
        "green", "yellow", or "red"
    """
    if min_val <= value <= max_val:
        return "green"
    
    # Check if within 20% of bounds
    range_size = max_val - min_val
    if range_size == 0:
        return "red" if value != min_val else "green"
    
    tolerance = range_size * 0.2
    
    if (min_val - tolerance) <= value <= (max_val + tolerance):
        return "yellow"
    
    return "red"


def score_compound_against_tpp(compound_id: UUID, tpp_id: UUID, db) -> Dict[str, Any]:
    """
    Score a compound against a Target Product Profile.
    
    Args:
        compound_id: UUID of the compound
        tpp_id: UUID of the TPP
        db: Database session
        
    Returns:
        Dictionary with compound_id, tpp_id, criteria_scores, and overall_score
    """
    compound = db.query(Compound).filter(Compound.id == compound_id).first()
    tpp = db.query(TargetProductProfile).filter(TargetProductProfile.id == tpp_id).first()
    
    if not compound:
        raise ValueError(f"Compound {compound_id} not found")
    if not tpp:
        raise ValueError(f"TPP {tpp_id} not found")
    
    criteria_scores = []
    weighted_scores = []
    total_weight = 0.0
    
    # Property mapping from TPP property names to Compound attributes
    property_map = {
        "molecular_weight": "molecular_weight",
        "mw": "molecular_weight",
        "logp": "logp",
        "log_p": "logp",
        "hbd": "hbd_count",
        "hbd_count": "hbd_count",
        "hba": "hba_count",
        "hba_count": "hba_count",
        "rotatable_bonds": "rotatable_bonds",
        "tpsa": "tpsa",  # Topological Polar Surface Area (if available)
    }
    
    for criterion in tpp.criteria:
        if not isinstance(criterion, dict):
            continue
        
        property_name = criterion.get("property", "")
        min_val = criterion.get("min")
        max_val = criterion.get("max")
        weight = criterion.get("weight", 1.0)
        unit = criterion.get("unit", "")
        
        if min_val is None or max_val is None:
            continue
        
        # Get compound property value
        prop_key = property_map.get(property_name.lower(), property_name.lower())
        compound_value = getattr(compound, prop_key, None)
        
        if compound_value is None:
            # Try alternative property names
            for alt_key in property_map.values():
                if hasattr(compound, alt_key):
                    compound_value = getattr(compound, alt_key)
                    break
        
        if compound_value is None:
            logger.warning("[CANDIDATE] Property %s not found for compound %s", property_name, compound_id)
            continue
        
        # Calculate traffic light
        traffic_light = get_traffic_light(float(compound_value), float(min_val), float(max_val))
        
        # Calculate score (green=100, yellow=50, red=0)
        score = 100 if traffic_light == "green" else (50 if traffic_light == "yellow" else 0)
        
        criteria_scores.append({
            "property": property_name,
            "value": compound_value,
            "min": min_val,
            "max": max_val,
            "weight": weight,
            "unit": unit,
            "traffic_light": traffic_light,
            "score": score,
        })
        
        # Add to weighted average
        weighted_scores.append(score * weight)
        total_weight += weight
    
    # Calculate overall score
    if total_weight > 0:
        overall_score = sum(weighted_scores) / total_weight
    else:
        overall_score = 0.0
    
    return {
        "compound_id": str(compound_id),
        "tpp_id": str(tpp_id),
        "criteria_scores": criteria_scores,
        "overall_score": round(overall_score, 2),
    }


def nominate_candidate(
    compound_id: UUID,
    tpp_id: UUID,
    user_id: Optional[UUID],
    notes: Optional[str],
    db
) -> CandidateNomination:
    """
    Nominate a compound candidate against a TPP.
    
    Args:
        compound_id: UUID of the compound to nominate
        tpp_id: UUID of the TPP
        user_id: UUID of the user making the nomination
        notes: Optional notes
        db: Database session
        
    Returns:
        Created CandidateNomination object
    """
    # Score compound against TPP
    scoring_result = score_compound_against_tpp(compound_id, tpp_id, db)
    
    # Build scores dictionary
    scores_dict = {}
    for criterion in scoring_result["criteria_scores"]:
        property_name = criterion["property"]
        scores_dict[property_name] = {
            "value": criterion["value"],
            "score": criterion["score"],
            "traffic_light": criterion["traffic_light"],
        }
    
    # Create nomination
    nomination = CandidateNomination(
        compound_id=compound_id,
        tpp_id=tpp_id,
        scores=scores_dict,
        overall_score=scoring_result["overall_score"],
        status="nominated",
        notes=notes,
        nominated_by_id=user_id,
    )
    
    db.add(nomination)
    db.commit()
    db.refresh(nomination)
    
    logger.info("[CANDIDATE] Nominated compound %s against TPP %s (score=%.2f)", compound_id, tpp_id, scoring_result["overall_score"])
    
    return nomination
