"""Multi-objective compound scoring and ranking."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

MAX_COMPOUNDS = 500

# Weight presets for different optimization strategies
PRESETS = {
    "balanced": {"potency": 0.35, "herg": 0.20, "alerts": 0.20, "logs": 0.15, "logp": 0.10},
    "potency_first": {"potency": 0.60, "herg": 0.10, "alerts": 0.10, "logs": 0.10, "logp": 0.10},
    "safety_first": {"potency": 0.15, "herg": 0.25, "alerts": 0.25, "logs": 0.20, "logp": 0.15},
    "cns_optimized": {"potency": 0.30, "herg": 0.20, "alerts": 0.15, "logs": 0.15, "logp": 0.20},
}


@dataclass
class ObjectiveScore:
    """Individual objective score for compound ranking."""
    
    name: str          # "potency", "herg", "alerts", "logs", "logp"
    raw_value: float
    normalized: float  # 0-1 (higher = better)
    weight: float
    confidence: Optional[float] = None


@dataclass
class CompoundRanking:
    """Complete ranking result for a compound."""
    
    compound_id: str
    smiles: str
    objectives: List[ObjectiveScore]
    weighted_score: float
    pareto_rank: int
    rank: int  # Overall rank by weighted_score


def normalize_potency(value: float, unit: str) -> float:
    """
    Normalize potency value to 0-1 scale (higher = better).
    
    Converts to nM, calculates pIC50, then normalizes.
    Range: 3 (1mM) to 10 (0.1pM) pIC50 scale.
    
    Args:
        value: Potency value
        unit: Unit string ("nM", "µM", "uM", "pM", "M")
        
    Returns:
        Normalized score 0-1 (higher = more potent)
    """
    if value <= 0:
        return 0.0
    
    # Convert to nM
    unit_lower = unit.lower().strip()
    if unit_lower in ("nm", "nanomolar"):
        value_nm = float(value)
    elif unit_lower in ("µm", "um", "micromolar", "μm"):
        value_nm = float(value) * 1000
    elif unit_lower in ("pm", "picomolar"):
        value_nm = float(value) / 1000
    elif unit_lower in ("m", "molar"):
        value_nm = float(value) * 1e9
    elif unit_lower in ("mm", "millimolar"):
        value_nm = float(value) * 1e6
    else:
        logger.warning("Unknown potency unit '%s', assuming nM", unit)
        value_nm = float(value)
    
    # Convert to pIC50: -log10(concentration in M)
    try:
        pic50 = -math.log10(value_nm * 1e-9)
    except (ValueError, OverflowError):
        return 0.0
    
    # Normalize: 3 (1mM) to 10 (0.1pM) -> 0 to 1
    normalized = (pic50 - 3.0) / 7.0
    return max(0.0, min(1.0, normalized))


def normalize_herg(pred_value: float) -> float:
    """
    Normalize hERG prediction to 0-1 scale (higher = safer).
    
    Lower hERG blocking = safer, so we invert the prediction.
    
    Args:
        pred_value: hERG prediction (0-1 scale or pIC50)
        
    Returns:
        Normalized score 0-1 (higher = safer)
    """
    # Assume pred_value is 0-1 probability of hERG blocking
    # If it's pIC50, we'd need different handling
    if pred_value < 0:
        pred_value = 0.0
    elif pred_value > 1:
        pred_value = 1.0
    
    # Invert: lower blocking = higher safety score
    return 1.0 - pred_value


def normalize_logs(pred_value: float) -> float:
    """
    Normalize logS prediction to 0-1 scale (higher = better solubility).
    
    Typical logS range: -6 (poorly soluble) to 0 (highly soluble).
    
    Args:
        pred_value: logS prediction
        
    Returns:
        Normalized score 0-1 (higher = more soluble)
    """
    # Normalize from typical range -6 to 0
    normalized = (pred_value + 6.0) / 6.0
    return max(0.0, min(1.0, normalized))


def normalize_logp(pred_value: float) -> float:
    """
    Normalize logP prediction to 0-1 scale (higher = better for oral drugs).
    
    Ideal range for oral drugs: 1-3. Penalize deviation from this range.
    
    Args:
        pred_value: logP prediction
        
    Returns:
        Normalized score 0-1 (higher = better for oral bioavailability)
    """
    # Ideal target: logP = 2
    # Penalty for deviation from ideal range
    deviation = abs(pred_value - 2.0)
    normalized = 1.0 - (deviation / 4.0)  # Penalty scales with deviation
    return max(0.0, min(1.0, normalized))


def normalize_alerts(alerts: List[Dict], alert_weights: Optional[Dict[str, float]] = None) -> float:
    """
    Normalize structural alerts to 0-1 scale (higher = safer).
    
    Args:
        alerts: List of alert dictionaries with 'type' and optional 'count'
        alert_weights: Custom weights per alert type
        
    Returns:
        Normalized score 0-1 (higher = fewer/less severe alerts)
    """
    if not alerts:
        return 1.0
    
    # Default alert weights
    if alert_weights is None:
        alert_weights = {"PAINS": 1.0, "BRENK": 0.5, "LILLY": 0.5}
    
    # Count weighted alerts
    weighted_sum = 0.0
    for alert in alerts:
        alert_type = alert.get("type", "").upper()
        count = alert.get("count", 1)
        weight = alert_weights.get(alert_type, 0.3)  # Default weight for unknown types
        weighted_sum += count * weight
    
    # Exponential decay: more alerts = lower score
    # Scale factor of 5 means ~5 PAINS alerts -> ~37% safety score
    safety_score = math.exp(-weighted_sum / 5.0)
    return max(0.0, min(1.0, safety_score))


def compute_weighted_score(objectives: List[ObjectiveScore]) -> float:
    """
    Compute weighted sum of normalized objective scores.
    
    Args:
        objectives: List of ObjectiveScore instances
        
    Returns:
        Weighted score (0-1 scale)
    """
    if not objectives:
        return 0.0
    
    total_weighted = sum(obj.normalized * obj.weight for obj in objectives)
    total_weight = sum(obj.weight for obj in objectives)
    
    if total_weight == 0:
        return 0.0
    
    return total_weighted / total_weight


def compute_pareto_ranks(compounds: List[Dict]) -> List[int]:
    """
    Compute Pareto ranks using non-dominated sorting.
    
    Args:
        compounds: List of dicts with {"id": str, "objectives": {name: normalized_value}}
        
    Returns:
        List of Pareto ranks (1 = Pareto front, 2 = second front, etc.)
    """
    if not compounds:
        return []
    
    n = len(compounds)
    if n == 1:
        return [1]
    
    # Extract objective values for each compound
    objectives_matrix = []
    for comp in compounds:
        obj_dict = comp.get("objectives", {})
        # Convert to list of values in consistent order
        obj_names = sorted(obj_dict.keys())
        obj_values = [obj_dict.get(name, 0.0) for name in obj_names]
        objectives_matrix.append(obj_values)
    
    # Non-dominated sorting
    ranks = [0] * n
    current_rank = 1
    remaining = list(range(n))
    
    while remaining:
        # Find non-dominated solutions in remaining set
        front = []
        
        for i in remaining:
            is_dominated = False
            for j in remaining:
                if i == j:
                    continue
                    
                # Check if j dominates i
                # j dominates i if j is >= i in all objectives and > i in at least one
                all_geq = True
                any_gt = False
                
                for k in range(len(objectives_matrix[i])):
                    if objectives_matrix[j][k] < objectives_matrix[i][k]:
                        all_geq = False
                        break
                    elif objectives_matrix[j][k] > objectives_matrix[i][k]:
                        any_gt = True
                
                if all_geq and any_gt:
                    is_dominated = True
                    break
            
            if not is_dominated:
                front.append(i)
        
        # Assign rank to current front
        for i in front:
            ranks[i] = current_rank
        
        # Remove front from remaining
        remaining = [i for i in remaining if i not in front]
        current_rank += 1
        
        # Safety check to avoid infinite loop
        if current_rank > n:
            logger.warning("Pareto ranking may have failed, assigning remaining compounds to rank %d", current_rank)
            for i in remaining:
                ranks[i] = current_rank
            break
    
    return ranks
