"""Multi-objective compound ranking service."""

from __future__ import annotations

from typing import Dict, List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.chemistry.alert_checker import StructuralAlertChecker
from amprenta_rag.database.models import BiochemicalResult, Compound
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.admet.predictor import get_admet_predictor
from amprenta_rag.ml.ranking.scorer import (
    MAX_COMPOUNDS,
    PRESETS,
    CompoundRanking,
    ObjectiveScore,
    compute_pareto_ranks,
    compute_weighted_score,
    normalize_alerts,
    normalize_herg,
    normalize_logp,
    normalize_logs,
    normalize_potency,
)

logger = get_logger(__name__)


def _get_best_potency(compound_id: UUID, db: Session) -> Optional[Dict]:
    """
    Get best (lowest) IC50/EC50 for compound from BiochemicalResult.
    
    Args:
        compound_id: UUID of compound
        db: Database session
        
    Returns:
        Dict with value, unit, type or None if no potency data
    """
    try:
        # Query all biochemical results for this compound
        results = (
            db.query(BiochemicalResult)
            .filter(BiochemicalResult.compound_id == compound_id)
            .all()
        )
        
        best_potency = None
        best_value = float('inf')
        
        for result in results:
            # Check IC50
            if result.ic50 is not None and result.ic50 > 0:
                if result.ic50 < best_value:
                    best_value = result.ic50
                    best_potency = {
                        "value": float(result.ic50),
                        "unit": getattr(result, "ic50_unit", "nM") or "nM",
                        "type": "ic50"
                    }
            
            # Check EC50
            if hasattr(result, "ec50") and result.ec50 is not None and result.ec50 > 0:
                if result.ec50 < best_value:
                    best_value = result.ec50
                    best_potency = {
                        "value": float(result.ec50),
                        "unit": getattr(result, "ec50_unit", "nM") or "nM", 
                        "type": "ec50"
                    }
            
            # Check Ki
            if hasattr(result, "ki") and result.ki is not None and result.ki > 0:
                if result.ki < best_value:
                    best_value = result.ki
                    best_potency = {
                        "value": float(result.ki),
                        "unit": getattr(result, "ki_unit", "nM") or "nM",
                        "type": "ki"
                    }
            
            # Check Kd
            if hasattr(result, "kd") and result.kd is not None and result.kd > 0:
                if result.kd < best_value:
                    best_value = result.kd
                    best_potency = {
                        "value": float(result.kd),
                        "unit": getattr(result, "kd_unit", "nM") or "nM",
                        "type": "kd"
                    }
        
        return best_potency
        
    except Exception as e:
        logger.warning("Failed to get potency for compound %s: %s", compound_id, e)
        return None


def _get_admet_predictions(smiles_list: List[str]) -> Dict[str, Dict]:
    """
    Batch predict ADMET for all compounds.
    
    Args:
        smiles_list: List of SMILES strings
        
    Returns:
        Dict keyed by SMILES with prediction results
    """
    if not smiles_list:
        return {}
    
    try:
        predictor = get_admet_predictor()
        results = predictor.predict(smiles_list, endpoints=["herg", "logs", "logp"])
        
        # Convert to dict keyed by SMILES
        predictions = {}
        for i, smiles in enumerate(smiles_list):
            if i < len(results.results):
                result = results.results[i]
                if result.error is None:
                    predictions[smiles] = {
                        "herg": result.predictions.get("herg", {}).get("mean"),
                        "logs": result.predictions.get("logs", {}).get("mean"), 
                        "logp": result.predictions.get("logp", {}).get("mean"),
                        "confidence": {
                            "herg": result.predictions.get("herg", {}).get("std"),
                            "logs": result.predictions.get("logs", {}).get("std"),
                            "logp": result.predictions.get("logp", {}).get("std"),
                        }
                    }
                else:
                    logger.warning("ADMET prediction failed for SMILES %s: %s", smiles, result.error)
                    predictions[smiles] = {"herg": None, "logs": None, "logp": None, "confidence": {}}
            else:
                logger.warning("Missing ADMET result for SMILES %s", smiles)
                predictions[smiles] = {"herg": None, "logs": None, "logp": None, "confidence": {}}
        
        return predictions
        
    except Exception as e:
        logger.warning("Batch ADMET prediction failed: %s", e)
        return {smiles: {"herg": None, "logs": None, "logp": None, "confidence": {}} for smiles in smiles_list}


def _get_structural_alerts(smiles_list: List[str]) -> Dict[str, List]:
    """
    Batch check structural alerts for all compounds.
    
    Args:
        smiles_list: List of SMILES strings
        
    Returns:
        Dict keyed by SMILES with alert lists
    """
    if not smiles_list:
        return {}
    
    try:
        checker = StructuralAlertChecker()
        alerts_dict = {}
        
        # Check if batch method exists
        if hasattr(checker, 'check_batch'):
            try:
                batch_results = checker.check_batch(smiles_list)
                for smiles, alerts in batch_results.items():
                    alerts_dict[smiles] = alerts if alerts else []
            except Exception as e:
                logger.warning("Batch alert checking failed, falling back to individual checks: %s", e)
                # Fallback to individual checks
                for smiles in smiles_list:
                    try:
                        alerts = checker.check_compound(smiles)
                        alerts_dict[smiles] = alerts if alerts else []
                    except Exception as e2:
                        logger.warning("Alert check failed for SMILES %s: %s", smiles, e2)
                        alerts_dict[smiles] = []
        else:
            # Individual checks
            for smiles in smiles_list:
                try:
                    alerts = checker.check_compound(smiles)
                    alerts_dict[smiles] = alerts if alerts else []
                except Exception as e:
                    logger.warning("Alert check failed for SMILES %s: %s", smiles, e)
                    alerts_dict[smiles] = []
        
        return alerts_dict
        
    except Exception as e:
        logger.warning("Structural alert checking failed: %s", e)
        return {smiles: [] for smiles in smiles_list}


def rank_compounds(
    compound_ids: List[UUID],
    weights: Dict[str, float],
    db: Session,
    include_pareto: bool = True
) -> List[CompoundRanking]:
    """
    Main ranking function for multi-objective compound optimization.
    
    Args:
        compound_ids: List of compound UUIDs to rank
        weights: Objective weights dict (e.g., {"potency": 0.4, "herg": 0.2, ...})
        db: Database session
        include_pareto: Whether to compute Pareto ranks
        
    Returns:
        List of CompoundRanking objects sorted by weighted score (best first)
    """
    # Validate batch size
    if len(compound_ids) > MAX_COMPOUNDS:
        logger.warning("Batch size %d exceeds limit %d, truncating", len(compound_ids), MAX_COMPOUNDS)
        compound_ids = compound_ids[:MAX_COMPOUNDS]
    
    if not compound_ids:
        return []
    
    try:
        # Fetch compounds from database
        compounds = (
            db.query(Compound)
            .filter(Compound.id.in_(compound_ids))
            .all()
        )
        
        if not compounds:
            logger.warning("No compounds found for provided IDs")
            return []
        
        # Extract SMILES list
        smiles_list = [c.smiles for c in compounds if c.smiles]
        
        if not smiles_list:
            logger.warning("No valid SMILES found for compounds")
            return []
        
        # Fetch potency data
        logger.info("Fetching potency data for %d compounds", len(compounds))
        potency_data = {}
        for compound in compounds:
            potency = _get_best_potency(compound.id, db)
            if potency:
                potency_data[str(compound.id)] = potency
        
        # Batch predict ADMET
        logger.info("Predicting ADMET properties for %d compounds", len(smiles_list))
        admet_predictions = _get_admet_predictions(smiles_list)
        
        # Batch check structural alerts
        logger.info("Checking structural alerts for %d compounds", len(smiles_list))
        alert_results = _get_structural_alerts(smiles_list)
        
        # Build rankings
        rankings = []
        
        for compound in compounds:
            compound_id_str = str(compound.id)
            smiles = compound.smiles
            
            if not smiles:
                logger.warning("Skipping compound %s: missing SMILES", compound_id_str)
                continue
            
            # Check for potency data (required)
            if compound_id_str not in potency_data:
                logger.info("Skipping compound %s: no potency data", compound_id_str)
                continue
            
            potency_info = potency_data[compound_id_str]
            admet_pred = admet_predictions.get(smiles, {})
            alerts = alert_results.get(smiles, [])
            
            # Build objective scores
            objectives = []
            
            # Potency objective
            potency_norm = normalize_potency(potency_info["value"], potency_info["unit"])
            objectives.append(ObjectiveScore(
                name="potency",
                raw_value=potency_info["value"],
                normalized=potency_norm,
                weight=weights.get("potency", 0.0),
                confidence=None
            ))
            
            # hERG objective
            herg_pred = admet_pred.get("herg")
            if herg_pred is not None:
                herg_norm = normalize_herg(herg_pred)
                herg_conf = admet_pred.get("confidence", {}).get("herg")
                objectives.append(ObjectiveScore(
                    name="herg",
                    raw_value=herg_pred,
                    normalized=herg_norm,
                    weight=weights.get("herg", 0.0),
                    confidence=herg_conf
                ))
            
            # logS objective
            logs_pred = admet_pred.get("logs")
            if logs_pred is not None:
                logs_norm = normalize_logs(logs_pred)
                logs_conf = admet_pred.get("confidence", {}).get("logs")
                objectives.append(ObjectiveScore(
                    name="logs",
                    raw_value=logs_pred,
                    normalized=logs_norm,
                    weight=weights.get("logs", 0.0),
                    confidence=logs_conf
                ))
            
            # logP objective
            logp_pred = admet_pred.get("logp")
            if logp_pred is not None:
                logp_norm = normalize_logp(logp_pred)
                logp_conf = admet_pred.get("confidence", {}).get("logp")
                objectives.append(ObjectiveScore(
                    name="logp",
                    raw_value=logp_pred,
                    normalized=logp_norm,
                    weight=weights.get("logp", 0.0),
                    confidence=logp_conf
                ))
            
            # Alerts objective
            alerts_norm = normalize_alerts(alerts)
            objectives.append(ObjectiveScore(
                name="alerts",
                raw_value=len(alerts),
                normalized=alerts_norm,
                weight=weights.get("alerts", 0.0),
                confidence=None
            ))
            
            # Compute weighted score
            weighted_score = compute_weighted_score(objectives)
            
            # Create ranking entry (pareto_rank will be filled later)
            ranking = CompoundRanking(
                compound_id=compound_id_str,
                smiles=smiles,
                objectives=objectives,
                weighted_score=weighted_score,
                pareto_rank=1,  # Placeholder
                rank=0  # Placeholder
            )
            rankings.append(ranking)
        
        # Sort by weighted score (descending - higher is better)
        rankings.sort(key=lambda r: r.weighted_score, reverse=True)
        
        # Assign overall ranks
        for i, ranking in enumerate(rankings):
            ranking.rank = i + 1
        
        # Compute Pareto ranks if requested
        if include_pareto and rankings:
            compound_objectives = []
            for ranking in rankings:
                obj_dict = {obj.name: obj.normalized for obj in ranking.objectives}
                compound_objectives.append({
                    "id": ranking.compound_id,
                    "objectives": obj_dict
                })
            
            pareto_ranks = compute_pareto_ranks(compound_objectives)
            for i, ranking in enumerate(rankings):
                ranking.pareto_rank = pareto_ranks[i]
        
        logger.info("Successfully ranked %d compounds", len(rankings))
        return rankings
        
    except Exception as e:
        logger.error("Compound ranking failed: %s", e)
        return []


def get_pareto_front(rankings: List[CompoundRanking]) -> List[CompoundRanking]:
    """
    Filter to Pareto rank 1 compounds only.
    
    Args:
        rankings: List of CompoundRanking objects
        
    Returns:
        List of compounds on the Pareto front (rank 1)
    """
    return [r for r in rankings if r.pareto_rank == 1]


def get_presets() -> List[Dict]:
    """
    Return available weight presets with descriptions.
    
    Returns:
        List of preset dictionaries with name, weights, and description
    """
    descriptions = {
        "balanced": "Balanced optimization across all objectives",
        "potency_first": "Prioritize potency over safety and ADMET",
        "safety_first": "Prioritize safety (hERG, alerts) over potency",
        "cns_optimized": "Optimized for CNS drug development (higher logP weight)"
    }
    
    return [
        {
            "name": name,
            "weights": weights,
            "description": descriptions.get(name, "Custom weight preset")
        }
        for name, weights in PRESETS.items()
    ]
