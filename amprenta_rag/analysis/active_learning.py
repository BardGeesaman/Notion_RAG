"""Active learning for compound screening prioritization."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional
from uuid import UUID

import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs, rdMolDescriptors

from amprenta_rag.database.models import Compound, MLModel
from amprenta_rag.database.session import db_session

logger = logging.getLogger(__name__)


@dataclass
class SuggestionResult:
    """Result of active learning compound suggestion."""
    
    compound_id: UUID
    smiles: str
    acquisition_score: float
    strategy_used: str
    rank: int
    explanation: str


@dataclass
class AcquisitionScore:
    """Acquisition score for a single compound."""
    
    compound_id: UUID
    score: float
    explanation: str


def suggest_next_compounds(
    screened: List[Dict[str, Any]],
    candidates: List[Dict[str, Any]],
    strategy: str,
    batch_size: int = 10,
    model_id: Optional[UUID] = None,
) -> List[SuggestionResult]:
    """
    Suggest next compounds to screen using active learning strategies.
    
    Args:
        screened: List of already screened compounds with results
                  [{"compound_id": UUID, "smiles": str, "activity": bool}, ...]
        candidates: List of candidate compounds to choose from
                    [{"compound_id": UUID, "smiles": str}, ...]
        strategy: Active learning strategy ("uncertainty" or "diversity")
        batch_size: Number of compounds to suggest
        model_id: Optional model ID for uncertainty-based strategies
        
    Returns:
        List of SuggestionResult objects ranked by acquisition score
    """
    logger.info(f"Suggesting {batch_size} compounds using {strategy} strategy from {len(candidates)} candidates")
    
    if not candidates:
        logger.warning("No candidate compounds provided")
        return []
    
    try:
        # Calculate acquisition scores for all candidates
        acquisition_scores = calculate_acquisition_scores(
            candidates=candidates,
            screened=screened,
            strategy=strategy,
            model_id=model_id,
        )
        
        if not acquisition_scores:
            logger.warning("No acquisition scores calculated")
            return []
        
        # Sort by score (descending - higher is better)
        acquisition_scores.sort(key=lambda x: x.score, reverse=True)
        
        # Take top batch_size compounds
        top_scores = acquisition_scores[:batch_size]
        
        # Convert to SuggestionResult objects
        suggestions = []
        for rank, score_obj in enumerate(top_scores, 1):
            # Find the compound info
            compound_info = next(
                (c for c in candidates if c["compound_id"] == score_obj.compound_id),
                None
            )
            
            if compound_info:
                suggestions.append(SuggestionResult(
                    compound_id=score_obj.compound_id,
                    smiles=compound_info["smiles"],
                    acquisition_score=score_obj.score,
                    strategy_used=strategy,
                    rank=rank,
                    explanation=score_obj.explanation,
                ))
        
        logger.info(f"Generated {len(suggestions)} suggestions using {strategy} strategy")
        return suggestions
        
    except Exception as e:
        logger.error(f"Error in suggest_next_compounds: {e}", exc_info=True)
        return []


def calculate_acquisition_scores(
    candidates: List[Dict[str, Any]],
    screened: List[Dict[str, Any]],
    strategy: str,
    model_id: Optional[UUID] = None,
) -> List[AcquisitionScore]:
    """
    Calculate acquisition scores for candidate compounds.
    
    Args:
        candidates: List of candidate compounds
        screened: List of already screened compounds
        strategy: Acquisition strategy ("uncertainty" or "diversity")
        model_id: Optional model ID for uncertainty-based strategies
        
    Returns:
        List of AcquisitionScore objects
    """
    logger.info(f"Calculating {strategy} acquisition scores for {len(candidates)} candidates")
    
    if strategy == "uncertainty":
        return _calculate_uncertainty_scores(candidates, model_id)
    elif strategy == "diversity":
        return _calculate_diversity_scores(candidates, screened)
    else:
        logger.error(f"Unknown acquisition strategy: {strategy}")
        return []


def _calculate_uncertainty_scores(
    candidates: List[Dict[str, Any]],
    model_id: Optional[UUID] = None,
) -> List[AcquisitionScore]:
    """Calculate uncertainty-based acquisition scores using model predictions."""
    scores = []
    
    try:
        if not model_id:
            logger.warning("No model_id provided for uncertainty sampling")
            return []
        
        # Load model from database
        with db_session() as db:
            ml_model = db.query(MLModel).filter(MLModel.id == model_id).first()
            if not ml_model:
                logger.error(f"Model {model_id} not found")
                return []
        
        # Import here to avoid circular imports
        from amprenta_rag.analysis.assay_predictor import predict_assay_outcome
        
        # Extract SMILES for prediction
        smiles_list = [c["smiles"] for c in candidates]
        
        # Get predictions with probabilities
        try:
            predictions = predict_assay_outcome(model_id, smiles_list)
            
            # Validate that predictions contain probability information
            if not predictions:
                logger.warning(f"No predictions returned from model {model_id}")
                return []
            
            # Check if first prediction has probability_active attribute
            first_prediction = predictions[0]
            if not hasattr(first_prediction, 'probability_active'):
                logger.error(f"Model {model_id} predictions do not contain probability information")
                return []
                
            # Validate probability range
            if not (0.0 <= first_prediction.probability_active <= 1.0):
                logger.error(f"Model {model_id} returned invalid probability: {first_prediction.probability_active}")
                return []
                
        except AttributeError as e:
            logger.error(f"Model {model_id} does not support probability predictions: {e}")
            return []
        except Exception as e:
            logger.error(f"Prediction failed for model {model_id}: {e}")
            return []
        
        # Calculate entropy (uncertainty) for each prediction
        for candidate, prediction in zip(candidates, predictions):
            prob_active = prediction.probability_active
            prob_inactive = 1.0 - prob_active
            
            # Calculate Shannon entropy
            # Entropy is maximized when probabilities are equal (maximum uncertainty)
            if prob_active > 0 and prob_inactive > 0:
                entropy = -(prob_active * np.log2(prob_active) + prob_inactive * np.log2(prob_inactive))
            else:
                entropy = 0.0  # No uncertainty if probability is 0 or 1
            
            # Normalize entropy to 0-1 scale (max entropy for binary is 1.0)
            normalized_entropy = entropy / 1.0
            
            explanation = f"Prediction uncertainty: {normalized_entropy:.3f} (prob_active: {prob_active:.3f})"
            
            scores.append(AcquisitionScore(
                compound_id=candidate["compound_id"],
                score=normalized_entropy,
                explanation=explanation,
            ))
        
        logger.info(f"Calculated uncertainty scores for {len(scores)} compounds")
        
    except Exception as e:
        logger.error(f"Error calculating uncertainty scores: {e}", exc_info=True)
    
    return scores


def _calculate_diversity_scores(
    candidates: List[Dict[str, Any]],
    screened: List[Dict[str, Any]],
) -> List[AcquisitionScore]:
    """Calculate diversity-based acquisition scores using Tanimoto distance."""
    scores = []
    
    try:
        # Get fingerprints for screened compounds
        screened_fps = []
        for screened_compound in screened:
            fp = _get_morgan_fingerprint(screened_compound["smiles"])
            if fp is not None:
                screened_fps.append(fp)
        
        if not screened_fps:
            logger.warning("No valid fingerprints from screened compounds, using random diversity")
            # If no screened compounds, assign random scores
            for candidate in candidates:
                scores.append(AcquisitionScore(
                    compound_id=candidate["compound_id"],
                    score=np.random.random(),
                    explanation="Random diversity (no screened compounds for comparison)",
                ))
            return scores
        
        # Calculate diversity scores for candidates
        for candidate in candidates:
            candidate_fp = _get_morgan_fingerprint(candidate["smiles"])
            
            if candidate_fp is None:
                # Invalid SMILES, assign low score
                scores.append(AcquisitionScore(
                    compound_id=candidate["compound_id"],
                    score=0.0,
                    explanation="Invalid SMILES - cannot calculate fingerprint",
                ))
                continue
            
            # Calculate Tanimoto similarities to all screened compounds
            similarities = []
            for screened_fp in screened_fps:
                similarity = DataStructs.TanimotoSimilarity(candidate_fp, screened_fp)
                similarities.append(similarity)
            
            # Diversity score = 1 - max_similarity (most diverse = least similar)
            max_similarity = max(similarities) if similarities else 0.0
            diversity_score = 1.0 - max_similarity
            
            explanation = f"Chemical diversity: {diversity_score:.3f} (max_similarity: {max_similarity:.3f})"
            
            scores.append(AcquisitionScore(
                compound_id=candidate["compound_id"],
                score=diversity_score,
                explanation=explanation,
            ))
        
        logger.info(f"Calculated diversity scores for {len(scores)} compounds")
        
    except Exception as e:
        logger.error(f"Error calculating diversity scores: {e}", exc_info=True)
    
    return scores


def _get_morgan_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048) -> Optional[Any]:
    """Get Morgan fingerprint for a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Generate Morgan fingerprint
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return fp
        
    except Exception as e:
        logger.warning(f"Failed to generate fingerprint for SMILES {smiles}: {e}")
        return None


def get_compound_suggestions_for_program(
    program_id: UUID,
    strategy: str = "uncertainty",
    batch_size: int = 10,
    model_id: Optional[UUID] = None,
) -> List[SuggestionResult]:
    """
    Get compound suggestions for a specific program.
    
    Convenience function that fetches screened and candidate compounds
    from the database for a given program.
    
    Args:
        program_id: UUID of the program
        strategy: Active learning strategy
        batch_size: Number of suggestions to return
        model_id: Optional model ID for uncertainty sampling
        
    Returns:
        List of SuggestionResult objects
    """
    try:
        with db_session() as db:
            # Get screened compounds (compounds with activity data)
            screened_query = (
                db.query(Compound)
                .join(Compound.biochemical_results)
                .filter(Compound.program_id == program_id)
                .distinct()
            )
            
            screened = []
            for compound in screened_query.all():
                # Determine activity from biochemical results
                activity = False
                for result in compound.biochemical_results:
                    if result.ic50 is not None and result.ic50 < 10.0:
                        activity = True
                        break
                    elif result.ki is not None and result.ki < 10.0:
                        activity = True
                        break
                    elif result.ec50 is not None and result.ec50 < 10.0:
                        activity = True
                        break
                
                screened.append({
                    "compound_id": compound.id,
                    "smiles": compound.smiles,
                    "activity": activity,
                })
            
            # Get candidate compounds (compounds without activity data)
            candidate_query = (
                db.query(Compound)
                .filter(Compound.program_id == program_id)
                .filter(~Compound.id.in_([s["compound_id"] for s in screened]))
                .filter(Compound.smiles.isnot(None))
            )
            
            candidates = []
            for compound in candidate_query.limit(1000).all():  # Limit for performance
                candidates.append({
                    "compound_id": compound.id,
                    "smiles": compound.smiles,
                })
            
            logger.info(f"Found {len(screened)} screened and {len(candidates)} candidate compounds for program {program_id}")
            
            # Generate suggestions
            return suggest_next_compounds(
                screened=screened,
                candidates=candidates,
                strategy=strategy,
                batch_size=batch_size,
                model_id=model_id,
            )
            
    except Exception as e:
        logger.error(f"Error getting suggestions for program {program_id}: {e}", exc_info=True)
        return []


def evaluate_strategy_performance(
    suggestions: List[SuggestionResult],
    true_activities: Dict[UUID, bool],
) -> Dict[str, float]:
    """
    Evaluate the performance of an active learning strategy.
    
    Args:
        suggestions: List of compound suggestions
        true_activities: Dictionary mapping compound_id to true activity
        
    Returns:
        Dictionary with performance metrics
    """
    if not suggestions or not true_activities:
        return {"hit_rate": 0.0, "enrichment": 0.0, "coverage": 0.0}
    
    # Calculate hit rate (fraction of suggested compounds that are active)
    suggested_ids = [s.compound_id for s in suggestions]
    hits = sum(1 for cid in suggested_ids if true_activities.get(cid, False))
    hit_rate = hits / len(suggested_ids) if suggested_ids else 0.0
    
    # Calculate enrichment vs random selection
    total_actives = sum(true_activities.values())
    total_compounds = len(true_activities)
    random_hit_rate = total_actives / total_compounds if total_compounds > 0 else 0.0
    enrichment = hit_rate / random_hit_rate if random_hit_rate > 0 else 0.0
    
    # Calculate coverage (fraction of suggested compounds with known activity)
    coverage = sum(1 for cid in suggested_ids if cid in true_activities) / len(suggested_ids)
    
    return {
        "hit_rate": hit_rate,
        "enrichment": enrichment,
        "coverage": coverage,
        "total_suggested": len(suggestions),
        "total_hits": hits,
    }


__all__ = [
    "SuggestionResult",
    "AcquisitionScore", 
    "suggest_next_compounds",
    "calculate_acquisition_scores",
    "get_compound_suggestions_for_program",
    "evaluate_strategy_performance",
]
