from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Set
from uuid import UUID

from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
from amprenta_rag.database.models import BiochemicalResult, Feature, dataset_feature_assoc
from amprenta_rag.database.session import db_session


@dataclass
class EvidenceContribution:
    feature_name: str
    value: float
    weight: float

    def asdict(self) -> Dict[str, float | str]:
        return asdict(self)


@dataclass
class MOACandidate:
    candidate_id: str
    type: str  # e.g., "target" or "pathway"
    probability: float
    rank: int
    contributions: List[EvidenceContribution]

    def asdict(self) -> Dict[str, object]:
        return {
            "candidate_id": self.candidate_id,
            "type": self.type,
            "probability": self.probability,
            "rank": self.rank,
            "contributions": [c.asdict() for c in self.contributions],
        }


def generate_moa_candidates(compound_id: UUID) -> List[str]:
    """Return candidate targets based on biochemical results."""
    with db_session() as db:
        results: List[BiochemicalResult] = (
            db.query(BiochemicalResult)
            .filter(BiochemicalResult.compound_id == compound_id)
            .all()
        )
    targets = []
    for r in results:
        if r.target:
            targets.append(r.target)
    return list({t for t in targets if t})


def score_bioactivity_evidence(compound_id: UUID, candidate: str) -> float:
    """Score potency evidence for candidate target (lower IC50/EC50 => higher score)."""
    with db_session() as db:
        results: List[BiochemicalResult] = (
            db.query(BiochemicalResult)
            .filter(BiochemicalResult.compound_id == compound_id, BiochemicalResult.target == candidate)
            .all()
        )
    if not results:
        return 0.0
    values: List[float] = []
    for r in results:
        for val in [r.ic50, r.ec50, r.kd, r.ki]:
            if val is not None and val > 0:
                values.append(val)
    if not values:
        return 0.1
    best = min(values)
    # Simple heuristic: convert potency to 0-1 score (log-scale-like)
    score = max(0.0, min(1.0, 1.0 / (1.0 + best / 1000.0)))
    return score


def _collect_features(dataset_ids: List[UUID]) -> Set[str]:
    with db_session() as db:
        feats: List[Feature] = (
            db.query(Feature)
            .join(dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id)
            .filter(dataset_feature_assoc.c.dataset_id.in_(dataset_ids))
            .all()
        )
    names = {f.name for f in feats if getattr(f, "name", None)}
    return names


def score_omics_concordance(dataset_ids: List[UUID], candidate: str) -> float:
    """Score if candidate target appears in omics features."""
    features = _collect_features(dataset_ids)
    if not features:
        return 0.0
    match = candidate in features
    return 0.6 if match else 0.0


def score_pathway_enrichment(features: Set[str], candidate: str) -> float:
    """Use pathway enrichment with candidate as seed feature."""
    if not features:
        return 0.0
    input_features = set(features)
    input_features.add(candidate)
    try:
        results = perform_pathway_enrichment(
            input_features=input_features,
            input_feature_types={"gene"},
            p_value_threshold=0.1,
        )
    except Exception:
        return 0.0
    if not results:
        return 0.0
    # Higher score if any significant pathway
    best = min(r.adjusted_p_value for r in results if r.adjusted_p_value is not None)
    score = max(0.0, min(1.0, 1.0 - best))
    return score


def fuse_evidence_scores(contributions: List[EvidenceContribution]) -> float:
    """Combine evidence scores via weighted average."""
    if not contributions:
        return 0.0
    weighted_sum = sum(c.value * c.weight for c in contributions)
    weight_total = sum(c.weight for c in contributions)
    if weight_total == 0:
        return 0.0
    return max(0.0, min(1.0, weighted_sum / weight_total))


def infer_moa(compound_id: UUID, dataset_ids: List[UUID]) -> List[MOACandidate]:
    features = _collect_features(dataset_ids)
    candidates = generate_moa_candidates(compound_id)
    ranked: List[MOACandidate] = []
    for cand in candidates:
        s_bio = score_bioactivity_evidence(compound_id, cand)
        s_omics = score_omics_concordance(dataset_ids, cand)
        s_path = score_pathway_enrichment(features, cand)
        contributions = [
            EvidenceContribution(feature_name="bioactivity", value=s_bio, weight=0.4),
            EvidenceContribution(feature_name="omics_concordance", value=s_omics, weight=0.3),
            EvidenceContribution(feature_name="pathway_enrichment", value=s_path, weight=0.3),
        ]
        probability = fuse_evidence_scores(contributions)
        ranked.append(
            MOACandidate(
                candidate_id=cand,
                type="target",
                probability=round(probability, 3),
                rank=0,
                contributions=contributions,
            )
        )
    ranked.sort(key=lambda x: x.probability, reverse=True)
    for idx, cand in enumerate(ranked, start=1):
        cand.rank = idx
    return ranked

