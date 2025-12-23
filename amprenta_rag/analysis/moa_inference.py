from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List, Set, Tuple
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


@dataclass
class BayesianMOACandidate:
    """MOA candidate with posterior uncertainty."""

    candidate_id: str
    type: str  # e.g., "target" or "pathway"
    probability: float  # posterior mean
    probability_ci: Tuple[float, float]  # 95% credible interval
    rank: int
    contributions: List[EvidenceContribution]

    def asdict(self) -> Dict[str, object]:
        return {
            "candidate_id": self.candidate_id,
            "type": self.type,
            "probability": self.probability,
            "probability_ci": self.probability_ci,
            "rank": self.rank,
            "contributions": [c.asdict() for c in self.contributions],
        }


def _require_pymc():
    try:
        import pymc as pm  # type: ignore
        import numpy as np  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError(
            "PyMC is required for Bayesian MOA inference. Install with: pymc>=5.0"
        ) from e
    return pm, np


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
                values.append(float(val))
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
    return {n for n in names if n is not None}


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
    for cand_name in candidates:
        s_bio = score_bioactivity_evidence(compound_id, cand_name)
        s_omics = score_omics_concordance(dataset_ids, cand_name)
        s_path = score_pathway_enrichment(features, cand_name)
        contributions = [
            EvidenceContribution(feature_name="bioactivity", value=s_bio, weight=0.4),
            EvidenceContribution(feature_name="omics_concordance", value=s_omics, weight=0.3),
            EvidenceContribution(feature_name="pathway_enrichment", value=s_path, weight=0.3),
        ]
        probability = fuse_evidence_scores(contributions)
        ranked.append(
            MOACandidate(
                candidate_id=cand_name,
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


def infer_moa_bayesian(
    compound_id: UUID,
    dataset_ids: List[UUID],
) -> List[BayesianMOACandidate]:
    """
    Bayesian logistic regression with evidence features.

    Implementation note (Phase 3 MVP):
    - We do not have ground-truth MOA labels; we treat the existing fusion score as a
      noisy observation of an underlying latent probability, and fit a Bayesian
      Beta regression with a logistic link over the 3 evidence features.
    - Outputs per-candidate posterior mean + 95% credible interval.
    """
    pm, np = _require_pymc()

    features = _collect_features(dataset_ids)
    candidates = generate_moa_candidates(compound_id)
    if not candidates:
        return []

    rows = []
    for cand_name in candidates:
        s_bio = score_bioactivity_evidence(compound_id, cand_name)
        s_omics = score_omics_concordance(dataset_ids, cand_name)
        s_path = score_pathway_enrichment(features, cand_name)
        contribs = [
            EvidenceContribution(feature_name="bioactivity", value=s_bio, weight=0.4),
            EvidenceContribution(feature_name="omics_concordance", value=s_omics, weight=0.3),
            EvidenceContribution(feature_name="pathway_enrichment", value=s_path, weight=0.3),
        ]
        y_fused = fuse_evidence_scores(contribs)
        rows.append((cand_name, contribs, [s_bio, s_omics, s_path], y_fused))

    # If too few candidates to fit a regression, return deterministic scores with a wide CI.
    if len(rows) < 3:
        out: List[BayesianMOACandidate] = []
        for cand_name, contribs, _, y in rows:
            lo = max(0.0, float(y) - 0.2)
            hi = min(1.0, float(y) + 0.2)
            out.append(
                BayesianMOACandidate(
                    candidate_id=cand_name,
                    type="target",
                    probability=round(float(y), 3),
                    probability_ci=(round(lo, 3), round(hi, 3)),
                    rank=0,
                    contributions=contribs,
                )
            )
        out.sort(key=lambda c: c.probability, reverse=True)
        for i, c in enumerate(out, start=1):
            c.rank = i
        return out

    X = np.asarray([r[2] for r in rows], dtype=float)  # (n, 3)
    y = np.asarray([r[3] for r in rows], dtype=float)
    y = np.clip(y, 1e-3, 1 - 1e-3)

    # Standardize features for sampling stability.
    X_mean = X.mean(axis=0, keepdims=True)
    X_std = np.clip(X.std(axis=0, keepdims=True), 1e-12, None)
    Xs = (X - X_mean) / X_std

    with pm.Model():
        beta0 = pm.Normal("beta0", mu=0.0, sigma=1.5)
        beta = pm.Normal("beta", mu=0.0, sigma=1.0, shape=(Xs.shape[1],))

        mu = pm.Deterministic("mu", pm.math.sigmoid(beta0 + pm.math.dot(Xs, beta)))
        kappa = pm.Exponential("kappa", lam=1.0)

        alpha = mu * kappa + 1e-3
        beta_param = (1 - mu) * kappa + 1e-3
        pm.Beta("y", alpha=alpha, beta=beta_param, observed=y)

        trace = pm.sample(
            draws=300,
            tune=300,
            chains=2,
            target_accept=0.9,
            progressbar=False,
            random_seed=42,
        )

    mu_samples = trace.posterior["mu"].values  # (chain, draw, n)
    mu_flat = mu_samples.reshape(-1, mu_samples.shape[-1])  # (samples, n)
    means = mu_flat.mean(axis=0)
    lo = np.quantile(mu_flat, 0.025, axis=0)
    hi = np.quantile(mu_flat, 0.975, axis=0)

    out2: List[BayesianMOACandidate] = []
    for idx, (cand_name, contribs, _, _) in enumerate(rows):
        out2.append(
            BayesianMOACandidate(
                candidate_id=cand_name,
                type="target",
                probability=round(float(means[idx]), 3),
                probability_ci=(round(float(lo[idx]), 3), round(float(hi[idx]), 3)),
                rank=0,
                contributions=contribs,
            )
        )

    out2.sort(key=lambda c: c.probability, reverse=True)
    for i, c in enumerate(out2, start=1):
        c.rank = i
    return out2

