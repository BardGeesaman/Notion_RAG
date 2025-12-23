"""Sphingolipid pathway imbalance scoring.

MVP constraints:
- Datasets do not store a unified numeric matrix in Postgres; we infer "signals" from feature names.
- If enzyme gene symbols are present in dataset features, we compute an enzyme-coverage score.
- Otherwise, we compute simple metabolite class ratios from lipid feature names (Cer/SM/HexCer).
"""

from __future__ import annotations

from typing import Any, Dict, List, Set
from uuid import UUID

from amprenta_rag.database.models import Feature, dataset_feature_assoc
from amprenta_rag.database.session import db_session
from amprenta_rag.ingestion.lipidomics.normalization import normalize_lipid_species

from .reactions import SPHINGOLIPID_ENZYMES, SPHINGOLIPID_REACTIONS


def _collect_feature_tokens(dataset_ids: List[UUID]) -> Set[str]:
    """Collect normalized tokens from Feature.name/normalized_name for datasets."""
    with db_session() as db:
        feats: List[Feature] = (
            db.query(Feature)
            .join(dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id)
            .filter(dataset_feature_assoc.c.dataset_id.in_(dataset_ids))
            .all()
        )

    tokens: Set[str] = set()
    for f in feats:
        for raw in (getattr(f, "name", None), getattr(f, "normalized_name", None)):
            if not raw:
                continue
            s = str(raw).strip()
            if not s:
                continue
            tokens.add(s.upper())
            # For lipids, add normalized class tokens
            norm = normalize_lipid_species(s)
            if norm:
                tokens.add(str(norm).upper())
                # Add lipid class token (Cer/SM/HexCer) when in canonical format
                cls = str(norm).split("(", 1)[0].strip()
                if cls:
                    tokens.add(cls.upper())
    return tokens


def _count_lipid_classes(tokens: Set[str]) -> Dict[str, int]:
    counts = {"CER": 0, "SM": 0, "HEXCER": 0, "S1P": 0, "SPHINGOSINE": 0}
    for t in tokens:
        u = t.upper()
        if u.startswith("CER(") or u == "CER":
            counts["CER"] += 1
        elif u.startswith("SM(") or u == "SM":
            counts["SM"] += 1
        elif u.startswith("HEXCER(") or u == "HEXCER":
            counts["HEXCER"] += 1
        elif u in {"S1P", "SPH-1-P"}:
            counts["S1P"] += 1
        elif u == "SPHINGOSINE":
            counts["SPHINGOSINE"] += 1
    return counts


def compute_pathway_imbalance(
    dataset_ids: List[UUID],
    pathway: str = "ceramide",
) -> Dict[str, Any]:
    """Compute a heuristic pathway imbalance score for sphingolipid metabolism.

    Returns a dict with:
    - score: float (0..1), higher means more ceramide-skewed imbalance (for ceramide pathway)
    - method: "enzymes" or "ratios"
    - matched: details of matched enzymes/metabolites
    - stats: counts and component ratios
    """
    if not dataset_ids:
        raise ValueError("dataset_ids must be non-empty")

    p = (pathway or "ceramide").lower().strip()
    tokens = _collect_feature_tokens(dataset_ids)

    # Enzyme-based scoring: fraction of key enzymes observed.
    matched_enzymes = sorted([e for e in SPHINGOLIPID_ENZYMES if e.upper() in tokens])
    if matched_enzymes:
        # Reaction coverage: how many reactions have at least one enzyme observed.
        covered_reactions = []
        for rid, rinfo in SPHINGOLIPID_REACTIONS.items():
            enz = [str(x).upper() for x in (rinfo.get("enzymes") or [])]  # type: ignore[union-attr]
            if any(e in tokens for e in enz):
                covered_reactions.append(rid)

        enzyme_coverage = len(matched_enzymes) / max(1, len(SPHINGOLIPID_ENZYMES))
        reaction_coverage = len(covered_reactions) / max(1, len(SPHINGOLIPID_REACTIONS))

        # Map to a 0..1 score; for now, average both coverages.
        score = max(0.0, min(1.0, 0.5 * enzyme_coverage + 0.5 * reaction_coverage))
        return {
            "pathway": p,
            "method": "enzymes",
            "score": float(score),
            "matched": {"enzymes": matched_enzymes, "covered_reactions": covered_reactions},
            "stats": {
                "enzyme_coverage": enzyme_coverage,
                "reaction_coverage": reaction_coverage,
                "n_tokens": len(tokens),
            },
        }

    # Fallback: metabolite ratios from lipid classes.
    counts = _count_lipid_classes(tokens)
    cer = counts["CER"]
    sm = counts["SM"]
    hexcer = counts["HEXCER"]

    # Ratio proxy: Cer abundance relative to (SM + HexCer)
    ratio = (cer + 1.0) / (sm + hexcer + 1.0)
    # Squash to 0..1
    score = ratio / (1.0 + ratio)

    return {
        "pathway": p,
        "method": "ratios",
        "score": float(score),
        "matched": {"lipid_class_counts": counts},
        "stats": {"ratio_cer_to_sm_hexcer": float(ratio), "n_tokens": len(tokens)},
    }


__all__ = ["compute_pathway_imbalance"]


