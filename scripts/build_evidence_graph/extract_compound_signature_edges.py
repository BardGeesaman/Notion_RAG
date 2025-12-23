"""Extract compound -> signature edges via feature overlap.

MVP approach:
- Build a signature feature index from `signature_feature` + Feature names.
- Build compound "target feature" sets from `activity_results` joined to biochemical_assays.target.
- Create edges where overlap > 0:
    compound(UUID) --[matches_signature]--> signature(UUID)

Confidence is computed as overlap / signature_feature_count (bounded [0,1]).

NOTE on ID Strategy:
    Both compound_id and signature_id are proper UUIDs from the Compound and
    Signature tables. Feature matching uses normalized names (uppercased) for
    comparison, but the resulting edges use the canonical entity UUIDs.
    Matched feature names are preserved in provenance for auditability.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from typing import Dict, List, Set, Tuple
from uuid import UUID

from amprenta_rag.database.models import ActivityResult, BiochemicalAssay, Feature, Signature, signature_feature_assoc
from amprenta_rag.database.session import db_session
from amprenta_rag.graph.edge_builder import EdgeBuilder


def _norm(s: str) -> str:
    return str(s or "").strip().upper()


def extract(
    limit_compounds: int = 50000,
    min_overlap: int = 2,
    max_edges_per_compound: int = 20,
) -> int:
    builder = EdgeBuilder()

    with db_session() as db:
        # Load signatures -> features and build inverted index
        sig_feats: Dict[UUID, Set[str]] = defaultdict(set)
        feat_to_sigs: Dict[str, Set[UUID]] = defaultdict(set)

        rows = (
            db.query(signature_feature_assoc.c.signature_id, Feature)
            .join(Feature, Feature.id == signature_feature_assoc.c.feature_id)
            .all()
        )
        for sig_id, feat in rows:
            nm = getattr(feat, "normalized_name", None) or getattr(feat, "name", None)
            if not nm:
                continue
            token = _norm(nm)
            if not token:
                continue
            sig_feats[sig_id].add(token)
            feat_to_sigs[token].add(sig_id)

        if not sig_feats:
            return 0

        # Preload signature name for provenance
        sig_names: Dict[UUID, str] = {s.id: s.name for s in db.query(Signature).all()}

        # Stream compound targets
        rows2 = (
            db.query(ActivityResult.compound_id, BiochemicalAssay.target)
            .join(BiochemicalAssay, BiochemicalAssay.id == ActivityResult.assay_id)
            .filter(BiochemicalAssay.target.isnot(None))
            .limit(limit_compounds)
            .all()
        )

        comp_targets: Dict[UUID, Set[str]] = defaultdict(set)
        for comp_id, target in rows2:
            tok = _norm(target)
            if tok:
                comp_targets[comp_id].add(tok)

        created = 0
        for comp_id, targets in comp_targets.items():
            counts: Dict[UUID, int] = defaultdict(int)
            matched_by_sig: Dict[UUID, Set[str]] = defaultdict(set)

            for t in targets:
                for sig_id in feat_to_sigs.get(t, set()):
                    counts[sig_id] += 1
                    matched_by_sig[sig_id].add(t)

            # Rank candidate signatures by overlap count
            ranked: List[Tuple[UUID, int]] = sorted(counts.items(), key=lambda kv: -kv[1])
            used = 0
            for sig_id, overlap in ranked:
                if overlap < min_overlap:
                    continue
                sig_size = len(sig_feats.get(sig_id, set())) or 1
                confidence = max(0.0, min(1.0, overlap / sig_size))
                builder.create_edge(
                    source_entity_type="compound",
                    source_entity_id=comp_id,
                    target_entity_type="signature",
                    target_entity_id=sig_id,
                    relationship_type="matches_signature",
                    confidence=confidence,
                    evidence_source="feature_overlap",
                    provenance={
                        "overlap": overlap,
                        "signature_feature_count": sig_size,
                        "compound_target_count": len(targets),
                        "matched_targets": sorted(list(matched_by_sig[sig_id]))[:200],
                        "signature_name": sig_names.get(sig_id),
                    },
                )
                created += 1
                used += 1
                if used >= max_edges_per_compound:
                    break

        return created


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit-compounds", type=int, default=50000)
    ap.add_argument("--min-overlap", type=int, default=2)
    ap.add_argument("--max-edges-per-compound", type=int, default=20)
    args = ap.parse_args()

    n = extract(
        limit_compounds=args.limit_compounds,
        min_overlap=args.min_overlap,
        max_edges_per_compound=args.max_edges_per_compound,
    )
    print(f"Created/updated {n} compound->signature edges.")


if __name__ == "__main__":
    main()


