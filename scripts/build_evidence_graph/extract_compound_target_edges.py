"""Extract compound -> target edges from chemistry assay/activity tables.

MVP sources:
- activity_results -> biochemical_assays.target

We represent targets as `Feature` records (feature_type="gene") and create edges:
  compound(UUID) --[activity_against]--> feature(UUID)

NOTE on Target ID Strategy:
    BiochemicalAssay.target is a string (gene name), but GraphEdge requires UUIDs.
    We resolve this by creating or reusing Feature records with feature_type="gene"
    for each unique target name. This allows proper UUID-based graph edges while
    preserving the target name in provenance metadata.

    Future improvement: Create a dedicated Target entity table with proper
    normalization (UniProt IDs, gene symbols, aliases).
"""

from __future__ import annotations

import argparse
from typing import Dict, List, Optional, Tuple
from uuid import UUID

from amprenta_rag.database.models import ActivityResult, BiochemicalAssay, Compound, Feature
from amprenta_rag.database.session import db_session
from amprenta_rag.graph.edge_builder import EdgeBuilder


def _normalize_to_nm(value: float, unit: str | None) -> Optional[float]:
    if value is None:
        return None
    u = (unit or "").strip().lower()
    if not u:
        return None
    if u == "nm":
        return float(value)
    if u in ("um", "µm"):
        return float(value) * 1_000.0
    if u == "mm":
        return float(value) * 1_000_000.0
    if u == "m":
        return float(value) * 1_000_000_000.0
    return None


def _get_or_create_target_feature(db, target_name: str) -> Feature:
    f = (
        db.query(Feature)
        .filter(Feature.feature_type == "gene", Feature.name == target_name)
        .first()
    )
    if f is not None:
        return f
    f = Feature(name=target_name, feature_type="gene", normalized_name=target_name)
    db.add(f)
    db.flush()
    return f


def extract(limit: int = 100000) -> int:
    builder = EdgeBuilder()
    created = 0

    with db_session() as db:
        # Collect per (compound_id, target) activity stats.
        rows = (
            db.query(
                ActivityResult.compound_id,
                BiochemicalAssay.target,
                ActivityResult.value,
                ActivityResult.assay_id,
                BiochemicalAssay.assay_type,
                BiochemicalAssay.unit,
            )
            .join(BiochemicalAssay, BiochemicalAssay.id == ActivityResult.assay_id)
            .filter(BiochemicalAssay.target.isnot(None))
            .limit(limit)
            .all()
        )

        pair_counts: Dict[Tuple[UUID, str], int] = {}
        pair_assay_ids: Dict[Tuple[UUID, str], set] = {}
        pair_values_nm: Dict[Tuple[UUID, str], List[float]] = {}
        pair_activity_types: Dict[Tuple[UUID, str], Dict[str, int]] = {}

        for compound_id, target, value, assay_id, assay_type, unit in rows:
            if not compound_id or not target or assay_id is None:
                continue
            key = (compound_id, str(target).strip())
            if not key[1]:
                continue
            pair_counts[key] = pair_counts.get(key, 0) + 1
            pair_assay_ids.setdefault(key, set()).add(str(assay_id))

            at = str(assay_type).strip() if assay_type is not None else ""
            if at:
                d = pair_activity_types.setdefault(key, {})
                d[at] = d.get(at, 0) + 1

            nm = _normalize_to_nm(float(value), unit)
            if nm is not None:
                pair_values_nm.setdefault(key, []).append(float(nm))

        for (compound_id, target_name), cnt in pair_counts.items():
            # Ensure compound exists
            comp = db.query(Compound).filter(Compound.id == compound_id).first()
            if comp is None:
                continue
            feat = _get_or_create_target_feature(db, target_name)
            db.commit()  # commit created features progressively

            edge = builder.create_edge(
                source_entity_type="compound",
                source_entity_id=compound_id,
                target_entity_type="feature",
                target_entity_id=feat.id,
                relationship_type="activity_against",
                confidence=min(1.0, 0.2 + 0.1 * cnt),
                evidence_source="activity_results",
                provenance={
                    "assays_count": cnt,
                    "target_name": target_name,
                    "assay_ids": sorted(list(pair_assay_ids.get((compound_id, target_name), set()))),
                    "activity_type": (
                        max((pair_activity_types.get((compound_id, target_name)) or {}).items(), key=lambda kv: kv[1])[0]
                        if pair_activity_types.get((compound_id, target_name))
                        else None
                    ),
                    "best_ic50_nm": (
                        min(pair_values_nm.get((compound_id, target_name), []))
                        if pair_values_nm.get((compound_id, target_name))
                        else None
                    ),
                },
            )
            if edge:
                created += 1

    return created


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=100000)
    args = ap.parse_args()

    n = extract(limit=args.limit)
    print(f"Created/updated {n} compound->target edges.")


if __name__ == "__main__":
    main()


