"""Extract target(feature) -> pathway edges from feature_pathway_map.

This uses the existing `feature_pathway_map` association table which maps
Feature UUIDs to Pathway UUIDs (pathway nodes may be backed by external sources).

Edges created:
  feature(UUID) --[in_pathway]--> pathway(UUID)

NOTE on ID Strategy:
    Both feature_id and pathway_id in feature_pathway_map are proper UUIDs
    referencing the Feature and Pathway tables respectively. No ID conversion
    is needed here since both ends are already normalized entities.
"""

from __future__ import annotations

import argparse
from typing import Any, Dict
from uuid import UUID

from amprenta_rag.database.models import Feature, feature_pathway_map
from amprenta_rag.database.session import db_session
from amprenta_rag.graph.edge_builder import EdgeBuilder


def extract(limit: int = 200000) -> int:
    builder = EdgeBuilder()
    created = 0

    with db_session() as db:
        rows = (
            db.query(
                feature_pathway_map.c.feature_id,
                feature_pathway_map.c.pathway_id,
                feature_pathway_map.c.meta,
            )
            .limit(limit)
            .all()
        )

        for feature_id, pathway_id, meta in rows:
            if not feature_id or not pathway_id:
                continue
            # Ensure feature exists
            feat = db.query(Feature).filter(Feature.id == feature_id).first()
            if feat is None:
                continue

            provenance: Dict[str, Any] = {"meta": meta} if meta is not None else {}
            builder.create_edge(
                source_entity_type="feature",
                source_entity_id=UUID(str(feature_id)),
                target_entity_type="pathway",
                target_entity_id=UUID(str(pathway_id)),
                relationship_type="in_pathway",
                confidence=0.8,
                evidence_source="feature_pathway_map",
                provenance=provenance,
            )
            created += 1

    return created


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=200000)
    args = ap.parse_args()

    n = extract(limit=args.limit)
    print(f"Created/updated {n} target->pathway edges.")


if __name__ == "__main__":
    main()


