"""Extract signature -> pathway edges using pathway enrichment.

We compute pathway enrichment on a signature's feature names/types and create edges:
  signature(UUID) --[enriched_in_pathway]--> pathway(UUID)

Because pathways are not currently persisted as a first-class DB model, we create
stable pathway UUIDs from the pathway (source, pathway_id) tuple using uuid5.

NOTE on Pathway ID Strategy:
    Pathways from enrichment analysis are external resources (KEGG, Reactome, etc.)
    without persistent DB records. We generate deterministic UUIDs using:
        uuid5(NAMESPACE, f"{source}:{pathway_id}")
    This ensures the same pathway always gets the same UUID across runs, enabling
    proper edge deduplication. Pathway metadata (name, source, ID) is preserved
    in the edge's provenance field for downstream queries.
"""

from __future__ import annotations

import argparse
from typing import Iterable, Optional, Set
from uuid import UUID, uuid5

from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
from amprenta_rag.database.models import Feature, Signature, SignatureComponent, signature_feature_assoc
from amprenta_rag.database.session import db_session
from amprenta_rag.graph.edge_builder import EdgeBuilder


PATHWAY_UUID_NAMESPACE = UUID("00000000-0000-0000-0000-000000000042")


def _stable_pathway_uuid(source: str, pathway_id: str) -> UUID:
    return uuid5(PATHWAY_UUID_NAMESPACE, f"{source}:{pathway_id}")


def _signature_feature_names(db, signature_id: UUID) -> tuple[Set[str], Set[str]]:
    names: Set[str] = set()
    types: Set[str] = set()

    # From explicit components (includes feature_type + possibly feature_name).
    comps = (
        db.query(SignatureComponent)
        .filter(SignatureComponent.signature_id == signature_id)
        .all()
    )
    for c in comps:
        if c.feature_type:
            types.add(str(c.feature_type).lower())
        if c.feature_name:
            names.add(str(c.feature_name).strip())
        elif c.feature is not None and getattr(c.feature, "name", None):
            names.add(str(c.feature.name).strip())

    # From signature_feature_assoc (Feature rows)
    feats = (
        db.query(Feature)
        .join(signature_feature_assoc, Feature.id == signature_feature_assoc.c.feature_id)
        .filter(signature_feature_assoc.c.signature_id == signature_id)
        .all()
    )
    for f in feats:
        ft = getattr(f, "feature_type", None)
        if ft:
            types.add(str(ft).lower())
        nm = getattr(f, "normalized_name", None) or getattr(f, "name", None)
        if nm:
            names.add(str(nm).strip())

    names = {n for n in names if n}
    # Default type if missing (enrichment code currently applies input_features for each type anyway)
    if not types:
        types = {"gene"}
    return names, types


def extract(
    limit: int = 200,
    pathway_sources: Optional[Iterable[str]] = None,
    p_value_threshold: float = 0.05,
    max_edges_per_signature: int = 30,
    signature_id: Optional[UUID] = None,
) -> int:
    builder = EdgeBuilder()
    created = 0
    sources = list(pathway_sources) if pathway_sources else None

    with db_session() as db:
        q = db.query(Signature).order_by(Signature.created_at.desc())
        if signature_id is not None:
            q = q.filter(Signature.id == signature_id)
        sigs = q.limit(limit).all()

        for sig in sigs:
            names, types = _signature_feature_names(db, sig.id)
            if not names:
                continue
            try:
                results = perform_pathway_enrichment(
                    input_features=names,
                    input_feature_types=types,
                    pathway_sources=sources,
                    p_value_threshold=p_value_threshold,
                )
            except Exception:
                # If enrichment fails (network, scipy, etc.), skip signature rather than crashing whole run.
                continue

            for r in results[: max_edges_per_signature]:
                p = r.pathway
                pathway_uuid = _stable_pathway_uuid(p.source, p.pathway_id)
                confidence = max(0.0, min(1.0, 1.0 - float(r.adjusted_p_value)))
                builder.create_edge(
                    source_entity_type="signature",
                    source_entity_id=sig.id,
                    target_entity_type="pathway",
                    target_entity_id=pathway_uuid,
                    relationship_type="enriched_in_pathway",
                    confidence=confidence,
                    evidence_source="pathway_enrichment",
                    provenance={
                        "pathway_source": p.source,
                        "pathway_id": p.pathway_id,
                        "pathway_name": p.name,
                        "p_value": r.p_value,
                        "adjusted_p_value": r.adjusted_p_value,
                        "enrichment_ratio": r.enrichment_ratio,
                        "matched_features": r.matched_features[:200],
                        "signature_name": sig.name,
                    },
                )
                created += 1

    return created


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=200)
    ap.add_argument("--signature-id", type=str, default=None)
    ap.add_argument("--pathway-sources", type=str, default="KEGG,Reactome")
    ap.add_argument("--p-threshold", type=float, default=0.05)
    ap.add_argument("--max-edges-per-signature", type=int, default=30)
    args = ap.parse_args()

    sig_id = UUID(args.signature_id) if args.signature_id else None
    sources = [s.strip() for s in args.pathway_sources.split(",") if s.strip()]
    n = extract(
        limit=args.limit,
        pathway_sources=sources,
        p_value_threshold=args.p_threshold,
        max_edges_per_signature=args.max_edges_per_signature,
        signature_id=sig_id,
    )
    print(f"Created/updated {n} signature->pathway edges.")


if __name__ == "__main__":
    main()


