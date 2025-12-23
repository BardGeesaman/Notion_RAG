"""Gene symbol -> Feature.id mapping for CRISPR results."""

from __future__ import annotations

from typing import Dict, List
from uuid import UUID

from amprenta_rag.database.models import Feature


def _normalize(s: str) -> str:
    return (s or "").strip().lower()


def map_genes_to_features(gene_symbols: List[str], db) -> Dict[str, UUID]:
    """Map gene symbols to Feature IDs (feature_type='gene').

    Lookup order:
    - Feature.name exact match
    - Feature.normalized_name match against lowercased symbol
    """
    symbols = [s for s in (gene_symbols or []) if s]
    if not symbols:
        return {}

    rows = (
        db.query(Feature)
        .filter(Feature.feature_type == "gene")
        .filter(Feature.name.in_(symbols))
        .all()
    )
    out: Dict[str, UUID] = {r.name: r.id for r in rows if r and r.name}

    remaining = [s for s in symbols if s not in out]
    if not remaining:
        return out

    norms = [_normalize(s) for s in remaining]
    rows2 = (
        db.query(Feature)
        .filter(Feature.feature_type == "gene")
        .filter(Feature.normalized_name.in_(norms))
        .all()
    )
    inv = {r.normalized_name: r.id for r in rows2 if r and r.normalized_name}
    for s in remaining:
        nid = inv.get(_normalize(s))
        if nid:
            out[s] = nid

    return out


__all__ = ["map_genes_to_features"]


