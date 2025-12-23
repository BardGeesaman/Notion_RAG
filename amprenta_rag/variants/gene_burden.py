"""Gene burden computation for variant sets."""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List
from uuid import UUID

from amprenta_rag.database.models import Feature, GeneBurden, Variant, VariantAnnotation


def _normalize(s: str) -> str:
    return (s or "").strip().lower()


def _map_genes_to_features(gene_symbols: List[str], db) -> Dict[str, UUID]:
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
        fid = inv.get(_normalize(s))
        if fid:
            out[s] = fid
    return out


def _classify_sig(sig: str | None) -> str | None:
    """Return pathogenic|vus|benign|other (best-effort)."""
    if not sig:
        return None
    s = sig.strip().lower()
    if "pathogenic" in s:
        return "pathogenic"
    if "uncertain" in s or "vus" in s:
        return "vus"
    if "benign" in s:
        return "benign"
    return "other"


def compute_gene_burden(variant_set_id: UUID, db) -> List[GeneBurden]:
    """Compute and persist per-gene burden metrics for a VariantSet.

    Counts:
    - n_variants: number of variants for gene_symbol
    - n_pathogenic: variants with ClinVar clinical_significance containing "pathogenic"
    - n_vus: variants with "uncertain significance" (VUS)
    - n_benign: variants with "benign"

    burden_score = (n_pathogenic * 10) + n_vus
    """
    vars_ = db.query(Variant).filter(Variant.variant_set_id == variant_set_id).all()
    if not vars_:
        db.query(GeneBurden).filter(GeneBurden.variant_set_id == variant_set_id).delete()
        db.commit()
        return []

    var_ids = [v.id for v in vars_]
    anns = db.query(VariantAnnotation).filter(VariantAnnotation.variant_id.in_(var_ids)).all()
    by_variant: Dict[UUID, List[VariantAnnotation]] = defaultdict(list)
    for a in anns:
        by_variant[a.variant_id].append(a)

    # Clear existing
    db.query(GeneBurden).filter(GeneBurden.variant_set_id == variant_set_id).delete()

    # Group by gene_symbol
    grouped: Dict[str, List[Variant]] = defaultdict(list)
    for v in vars_:
        g = (v.gene_symbol or "").strip()
        if not g:
            continue
        grouped[g].append(v)

    gene_to_feature = _map_genes_to_features(list(grouped.keys()), db)

    out: List[GeneBurden] = []
    for gene, vs in grouped.items():
        n_variants = len(vs)
        n_path, n_vus, n_ben, = 0, 0, 0
        for v in vs:
            sig = None
            # Prefer ClinVar annotation if present
            if by_variant.get(v.id):
                sig = by_variant[v.id][0].clinical_significance
            cls = _classify_sig(sig)
            if cls == "pathogenic":
                n_path += 1
            elif cls == "vus":
                n_vus += 1
            elif cls == "benign":
                n_ben += 1

        burden_score = (n_path * 10) + n_vus
        gb = GeneBurden(
            variant_set_id=variant_set_id,
            gene_symbol=gene,
            feature_id=gene_to_feature.get(gene),
            n_variants=n_variants,
            n_pathogenic=n_path,
            n_vus=n_vus,
            n_benign=n_ben,
            burden_score=float(burden_score),
        )
        db.add(gb)
        out.append(gb)

    db.commit()
    return out


__all__ = ["compute_gene_burden"]


