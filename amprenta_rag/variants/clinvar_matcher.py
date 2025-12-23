"""Match internal Variant rows to ClinVar lookups and create VariantAnnotation records."""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

from amprenta_rag.database.models import Variant, VariantAnnotation


def _norm_chr(ch: str | None) -> str | None:
    if not ch:
        return None
    return str(ch).strip().replace("chr", "")


def _norm_allele(a: str | None) -> str | None:
    if not a:
        return None
    return str(a).strip().upper()


def match_variants_to_clinvar(
    variants: List[Variant],
    clinvar_lookup: Dict[Tuple[str, int, str, str], Dict[str, Any]],
    *,
    annotation_source: str = "clinvar",
) -> List[VariantAnnotation]:
    """Create VariantAnnotation records for variants that match the ClinVar lookup."""
    out: List[VariantAnnotation] = []
    for v in variants or []:
        chrom = _norm_chr(getattr(v, "chromosome", None))
        pos = getattr(v, "position", None)
        ref = _norm_allele(getattr(v, "ref_allele", None))
        alt = _norm_allele(getattr(v, "alt_allele", None))
        if chrom is None or pos is None or ref is None or alt is None:
            continue
        key = (chrom, int(pos), ref, alt)
        hit = clinvar_lookup.get(key)
        if not hit:
            continue
        out.append(
            VariantAnnotation(
                variant_id=v.id,
                annotation_source=str(annotation_source),
                clinvar_id=hit.get("clinvar_id"),
                clinical_significance=hit.get("clinical_significance"),
                review_status=hit.get("review_status"),
                condition=hit.get("condition"),
                cadd_phred=None,
                revel_score=None,
                sift_prediction=None,
                polyphen_prediction=None,
            )
        )
    return out


__all__ = ["match_variants_to_clinvar"]


