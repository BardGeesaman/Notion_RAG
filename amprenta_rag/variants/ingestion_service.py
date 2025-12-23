"""Variant set ingestion service (VEP TSV + ClinVar matching + gene burden)."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional
from uuid import UUID

from amprenta_rag.database.models import Feature, Variant, VariantAnnotation, VariantSet
from amprenta_rag.variants.clinvar_loader import download_clinvar, parse_clinvar
from amprenta_rag.variants.clinvar_matcher import match_variants_to_clinvar
from amprenta_rag.variants.gene_burden import compute_gene_burden
from amprenta_rag.variants.vep_parser import parse_vep_tsv


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


def _infer_ref_alt(raw: dict, alt_from_allele: Optional[str]) -> tuple[str, str]:
    """Infer ref/alt from common VEP columns; fallback to alt_from_allele."""
    uploaded = raw.get("Uploaded_variation") or raw.get("uploaded_variation") or raw.get("Uploaded Variation")
    if uploaded:
        s = str(uploaded)
        # common: "1_123_A/T" or "1_123_A_G"
        if "_" in s and "/" in s:
            parts = s.split("_")
            if len(parts) >= 3:
                ref = parts[-1].split("/", 1)[0]
                alt = parts[-1].split("/", 1)[1]
                return ref, alt
        if "_" in s and len(s.split("_")) >= 4:
            parts = s.split("_")
            ref = parts[-2]
            alt = parts[-1]
            return ref, alt
        if ">" in s:
            # HGVS-like: g.123A>G
            after = s.split(">", 1)
            alt = after[1][-1] if after[1] else (alt_from_allele or "N")
            before = after[0]
            ref = before[-1] if before else "N"
            return ref, alt

    # fallbacks
    alt = (alt_from_allele or "N").strip()
    return "N", alt


def _norm_chr(ch: str | None) -> str | None:
    if not ch:
        return None
    return str(ch).strip().replace("chr", "")


def ingest_vep_tsv(file_path: str, name: str, db) -> VariantSet:
    """Ingest a VEP TSV and populate VariantSet + Variant + ClinVar annotations + GeneBurden."""
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    vs = VariantSet(
        name=name,
        description=None,
        source_file=str(p),
        source_type="vep_tsv",
        n_variants=None,
        n_genes=None,
        status="processing",
    )
    db.add(vs)
    db.commit()
    db.refresh(vs)

    try:
        parsed = parse_vep_tsv(str(p))
        gene_symbols = []
        for r in parsed:
            gs = r.get("gene_symbol") or r.get("gene")
            if gs:
                gene_symbols.append(str(gs))
        gene_to_feature = _map_genes_to_features(list(set(gene_symbols)), db)

        variants: List[Variant] = []
        for r in parsed:
            chrom = _norm_chr(r.get("chromosome"))
            pos = r.get("position")
            if chrom is None or pos is None:
                continue
            raw = r.get("raw") or {}
            ref, alt = _infer_ref_alt(raw, r.get("alt_allele"))
            gene_symbol = r.get("gene_symbol") or r.get("gene")
            gene_symbol = str(gene_symbol).strip() if gene_symbol else None
            variants.append(
                Variant(
                    variant_set_id=vs.id,
                    chromosome=str(chrom),
                    position=int(pos),
                    ref_allele=str(ref),
                    alt_allele=str(alt),
                    rs_id=(raw.get("Existing_variation") or raw.get("rsid") or raw.get("RSID")),
                    hgvs_genomic=(raw.get("HGVSg") or raw.get("HGVSG")),
                    hgvs_coding=(raw.get("HGVSc") or raw.get("HGVSC")),
                    hgvs_protein=(raw.get("HGVSp") or raw.get("HGVSP")),
                    gene_symbol=gene_symbol,
                    feature_id=gene_to_feature.get(gene_symbol) if gene_symbol else None,
                    consequence=r.get("consequence"),
                    impact=r.get("impact"),
                    gnomad_af=r.get("gnomad_af"),
                )
            )

        if variants:
            db.bulk_save_objects(variants)
        db.commit()

        # Refresh variants with IDs for matching
        db_vars = db.query(Variant).filter(Variant.variant_set_id == vs.id).all()

        # ClinVar: download/cache + parse
        clin_dir = Path("data") / "clinvar"
        clin_dir.mkdir(parents=True, exist_ok=True)
        clin_txt = clin_dir / "variant_summary.txt"
        if not clin_txt.exists():
            download_clinvar(str(clin_dir))
        clin_lookup = parse_clinvar(str(clin_txt))

        annotations: List[VariantAnnotation] = match_variants_to_clinvar(db_vars, clin_lookup)
        if annotations:
            db.bulk_save_objects(annotations)
        db.commit()

        # Compute gene burden (persists)
        compute_gene_burden(vs.id, db)

        # Update counts + status
        vs.n_variants = len(db_vars)
        vs.n_genes = len({v.gene_symbol for v in db_vars if v.gene_symbol})
        vs.status = "completed"
        db.add(vs)
        db.commit()
        db.refresh(vs)
        return vs
    except Exception:
        vs.status = "failed"
        db.add(vs)
        db.commit()
        raise


__all__ = ["ingest_vep_tsv"]


