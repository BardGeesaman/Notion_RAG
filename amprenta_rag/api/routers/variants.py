"""Variant interpretation API endpoints."""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field
from sqlalchemy.orm import aliased

from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
from amprenta_rag.database.models import GeneBurden, Variant, VariantAnnotation, VariantSet
from amprenta_rag.database.session import db_session
from amprenta_rag.variants.ingestion_service import ingest_vep_tsv


router = APIRouter(prefix="/variants", tags=["Variants"])


class IngestVariantSetRequest(BaseModel):
    name: str
    vep_tsv_path: str


class VariantSetResponse(BaseModel):
    id: UUID
    name: str
    description: Optional[str] = None
    source_file: Optional[str] = None
    source_type: Optional[str] = None
    n_variants: Optional[int] = None
    n_genes: Optional[int] = None
    status: Optional[str] = None

    class Config:
        from_attributes = True


class VariantResponse(BaseModel):
    id: UUID
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    rs_id: Optional[str] = None
    gene_symbol: Optional[str] = None
    consequence: Optional[str] = None
    impact: Optional[str] = None
    gnomad_af: Optional[float] = None

    clinical_significance: Optional[str] = None
    review_status: Optional[str] = None
    clinvar_id: Optional[str] = None
    condition: Optional[str] = None


class GeneBurdenResponse(BaseModel):
    id: UUID
    gene_symbol: Optional[str] = None
    feature_id: Optional[UUID] = None
    n_variants: Optional[int] = None
    n_pathogenic: Optional[int] = None
    n_vus: Optional[int] = None
    n_benign: Optional[int] = None
    burden_score: Optional[float] = None

    class Config:
        from_attributes = True


class PathwayEnrichRequest(BaseModel):
    min_burden_score: float = Field(5.0, ge=0.0)


@router.post("/sets", response_model=VariantSetResponse)
def ingest_variant_set(payload: IngestVariantSetRequest) -> VariantSetResponse:
    with db_session() as db:
        try:
            vs = ingest_vep_tsv(payload.vep_tsv_path, payload.name, db)
        except Exception as e:  # noqa: BLE001
            raise HTTPException(status_code=400, detail=str(e))
        return VariantSetResponse.model_validate(vs)


@router.get("/sets", response_model=List[VariantSetResponse])
def list_variant_sets(limit: int = 100) -> List[VariantSetResponse]:
    with db_session() as db:
        rows = db.query(VariantSet).order_by(VariantSet.created_at.desc()).limit(limit).all()
        return [VariantSetResponse.model_validate(r) for r in rows]


@router.get("/sets/{variant_set_id}", response_model=VariantSetResponse)
def get_variant_set(variant_set_id: UUID) -> VariantSetResponse:
    with db_session() as db:
        vs = db.query(VariantSet).filter(VariantSet.id == variant_set_id).first()
        if not vs:
            raise HTTPException(status_code=404, detail="VariantSet not found")
        return VariantSetResponse.model_validate(vs)


@router.get("/sets/{variant_set_id}/variants", response_model=List[VariantResponse])
def list_variants(
    variant_set_id: UUID,
    gene: Optional[str] = None,
    impact: Optional[str] = None,
    significance: Optional[str] = None,
    limit: int = Query(100, ge=1, le=1000),
) -> List[VariantResponse]:
    with db_session() as db:
        ann = aliased(VariantAnnotation)
        q = (
            db.query(Variant, ann)
            .outerjoin(ann, (ann.variant_id == Variant.id) & (ann.annotation_source == "clinvar"))
            .filter(Variant.variant_set_id == variant_set_id)
        )
        if gene:
            q = q.filter(Variant.gene_symbol == gene)
        if impact:
            q = q.filter(Variant.impact == impact)
        if significance:
            q = q.filter(ann.clinical_significance.ilike(f"%{significance}%"))

        rows = q.order_by(Variant.chromosome.asc(), Variant.position.asc()).limit(limit).all()
        out: List[VariantResponse] = []
        for v, a in rows:
            out.append(
                VariantResponse(
                    id=v.id,
                    chromosome=v.chromosome,
                    position=v.position,
                    ref_allele=v.ref_allele,
                    alt_allele=v.alt_allele,
                    rs_id=v.rs_id,
                    gene_symbol=v.gene_symbol,
                    consequence=v.consequence,
                    impact=v.impact,
                    gnomad_af=v.gnomad_af,
                    clinical_significance=getattr(a, "clinical_significance", None) if a else None,
                    review_status=getattr(a, "review_status", None) if a else None,
                    clinvar_id=getattr(a, "clinvar_id", None) if a else None,
                    condition=getattr(a, "condition", None) if a else None,
                )
            )
        return out


@router.get("/sets/{variant_set_id}/burdens", response_model=List[GeneBurdenResponse])
def list_gene_burdens(
    variant_set_id: UUID,
    min_score: float = 5.0,
    limit: int = Query(50, ge=1, le=1000),
) -> List[GeneBurdenResponse]:
    with db_session() as db:
        rows = (
            db.query(GeneBurden)
            .filter(GeneBurden.variant_set_id == variant_set_id)
            .filter((GeneBurden.burden_score.is_(None)) | (GeneBurden.burden_score >= float(min_score)))
            .order_by(GeneBurden.burden_score.desc().nullslast())
            .limit(limit)
            .all()
        )
        return [GeneBurdenResponse.model_validate(r) for r in rows]


@router.post("/sets/{variant_set_id}/pathway-enrichment")
def pathway_enrichment(variant_set_id: UUID, payload: PathwayEnrichRequest) -> List[Dict[str, Any]]:
    with db_session() as db:
        rows = (
            db.query(GeneBurden)
            .filter(GeneBurden.variant_set_id == variant_set_id)
            .filter(GeneBurden.burden_score >= float(payload.min_burden_score))
            .all()
        )
        genes = {r.gene_symbol for r in rows if r.gene_symbol}
        if not genes:
            return []
        results = perform_pathway_enrichment(genes, {"gene"}, p_value_threshold=0.05)
        out: List[Dict[str, Any]] = []
        for r in results[:200]:
            out.append(
                {
                    "pathway_id": r.pathway.pathway_id,
                    "name": r.pathway.name,
                    "source": r.pathway.source,
                    "adjusted_p_value": r.adjusted_p_value,
                    "enrichment_ratio": r.enrichment_ratio,
                    "matched_features": r.matched_features,
                }
            )
        return out


__all__ = ["router"]


