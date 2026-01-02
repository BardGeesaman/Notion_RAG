"""Variant annotation service using VEP."""

from __future__ import annotations

import logging
from typing import Optional
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Variant, VariantAnnotation
from amprenta_rag.ingestion.genomics.vep_client import (
    annotate_variant as vep_annotate,
    annotate_variants_batch as vep_annotate_batch,
    get_vep_version,
    VEPAnnotation,
)

logger = logging.getLogger(__name__)


def annotate_variant_by_id(variant_id: UUID) -> Optional[VariantAnnotation]:
    """
    Annotate a single variant by ID using VEP.
    
    Args:
        variant_id: UUID of the variant to annotate
        
    Returns:
        Created VariantAnnotation record or None if failed
    """
    with db_session() as db:
        variant = db.query(Variant).filter(Variant.id == variant_id).first()
        
        if not variant:
            logger.warning(f"Variant {variant_id} not found")
            return None
        
        # Check required fields
        if not all([variant.chromosome, variant.position, variant.alt_allele]):
            logger.warning(f"Variant {variant_id} missing required fields for annotation")
            return None
        
        # Call VEP
        vep_result = vep_annotate(
            chromosome=variant.chromosome,
            position=variant.position,
            ref=variant.ref_allele or "N",
            alt=variant.alt_allele,
        )
        
        if not vep_result:
            logger.warning(f"VEP annotation failed for variant {variant_id}")
            return None
        
        # Create annotation record
        annotation = _create_annotation_from_vep(db, variant_id, vep_result)
        
        db.expunge(annotation)
        return annotation


def annotate_variants_bulk(
    variant_ids: list[UUID],
    batch_size: int = 200,
) -> dict:
    """
    Annotate multiple variants using VEP batch endpoint.
    
    Args:
        variant_ids: List of variant UUIDs to annotate
        batch_size: Number of variants per VEP request (max 200)
        
    Returns:
        Dict with counts: {annotated, failed, skipped}
    """
    results = {"annotated": 0, "failed": 0, "skipped": 0}
    
    with db_session() as db:
        # Fetch all variants
        variants = db.query(Variant).filter(Variant.id.in_(variant_ids)).all()
        
        # Filter to those with required fields
        valid_variants = []
        for v in variants:
            if all([v.chromosome, v.position, v.alt_allele]):
                valid_variants.append(v)
            else:
                results["skipped"] += 1
        
        if not valid_variants:
            return results
        
        # Process in batches
        for i in range(0, len(valid_variants), batch_size):
            batch = valid_variants[i:i + batch_size]
            
            # Prepare input for VEP
            vep_input = [
                {
                    "chromosome": v.chromosome,
                    "position": v.position,
                    "ref": v.ref_allele or "N",
                    "alt": v.alt_allele,
                }
                for v in batch
            ]
            
            # Call VEP batch
            vep_results = vep_annotate_batch(vep_input)
            
            # Create annotation records
            for variant, vep_result in zip(batch, vep_results):
                if vep_result:
                    _create_annotation_from_vep(db, variant.id, vep_result)
                    results["annotated"] += 1
                else:
                    results["failed"] += 1
        
        db.commit()
    
    logger.info(f"Bulk annotation complete: {results}")
    return results


def get_unannotated_variants(limit: int = 1000) -> list[UUID]:
    """
    Find variants that don't have VEP annotations.
    
    Args:
        limit: Maximum number of variant IDs to return
        
    Returns:
        List of variant UUIDs needing annotation
    """
    with db_session() as db:
        # Find variants with VCF fields but no VEP annotation
        subquery = db.query(VariantAnnotation.variant_id).filter(
            VariantAnnotation.annotation_source == "VEP"
        ).subquery()
        
        unannotated = db.query(Variant.id).filter(
            Variant.chromosome.isnot(None),
            Variant.position.isnot(None),
            Variant.alt_allele.isnot(None),
            ~Variant.id.in_(subquery),
        ).limit(limit).all()
        
        return [v[0] for v in unannotated]


def _create_annotation_from_vep(
    db,
    variant_id: UUID,
    vep_result: VEPAnnotation,
) -> VariantAnnotation:
    """Create VariantAnnotation record from VEP result."""
    vep_version = get_vep_version()
    
    annotation = VariantAnnotation(
        variant_id=variant_id,
        annotation_source="VEP",
        consequence=vep_result.consequence,
        impact=vep_result.impact,
        symbol=vep_result.symbol,
        gene_id=vep_result.gene_id,
        transcript_id=vep_result.transcript_id,
        biotype=vep_result.biotype,
        amino_acids=vep_result.amino_acids,
        codons=vep_result.codons,
        protein_position=vep_result.protein_position,
        sift_prediction=vep_result.sift_prediction,
        sift_score=vep_result.sift_score,
        polyphen_prediction=vep_result.polyphen_prediction,
        polyphen_score=vep_result.polyphen_score,
        clin_sig=vep_result.clin_sig,
        source_version=vep_version,
    )
    
    db.add(annotation)
    db.flush()
    
    return annotation
