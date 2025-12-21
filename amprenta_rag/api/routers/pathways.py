"""
Pathway analysis API routes.

These endpoints provide pathway enrichment for datasets and cross-omics pathway
analysis for programs, plus convenience helpers for mapping features to a
pathway context.
"""

from __future__ import annotations

from typing import List, Set
from uuid import UUID

from fastapi import APIRouter, HTTPException

from amprenta_rag.analysis.cross_omics_pathways import (
    combine_omics_features,
    get_cross_omics_enrichment,
)
from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
from amprenta_rag.api import schemas
from amprenta_rag.database.models import Feature, dataset_feature_assoc
from amprenta_rag.database.session import db_session

router = APIRouter()


def _load_dataset_features(dataset_id: UUID) -> tuple[Set[str], Set[str]]:
    """
    Load a dataset's features and derive a canonical set of feature type tokens.

    Returns:
        Tuple of (feature_names, type_tokens) where type_tokens are normalized to
        {"gene","protein","metabolite","lipid"}.
    """
    with db_session() as db:
        feats: List[Feature] = (
            db.query(Feature)
            .join(dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id)
            .filter(dataset_feature_assoc.c.dataset_id == dataset_id)
            .all()
        )
    features: Set[str] = set()
    types: Set[str] = set()
    for f in feats:
        name = getattr(f, "name", None)
        ftype = getattr(f, "feature_type", None) or getattr(f, "type", None)
        if name:
            features.add(name)
        if ftype:
            # Map to canonical type tokens
            if ftype in ("gene", "rna", "transcript"):
                types.add("gene")
            elif ftype in ("protein", "proteomics"):
                types.add("protein")
            elif ftype in ("metabolite", "metabolomics"):
                types.add("metabolite")
            elif ftype in ("lipid", "lipidomics"):
                types.add("lipid")
    return features, types


@router.post(
    "/pathways/enrich",
    summary="Run pathway enrichment for a dataset",
    response_model=List[schemas.ConvergentPathway],
)
def pathway_enrich(dataset_id: UUID) -> List[schemas.ConvergentPathway]:
    """Run pathway enrichment for a dataset and return convergent pathway results."""
    features, types = _load_dataset_features(dataset_id)
    if not features:
        raise HTTPException(status_code=404, detail="No features for dataset")
    results = perform_pathway_enrichment(features, types, p_value_threshold=0.05)
    # Convert to convergent pathways with simple matched_by_omics from type mapping
    conv: List[schemas.ConvergentPathway] = []
    for r in results:
        conv.append(
            schemas.ConvergentPathway(
                pathway_id=r.pathway.pathway_id,
                name=r.pathway.name,
                source=r.pathway.source,
                adjusted_p_value=r.adjusted_p_value,
                enrichment_ratio=r.enrichment_ratio,
                matched_by_omics={"all": r.matched_features},
            )
        )
    return conv


@router.get(
    "/programs/{program_id}/pathway-analysis",
    summary="Cross-omics pathway analysis for program",
    response_model=schemas.CrossOmicsEnrichmentResult,
)
def program_pathway_analysis(program_id: UUID) -> schemas.CrossOmicsEnrichmentResult:
    """Run cross-omics pathway analysis for a program."""
    result = get_cross_omics_enrichment(program_id)
    return schemas.CrossOmicsEnrichmentResult.model_validate(result.asdict())


@router.get(
    "/pathways/{pathway_id}/features",
    summary="Features mapped to a pathway by omics type",
    response_model=schemas.PathwayFeatures,
)
def pathway_features(pathway_id: str, dataset_ids: List[UUID]) -> schemas.PathwayFeatures:
    """Return features grouped by omics type for the provided datasets (placeholder mapping)."""
    # Combine features across provided datasets
    combined = combine_omics_features(dataset_ids)
    # Simple filter: return all combined features grouped by omics for now
    return schemas.PathwayFeatures(
        pathway_id=pathway_id,
        pathway_name=pathway_id,
        source="mixed",
        features_by_omics=combined.asdict(),
    )

