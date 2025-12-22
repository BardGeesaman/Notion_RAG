from __future__ import annotations

from typing import List, Optional, Set, Any, cast
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query

from amprenta_rag.api import schemas
from amprenta_rag.database.models import Dataset, Feature, Signature
from amprenta_rag.database.session import db_session
from amprenta_rag.signatures.signature_loader import Signature as SigDef
from amprenta_rag.signatures.signature_loader import SignatureComponent
from amprenta_rag.signatures.signature_scoring import score_signature

router = APIRouter()


@router.get(
    "/signatures/{signature_id}/explain/{dataset_id}",
    summary="Explain signature match",
    response_model=schemas.SignatureExplanationResponse,
)
def explain_signature_match(
    signature_id: UUID,
    dataset_id: UUID,
    top_n: int = Query(10, ge=1, le=50),
) -> schemas.SignatureExplanationResponse:
    """Get per-feature explanation of why a signature matches a dataset."""
    with db_session() as db:
        # Load signature
        sig_model: Optional[Signature] = (
            db.query(Signature).filter(Signature.id == signature_id).first()
        )
        if not sig_model:
            raise HTTPException(404, "Signature not found")

        # Load dataset
        dataset: Optional[Dataset] = (
            db.query(Dataset).filter(Dataset.id == dataset_id).first()
        )
        if not dataset:
            raise HTTPException(404, "Dataset not found")

        # Get dataset features via relationship
        features: List[Feature] = list(dataset.features)  # type: ignore[call-overload]
        dataset_species: Set[str] = {str(f.name) for f in features if getattr(f, "name", None)}

        # Build signature definition
        components: List[SignatureComponent] = []
        if sig_model.components:
            components_source: List[dict[str, Any]] = cast(List[dict[str, Any]], list(sig_model.components or []))
            for comp_data in components_source:
                if isinstance(comp_data, dict):
                    components.append(
                        SignatureComponent(
                            feature_name=comp_data.get("feature_name", ""),
                            feature_type=comp_data.get("feature_type"),
                            direction=comp_data.get("direction"),
                            weight=comp_data.get("weight"),
                        )
                    )

        sig_def = SigDef(name=sig_model.name or str(sig_model.id), components=components)

        # Score
        result = score_signature(sig_def, dataset_species)

        expected_dir_map = {
            c.feature_name: c.direction for c in components if c.feature_name
        }

        # Build response
        matches = list(result.component_matches)

        all_contribs = [
            schemas.FeatureContribution(
                feature_name=m.signature_species,
                matched_to=m.matched_dataset_species,
                direction_expected=expected_dir_map.get(m.signature_species),
                direction_actual=None,  # Dataset direction data not yet available
                weight=m.weight,
                contribution=m.contribution,
                match_type=m.match_type,
                direction_match=m.direction_match,
            )
            for m in matches
        ]

        top = result.get_top_contributors(top_n)

        return schemas.SignatureExplanationResponse(
            signature_id=signature_id,
            signature_name=sig_model.name or str(sig_model.id),
            dataset_id=dataset_id,
            dataset_name=dataset.name or str(dataset.id),
            total_score=result.total_score,
            direction_concordance=result.get_direction_concordance(),
            top_positive=[
                schemas.FeatureContribution(
                    feature_name=m.signature_species,
                    matched_to=m.matched_dataset_species,
                    direction_expected=expected_dir_map.get(m.signature_species),
                    direction_actual=None,
                    weight=m.weight,
                    contribution=m.contribution,
                    match_type=m.match_type,
                    direction_match=m.direction_match,
                )
                for m in top["positive"]
            ],
            top_negative=[
                schemas.FeatureContribution(
                    feature_name=m.signature_species,
                    matched_to=m.matched_dataset_species,
                    direction_expected=expected_dir_map.get(m.signature_species),
                    direction_actual=None,
                    weight=m.weight,
                    contribution=m.contribution,
                    match_type=m.match_type,
                    direction_match=m.direction_match,
                )
                for m in top["negative"]
            ],
            all_contributions=all_contribs,
        )

