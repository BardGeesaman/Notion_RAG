from typing import List, cast, Any
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, Signature
from amprenta_rag.domain.analytics import SignatureValidationMetrics, SignatureValidationResult
from typing import Any
from amprenta_rag.ingestion.signature_matching.matching import score_signature_against_dataset


def validate_signature_against_all_datasets(
    signature_id: UUID, coverage_threshold: float = 0.1
) -> SignatureValidationResult:
    with db_session() as db:
        sig = db.query(Signature).filter(Signature.id == signature_id).first()
        datasets = db.query(Dataset).all()
        matched_ids: List[UUID] = []
        scores: List[float] = []
        for ds in datasets:
            if ds.id is None or sig is None:
                continue
            res = score_signature_against_dataset(cast(Any, sig), {str(ds.id)})
            overlap = getattr(res, "overlap", 0)
            score_val = getattr(res, "score", None)
            if res and overlap is not None and overlap >= coverage_threshold:
                matched_ids.append(cast(UUID, ds.id))
                if score_val is not None:
                    scores.append(float(score_val))
        num_match = len(matched_ids)
        num_total = len(datasets)
        coverage = num_match / num_total if num_total else 0.0
        mean_score = sum(scores) / num_match if num_match else None
        metrics = SignatureValidationMetrics(
            signature_id=signature_id,
            num_matched_datasets=num_match,
            num_total_datasets=num_total,
            coverage=coverage,
            specificity=None,
            reproducibility=None,
            mean_score=mean_score,
            p_value=None,
        )
        summary = f"Matched {num_match}/{num_total} datasets (coverage: {coverage:.2f}, mean_score: {mean_score if mean_score else 0:.2f})"
        return SignatureValidationResult(
            signature_id=signature_id,
            metrics=metrics,
            matched_dataset_ids=matched_ids,
            summary=summary,
        )


def validate_all_signatures(coverage_threshold: float = 0.1) -> List[SignatureValidationResult]:
    with db_session() as db:
        sigs = db.query(Signature).all()
        return [
            validate_signature_against_all_datasets(cast(UUID, s.id), coverage_threshold)
            for s in sigs
            if s.id is not None
        ]
