from typing import List
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, Signature
from amprenta_rag.domain.analytics import SignatureValidationMetrics, SignatureValidationResult
from amprenta_rag.ingestion.signature_matching.matching import score_signature_against_dataset


def validate_signature_against_all_datasets(
    signature_id: UUID, coverage_threshold: float = 0.1
) -> SignatureValidationResult:
    with db_session() as db:
        db.query(Signature).filter(Signature.id == signature_id).first()
        datasets = db.query(Dataset).all()
        matched_ids = []
        scores = []
        for ds in datasets:
            res = score_signature_against_dataset(signature_id, ds.id, db=db)
            if res and (res.get("overlap", 0) >= coverage_threshold):
                matched_ids.append(ds.id)
                scores.append(res.get("score", 0))
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
        return [validate_signature_against_all_datasets(s.id, coverage_threshold) for s in sigs]
