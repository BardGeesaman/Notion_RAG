from datetime import datetime
from uuid import uuid4

from amprenta_rag.domain.analytics import (
    DatasetCluster,
    DatasetComparisonResult,
    DatasetFeatureProfile,
    EvidenceReport,
    EvidenceReportRequest,
    EvidenceSection,
    SignatureValidationMetrics,
    SignatureValidationResult,
)


def test_dataset_feature_profile_allows_optional():
    profile = DatasetFeatureProfile(
        dataset_id=uuid4(),
        omics_type="gene",
        diseases=None,
        matrix=None,
        features_by_type={"gene": {"A"}},
    )
    assert profile.features_by_type["gene"] == {"A"}


def test_dataset_comparison_result_sets_fields():
    comp = DatasetComparisonResult(
        dataset_id_1=uuid4(),
        dataset_id_2=uuid4(),
        jaccard_by_type={"gene": 0.5},
        overall_score=0.6,
        shared_features={"A"},
        unique_to_1={"B"},
        unique_to_2={"C"},
    )
    assert comp.overall_score == 0.6


def test_dataset_cluster_defaults():
    cluster = DatasetCluster(cluster_id=1, member_dataset_ids=[uuid4()])
    assert cluster.linkage_score is None


def test_evidence_report_round_trip():
    section = EvidenceSection(title="t", summary_text="s", supporting_datasets=[], key_features=None)
    report = EvidenceReport(
        entity_id=uuid4(),
        entity_type="dataset",
        generated_at=datetime.utcnow(),
        sections=[section],
    )
    assert report.sections[0].title == "t"


def test_signature_validation_result():
    metrics = SignatureValidationMetrics(
        signature_id=uuid4(),
        num_matched_datasets=1,
        num_total_datasets=2,
        coverage=0.5,
    )
    result = SignatureValidationResult(
        signature_id=metrics.signature_id,
        metrics=metrics,
        matched_dataset_ids=[uuid4()],
        summary="ok",
    )
    assert result.metrics.coverage == 0.5


def test_evidence_report_request_defaults():
    req = EvidenceReportRequest(entity_type="dataset", id=uuid4())
    assert req.include_sections is None

