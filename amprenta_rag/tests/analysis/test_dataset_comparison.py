from __future__ import annotations

import pytest

from amprenta_rag.analysis import dataset_comparison as dc
from amprenta_rag.analysis.dataset_comparison import (
    DatasetCluster,
    DatasetComparison,
    compute_jaccard_similarity,
)


def test_compute_jaccard_similarity_basic() -> None:
    assert compute_jaccard_similarity(set(), set()) == 1.0
    assert compute_jaccard_similarity({"a"}, set()) == 0.0
    assert compute_jaccard_similarity({"a", "b"}, {"b", "c"}) == pytest.approx(1 / 3)


def test_compare_datasets_happy_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(dc, "_get_dataset_name", lambda ds_id: f"Name-{ds_id}")

    features = {
        "d1": {"gene": {"A", "B"}},
        "d2": {"gene": {"A", "C"}},
    }

    monkeypatch.setattr(
        dc,
        "extract_dataset_features_by_type",
        lambda ds_id, use_cache=True: features[ds_id],
    )

    result = dc.compare_datasets("d1", "d2", use_cache=False)

    assert result.dataset1_name == "Name-d1"
    assert result.dataset2_name == "Name-d2"
    assert result.shared_features == {"gene": {"A"}}
    assert result.unique_to_dataset1 == {"gene": {"B"}}
    assert result.unique_to_dataset2 == {"gene": {"C"}}
    assert result.similarity_by_omics["gene"] == pytest.approx(1 / 3)
    assert result.overall_similarity == pytest.approx(1 / 3)
    assert result.jaccard_similarity == pytest.approx(1 / 3)


def test_compare_multiple_datasets_pairs(monkeypatch: pytest.MonkeyPatch) -> None:
    comparisons_created: list[tuple[str, str]] = []

    def _fake_compare(d1: str, d2: str, use_cache: bool = True) -> DatasetComparison:
        comparisons_created.append((d1, d2))
        return DatasetComparison(
            dataset1_id=d1,
            dataset2_id=d2,
            dataset1_name=f"Name-{d1}",
            dataset2_name=f"Name-{d2}",
            overall_similarity=0.5,
            similarity_by_omics={},
            shared_features={},
            unique_to_dataset1={},
            unique_to_dataset2={},
            jaccard_similarity=0.0,
        )

    monkeypatch.setattr(dc, "compare_datasets", _fake_compare)

    results = dc.compare_multiple_datasets(["a", "b", "c"], use_cache=False)

    assert len(results) == 3  # (a,b), (a,c), (b,c)
    assert comparisons_created == [("a", "b"), ("a", "c"), ("b", "c")]


def test_cluster_datasets_by_similarity() -> None:
    comps = [
        DatasetComparison(
            dataset1_id="d1",
            dataset2_id="d2",
            dataset1_name="D1",
            dataset2_name="D2",
            overall_similarity=0.8,
            similarity_by_omics={},
            shared_features={},
            unique_to_dataset1={},
            unique_to_dataset2={},
            jaccard_similarity=0.0,
        ),
        DatasetComparison(
            dataset1_id="d2",
            dataset2_id="d3",
            dataset1_name="D2",
            dataset2_name="D3",
            overall_similarity=0.7,
            similarity_by_omics={},
            shared_features={},
            unique_to_dataset1={},
            unique_to_dataset2={},
            jaccard_similarity=0.0,
        ),
        DatasetComparison(
            dataset1_id="d1",
            dataset2_id="d3",
            dataset1_name="D1",
            dataset2_name="D3",
            overall_similarity=0.2,
            similarity_by_omics={},
            shared_features={},
            unique_to_dataset1={},
            unique_to_dataset2={},
            jaccard_similarity=0.0,
        ),
    ]

    clusters = dc.cluster_datasets_by_similarity(comps, similarity_threshold=0.5)

    assert len(clusters) == 1
    cluster = clusters[0]
    assert set(cluster.dataset_ids) == {"d1", "d2", "d3"}
    # Average includes all pair similarities
    assert cluster.average_similarity == pytest.approx((0.8 + 0.7 + 0.2) / 3)
    assert cluster.representative_dataset in {"d1", "d2", "d3"}


def test_generate_comparison_report_contains_sections() -> None:
    comparison = DatasetComparison(
        dataset1_id="d1",
        dataset2_id="d2",
        dataset1_name="Dataset One",
        dataset2_name="Dataset Two",
        overall_similarity=0.5,
        similarity_by_omics={"gene": 0.4},
        shared_features={"gene": {"A", "B"}},
        unique_to_dataset1={"gene": {"C"}},
        unique_to_dataset2={"gene": {"D"}},
        jaccard_similarity=0.3,
    )

    report = dc.generate_comparison_report(comparison)

    assert "Dataset One vs Dataset Two" in report
    assert "Overall Similarity" in report
    assert "Shared Features" in report
    assert "Unique to Dataset One" in report
    assert "Unique to Dataset Two" in report


def test_generate_clustering_report_lists_clusters() -> None:
    clusters = [
        DatasetCluster(
            cluster_id=0,
            dataset_ids=["d1", "d2"],
            dataset_names=["One", "Two"],
            average_similarity=0.6,
            representative_dataset="d1",
        )
    ]

    report = dc.generate_clustering_report(clusters)

    assert "Cluster 1" in report
    assert "Number of datasets" in report
    assert "One (`d1`)" in report

