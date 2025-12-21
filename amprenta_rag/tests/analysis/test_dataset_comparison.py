from __future__ import annotations

import pytest

from amprenta_rag.analysis import dataset_comparison as dc
from amprenta_rag.analysis.dataset_comparison import DatasetCluster, DatasetComparison


def test_compute_jaccard_similarity_basic():
    assert dc.compute_jaccard_similarity(set(), set()) == 1.0
    assert dc.compute_jaccard_similarity({"a"}, {"a"}) == 1.0
    assert dc.compute_jaccard_similarity({"a"}, {"b"}) == 0.0
    assert dc.compute_jaccard_similarity({"a", "b"}, {"b", "c"}) == pytest.approx(1 / 3)


def test_compare_datasets_happy(monkeypatch):
    monkeypatch.setattr(dc, "_get_dataset_name", lambda did: f"name-{did}")

    def fake_extract(dataset_id, use_cache=True):
        if dataset_id == "d1":
            return {"gene": {"a", "b"}, "protein": {"p1"}}
        return {"gene": {"b", "c"}, "metabolite": {"m1"}}

    monkeypatch.setattr(dc, "extract_dataset_features_by_type", fake_extract)

    res = dc.compare_datasets("d1", "d2", use_cache=False)
    assert res.dataset1_name == "name-d1"
    assert res.dataset2_name == "name-d2"
    assert res.shared_features["gene"] == {"b"}
    assert pytest.approx(res.similarity_by_omics["gene"], 0.01) == 1 / 3
    assert pytest.approx(res.overall_similarity, 0.01) == 4 / 18  # weighted by feature counts
    assert pytest.approx(res.jaccard_similarity, 0.01) == 0.2


def test_compare_datasets_no_features(monkeypatch):
    monkeypatch.setattr(dc, "_get_dataset_name", lambda did: did)
    monkeypatch.setattr(dc, "extract_dataset_features_by_type", lambda did, use_cache=True: {})

    res = dc.compare_datasets("d1", "d2")
    assert res.overall_similarity == 0.0
    assert res.jaccard_similarity == 1.0


def test_compare_multiple_datasets_pairs(monkeypatch):
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
    assert all(isinstance(r, DatasetComparison) for r in results)


def test_cluster_datasets_by_similarity():
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
        DatasetComparison(
            dataset1_id="d2",
            dataset2_id="d3",
            dataset1_name="D2",
            dataset2_name="D3",
            overall_similarity=0.5,
            similarity_by_omics={},
            shared_features={},
            unique_to_dataset1={},
            unique_to_dataset2={},
            jaccard_similarity=0.0,
        ),
    ]

    clusters = dc.cluster_datasets_by_similarity(comps, similarity_threshold=0.6)

    assert len(clusters) == 2
    assert any(isinstance(c, DatasetCluster) for c in clusters)
    # One cluster should contain the high-similarity pair d1/d2
    assert any({"d1", "d2"}.issubset(set(cluster.dataset_ids)) for cluster in clusters)

