from __future__ import annotations

from amprenta_rag.query import trust_scoring as ts


def test_get_source_trust_variations() -> None:
    internal_validated = ts.get_source_trust("Dataset", {"internal": True, "validated": True})
    curated = ts.get_source_trust("geo", {"curated": True})
    literature = ts.get_source_trust("literature", {})
    general = ts.get_source_trust("note", {"quality_score": 0.1})

    assert internal_validated == 1.0
    assert curated > literature
    assert general < literature


def test_get_source_trust_defaults() -> None:
    assert ts.get_source_trust("unknown") >= 0.0
    assert ts.get_source_trust("", None) >= 0.0


def test_weight_results_by_trust_orders_by_final_score() -> None:
    matches = [
        {"score": 0.5, "metadata": {"source_type": "literature"}},
        {"score": 0.4, "metadata": {"source_type": "Dataset", "validated": True, "internal": True}},
    ]

    weighted = ts.weight_results_by_trust(matches, alpha=0.5)

    assert weighted[0]["final_score"] >= weighted[1]["final_score"]
    assert "trust_score" in weighted[0]


def test_weight_results_by_trust_empty() -> None:
    assert ts.weight_results_by_trust([]) == []


def test_get_trust_summary_counts_levels() -> None:
    matches = [
        {"score": 0.5, "metadata": {"source_type": "Dataset", "validated": True}},
        {"score": 0.2, "metadata": {"source_type": "general"}},
    ]

    summary = ts.get_trust_summary(matches)

    assert summary["levels"]["internal_validated"] >= 1
    assert summary["high_trust_count"] >= 1
    assert summary["average_trust"] > 0


def test_get_trust_summary_empty() -> None:
    summary = ts.get_trust_summary([])
    assert summary["levels"] == {}
    assert summary["average_trust"] == 0.0

