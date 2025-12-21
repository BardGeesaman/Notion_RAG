from __future__ import annotations

import types
import sys

from amprenta_rag.analysis import enrichment


def test_enrich_dataset_pathways_no_features(monkeypatch):
    monkeypatch.setitem(
        enrichment.__dict__, "extract_dataset_features_by_type", lambda *a, **k: {"rna": []}
    )
    res = enrichment.enrich_dataset_pathways("ds1")
    assert res == []


def test_enrich_dataset_pathways_calls_enrichment(monkeypatch):
    monkeypatch.setitem(
        enrichment.__dict__,
        "extract_dataset_features_by_type",
        lambda *a, **k: {"rna": ["A"], "protein": ["P"]},
    )
    called = {"yes": False}

    def fake_perform(input_features, input_feature_types, pathway_sources, p_value_threshold):
        called["yes"] = True
        return ["ok"]

    monkeypatch.setitem(enrichment.__dict__, "perform_pathway_enrichment", fake_perform)
    res = enrichment.enrich_dataset_pathways("ds1")
    assert res == ["ok"]
    assert called["yes"]


def test_enrich_signature_pathways_missing(monkeypatch):
    monkeypatch.setitem(enrichment.__dict__, "fetch_all_signatures_from_notion", lambda: [])
    res = enrichment.enrich_signature_pathways("sig1")
    assert res == []


def test_enrich_signature_pathways_loads(monkeypatch):
    comps = [type("C", (), {"feature_name": "A", "feature_type": "rna"})()]
    sig = type("S", (), {"components": comps})()

    # Provide import stubs for signature_loader module
    fake_loader = types.SimpleNamespace(
        fetch_all_signatures_from_notion=lambda: [{"id": "sig1"}],
        load_signature_from_notion_page=lambda p: sig if p.get("id") == "sig1" else None,
    )
    sys.modules["amprenta_rag.ingestion.signature_matching.signature_loader"] = fake_loader

    monkeypatch.setitem(
        enrichment.__dict__,
        "perform_pathway_enrichment",
        lambda input_features, input_feature_types, pathway_sources, p_value_threshold: ["ok"],
    )

    res = enrichment.enrich_signature_pathways("sig1")
    assert res == ["ok"]

