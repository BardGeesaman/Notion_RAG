from __future__ import annotations

from amprenta_rag.analysis import enrichment


class FakeComp:
    def __init__(self, feature_name: str, feature_type: str | None):
        self.feature_name = feature_name
        self.feature_type = feature_type


def test_enrich_dataset_pathways_no_features(monkeypatch):
    def fake_extract(dataset_id, use_cache=True):
        return {"rna": set(), "prot": set()}

    monkeypatch.setattr(enrichment, "perform_pathway_enrichment", lambda **kwargs: ["hit"])
    monkeypatch.setattr(
        "amprenta_rag.ingestion.multi_omics_scoring.extract_dataset_features_by_type",
        fake_extract,
    )

    assert enrichment.enrich_dataset_pathways("ds1") == []


def test_enrich_dataset_pathways_runs(monkeypatch):
    def fake_extract(dataset_id, use_cache=True):
        return {"rna": {"G1", "G2"}}

    monkeypatch.setattr(enrichment, "perform_pathway_enrichment", lambda **kwargs: ["hit"])
    monkeypatch.setattr(
        "amprenta_rag.ingestion.multi_omics_scoring.extract_dataset_features_by_type",
        fake_extract,
    )

    assert enrichment.enrich_dataset_pathways("ds1") == ["hit"]

