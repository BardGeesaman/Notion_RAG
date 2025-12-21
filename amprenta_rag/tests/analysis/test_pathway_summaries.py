from __future__ import annotations

from amprenta_rag.analysis import pathway_summaries as ps


class FakePathway:
    def __init__(self, pathway_id: str, name: str, source: str, description: str = "", feature_types=None):
        self.pathway_id = pathway_id
        self.name = name
        self.source = source
        self.description = description
        self.feature_types = feature_types or []


class FakeResult:
    def __init__(
        self,
        pathway: FakePathway,
        p_value: float,
        adjusted_p_value: float,
        enrichment_ratio: float,
        input_features: int,
        pathway_size: int,
        matched_features=None,
    ):
        self.pathway = pathway
        self.p_value = p_value
        self.adjusted_p_value = adjusted_p_value
        self.enrichment_ratio = enrichment_ratio
        self.input_features = input_features
        self.pathway_size = pathway_size
        self.matched_features = matched_features or []


def test_generate_pathway_aware_dataset_summary(monkeypatch) -> None:
    results = [
        FakeResult(
            pathway=FakePathway("p1", "Pathway One", "KEGG", "Desc", feature_types=["gene"]),
            p_value=0.001,
            adjusted_p_value=0.005,
            enrichment_ratio=2.5,
            input_features=5,
            pathway_size=10,
            matched_features=["A", "B", "C"],
        ),
        FakeResult(
            pathway=FakePathway("p2", "Pathway Two", "Reactome", "", feature_types=["protein"]),
            p_value=0.01,
            adjusted_p_value=0.02,
            enrichment_ratio=1.8,
            input_features=3,
            pathway_size=12,
            matched_features=["X", "Y"],
        ),
    ]

    monkeypatch.setattr(ps, "enrich_dataset_pathways", lambda dataset_page_id, p_value_threshold: results)

    summary = ps.generate_pathway_aware_dataset_summary("ds1", top_pathways=1, p_value_threshold=0.05)

    assert "Pathway Enrichment Analysis" in summary
    assert "Top 1 Enriched Pathways" in summary
    assert "Pathway One" in summary
    assert "Cross-Database Pathway Convergence" in summary


def test_generate_pathway_aware_dataset_summary_no_results(monkeypatch) -> None:
    monkeypatch.setattr(ps, "enrich_dataset_pathways", lambda dataset_page_id, p_value_threshold: [])

    summary = ps.generate_pathway_aware_dataset_summary("ds1")

    assert "No significantly enriched pathways" in summary


def test_generate_pathway_aware_signature_summary(monkeypatch) -> None:
    results = [
        FakeResult(
            pathway=FakePathway("p1", "Pathway One", "KEGG", feature_types=["gene", "protein"]),
            p_value=0.001,
            adjusted_p_value=0.002,
            enrichment_ratio=3.0,
            input_features=4,
            pathway_size=20,
            matched_features=["A", "B"],
        )
    ]

    monkeypatch.setattr(ps, "enrich_signature_pathways", lambda signature_page_id, p_value_threshold: results)

    summary = ps.generate_pathway_aware_signature_summary("sig1")

    assert "Pathway Enrichment Analysis" in summary
    assert "Multi-Omics Pathway Coverage" in summary
    assert "gene, protein" in summary


def test_integrate_pathway_summary_into_cross_omics() -> None:
    combined = ps.integrate_pathway_summary_into_cross_omics("base", "extra")
    assert combined.startswith("base")
    assert "extra" in combined

    unchanged = ps.integrate_pathway_summary_into_cross_omics("base text", "")
    assert unchanged == "base text"

