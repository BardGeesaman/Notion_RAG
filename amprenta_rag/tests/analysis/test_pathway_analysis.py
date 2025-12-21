from __future__ import annotations

from amprenta_rag.analysis import pathway_analysis as pa


class FakeResp:
    def __init__(self, text="", status=200, json_data=None):
        self.text = text
        self.status_code = status
        self._json = json_data or {}

    def json(self):
        return self._json


def test_map_features_to_kegg_pathways_invalid_type(monkeypatch):
    monkeypatch.setattr(pa, "requests", None, raising=False)
    res = pa.map_features_to_kegg_pathways({"A"}, feature_type="unknown")
    assert res == {}


def test_map_features_to_kegg_pathways_parses(monkeypatch):
    def fake_map_gene(feature, organism="hsa"):
        return "hsa:1"

    class FakeRequests:
        RequestException = Exception
        ConnectionError = Exception
        TimeoutError = Exception

        @staticmethod
        def get(*a, **k):
            return FakeResp("path:hsa00010\thsa:1", 200)

    monkeypatch.setattr(pa, "map_gene_to_kegg", fake_map_gene)
    monkeypatch.setattr(pa, "_fetch_kegg_pathway_info", lambda pid: {"name": "p"})
    monkeypatch.setattr(pa, "requests", FakeRequests(), raising=False)
    monkeypatch.setattr(pa, "time", type("T", (), {"sleep": lambda s: None})(), raising=False)

    res = pa.map_features_to_kegg_pathways({"A"}, feature_type="gene")
    assert "path:hsa00010" in res


def test_map_features_to_reactome_pathways_handles_missing(monkeypatch):
    monkeypatch.setattr(pa, "map_gene_to_reactome", lambda f, organism="human": None)
    res = pa.map_features_to_reactome_pathways({"A"}, feature_type="gene")
    assert res == {}


def test_map_features_to_reactome_pathways_parses(monkeypatch):
    monkeypatch.setattr(pa, "map_gene_to_reactome", lambda f, organism="human": "R-1")

    def fake_get(url, timeout=10, **kwargs):
        if "mapping" in url:
            return FakeResp(json_data={"pathways": ["R-A"]})
        return FakeResp(status=404)

    class FakeRequests:
        RequestException = Exception
        ConnectionError = Exception
        TimeoutError = Exception

        @staticmethod
        def get(url, timeout=10, **kwargs):
            return fake_get(url, timeout=timeout, **kwargs)

    monkeypatch.setattr(pa, "requests", FakeRequests(), raising=False)

    res = pa.map_features_to_reactome_pathways({"A"}, feature_type="gene")
    assert "R-A" in res


def test_perform_pathway_enrichment_empty():
    res = pa.perform_pathway_enrichment(set(), [])
    assert res == []

