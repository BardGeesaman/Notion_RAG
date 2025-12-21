from __future__ import annotations

from typing import Dict

from amprenta_rag.analysis.pathway import mapping
from amprenta_rag.analysis.pathway.models import Pathway


def test_map_features_to_kegg_pathways_unsupported_type(monkeypatch):
    monkeypatch.setattr(mapping, "requests", None, raising=False)
    res = mapping.map_features_to_kegg_pathways({"A"}, feature_type="unknown")
    assert res == {}


def test_map_features_to_kegg_pathways_happy_path(monkeypatch):
    calls: Dict[str, int] = {"info": 0}

    def fake_map_gene(feature, organism="hsa"):
        return "hsa:1234"

    class FakeResp:
        status_code = 200
        text = "path:hsa00010\thsa:1234\npath:hsa00020\thsa:1234\n"

    def fake_get(url, timeout=10):
        return FakeResp()

    def fake_fetch(pathway_id):
        calls["info"] += 1
        return {"name": "p", "description": "d"}

    monkeypatch.setattr(mapping, "map_gene_to_kegg", fake_map_gene)
    monkeypatch.setattr(mapping, "_fetch_kegg_pathway_info", fake_fetch)
    monkeypatch.setattr(mapping, "requests", type("R", (), {"get": fake_get})())
    monkeypatch.setattr(mapping, "time", type("T", (), {"sleep": lambda s: None})())

    res = mapping.map_features_to_kegg_pathways({"A"}, feature_type="gene")
    assert "path:hsa00010" in res and isinstance(res["path:hsa00010"], Pathway)
    assert res["path:hsa00010"].features == {"A"}
    assert calls["info"] == 2  # one per unique pathway


def test_map_features_to_reactome_pathways_handles_errors(monkeypatch):
    def fake_map_gene(feature, organism="human"):
        raise ValueError("bad")

    monkeypatch.setattr(mapping, "map_gene_to_reactome", fake_map_gene)
    monkeypatch.setattr(mapping, "requests", type("R", (), {"get": lambda *a, **k: type("Resp", (), {"status_code": 500, "json": lambda self: {}})()})())
    res = mapping.map_features_to_reactome_pathways({"A"}, feature_type="gene")
    assert res == {}


def test_map_features_to_reactome_pathways_happy(monkeypatch):
    def fake_map_gene(feature, organism="human"):
        return "R-HSA-1"

    class FakeResp:
        def __init__(self, data):
            self._data = data
            self.status_code = 200

        def json(self):
            return self._data

    def fake_get(url, timeout=10):
        if "pathways" in url:
            return FakeResp([{"stId": "R-TEST", "displayName": "Pathway"}])
        return FakeResp({})

    monkeypatch.setattr(mapping, "map_gene_to_reactome", fake_map_gene)
    monkeypatch.setattr(mapping, "requests", type("R", (), {"get": fake_get})())

    res = mapping.map_features_to_reactome_pathways({"A"}, feature_type="gene")
    assert "R-TEST" in res
    assert res["R-TEST"].features == {"A"}


def test_fetch_kegg_pathway_info_parses(monkeypatch):
    class FakeResp:
        status_code = 200
        text = "NAME  Glycolysis / Gluconeogenesis\nDESCRIPTION  Energy metabolism\n"

    monkeypatch.setattr(mapping, "requests", type("R", (), {"get": lambda *a, **k: FakeResp()})())
    info = mapping._fetch_kegg_pathway_info("hsa00010")
    assert info["name"] == "Glycolysis / Gluconeogenesis"
    assert "Energy metabolism" in info["description"]

