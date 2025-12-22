from __future__ import annotations

from types import SimpleNamespace

import pytest

from amprenta_rag.analysis import id_mapping as idm


def _reset_cache():
    idm._id_mapping_cache.clear()


def test_map_protein_to_uniprot_cached(monkeypatch):
    _reset_cache()
    idm._id_mapping_cache[("uniprot", "X", "")] = "X12345"
    monkeypatch.setattr(idm, "requests", None)
    assert idm.map_protein_to_uniprot("X") == "X12345"


def test_map_protein_to_uniprot_direct_pattern():
    _reset_cache()
    assert idm.map_protein_to_uniprot("P12345") == "P12345"


def test_map_protein_to_uniprot_api_success(monkeypatch):
    _reset_cache()
    calls = {"post": 0, "get": 0}

    class FakeResp:
        def __init__(self, status_code, payload=None, text=""):
            self.status_code = status_code
            self._payload = payload or {}
            self.text = text

        def json(self):
            return self._payload

    def fake_post(url, data=None, timeout=10):
        calls["post"] += 1
        return FakeResp(200, {"jobId": "job1"})

    def fake_get(url, timeout=10):
        calls["get"] += 1
        return FakeResp(200, {"results": [{"to": {"primaryAccession": "U12345"}}]})

    monkeypatch.setattr(idm, "requests", SimpleNamespace(post=fake_post, get=fake_get))
    monkeypatch.setattr(idm.time, "sleep", lambda *_: None)

    assert idm.map_protein_to_uniprot("TP53") == "U12345"
    assert calls["post"] >= 1
    assert calls["get"] >= 1


def test_map_gene_to_kegg_success(monkeypatch):
    _reset_cache()

    class FakeResp:
        def __init__(self, text):
            self.status_code = 200
            self.text = text

    monkeypatch.setattr(idm, "requests", SimpleNamespace(get=lambda url, timeout=10: FakeResp("hsa:123\tTP53")))
    assert idm.map_gene_to_kegg("TP53") == "hsa:123"


def test_map_protein_to_kegg_via_uniprot(monkeypatch):
    _reset_cache()
    monkeypatch.setattr(idm, "map_protein_to_uniprot", lambda pid: "U123")

    class FakeResp:
        def __init__(self, text):
            self.status_code = 200
            self.text = text

    monkeypatch.setattr(idm, "requests", SimpleNamespace(get=lambda url, timeout=10: FakeResp("uniprot:U123\thsa:999")))
    assert idm.map_protein_to_kegg("anything") == "hsa:999"


def test_map_metabolite_to_kegg_exact(monkeypatch):
    _reset_cache()

    class FakeResp:
        def __init__(self, text):
            self.status_code = 200
            self.text = text

    monkeypatch.setattr(
        idm,
        "requests",
        SimpleNamespace(get=lambda url, timeout=10: FakeResp("cpd:C00031\tglucose; d-glucose")),
    )
    assert idm.map_metabolite_to_kegg("glucose") == "cpd:C00031"


def test_map_gene_to_reactome_valid_and_invalid():
    assert idm.map_gene_to_reactome("TP53") == "TP53"
    assert idm.map_gene_to_reactome("bad symbol!") is None


def test_batch_map_features_handles_types(monkeypatch):
    _reset_cache()
    monkeypatch.setattr(idm, "map_gene_to_kegg", lambda f, organism="hsa": f"kegg:{f}")
    monkeypatch.setattr(idm, "map_protein_to_kegg", lambda f, organism="hsa": f"kegg_p:{f}")
    monkeypatch.setattr(idm, "map_metabolite_to_kegg", lambda f: f"kegg_m:{f}")
    monkeypatch.setattr(idm, "map_gene_to_reactome", lambda f: f"re:{f}")
    monkeypatch.setattr(idm, "map_protein_to_reactome", lambda f: f"re_p:{f}")

    genes = idm.batch_map_features_to_pathway_ids({"G1"}, "gene", pathway_source="KEGG")
    prots = idm.batch_map_features_to_pathway_ids({"P1"}, "protein", pathway_source="Reactome")
    mets = idm.batch_map_features_to_pathway_ids({"M1"}, "metabolite", pathway_source="KEGG")

    assert genes["G1"] == "kegg:G1"
    assert prots["P1"] == "re_p:P1"
    assert mets["M1"] == "kegg_m:M1"

