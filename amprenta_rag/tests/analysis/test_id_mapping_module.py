from __future__ import annotations

import types

from amprenta_rag.analysis import id_mapping as im


def setup_function(function):
    im.clear_id_mapping_cache()


def test_map_protein_to_uniprot_regex_and_cache(monkeypatch):
    monkeypatch.setattr(im, "requests", None, raising=False)
    pid = "P12345"
    first = im.map_protein_to_uniprot(pid)
    assert first == pid
    # Cached path
    second = im.map_protein_to_uniprot(pid)
    assert second == pid


def test_map_gene_to_kegg_parses_and_caches(monkeypatch):
    class FakeResp:
        status_code = 200
        text = "hsa:1234\tTP53"

    class FakeRequests:
        @staticmethod
        def get(url, timeout=10):
            return FakeResp()

        RequestException = Exception

    monkeypatch.setattr(im, "requests", FakeRequests())
    res = im.map_gene_to_kegg("TP53")
    assert res == "hsa:1234"
    # cached
    monkeypatch.setattr(im, "requests", None, raising=False)
    assert im.map_gene_to_kegg("TP53") == "hsa:1234"


def test_map_metabolite_to_kegg_exact_and_fuzzy(monkeypatch):
    class FakeResp:
        status_code = 200
        text = "cpd:C00001\tglucose; sugar\ncpd:C00002\tother"

    class FakeRequests:
        @staticmethod
        def get(url, timeout=10):
            return FakeResp()

        RequestException = Exception

    monkeypatch.setattr(im, "requests", FakeRequests())
    exact = im.map_metabolite_to_kegg("glucose")
    assert exact == "cpd:C00001"
    im.clear_id_mapping_cache()
    # fuzzy fallback uses first line
    fuzzy = im.map_metabolite_to_kegg("nomatch")
    assert fuzzy == "cpd:C00001"


def test_batch_map_features_kegg_and_reactome(monkeypatch):
    monkeypatch.setattr(im, "time", types.SimpleNamespace(sleep=lambda s: None))
    monkeypatch.setattr(im, "map_gene_to_kegg", lambda f, organism="hsa": "kegg-" + f)
    monkeypatch.setattr(im, "map_protein_to_kegg", lambda f, organism="hsa": None)
    monkeypatch.setattr(im, "map_metabolite_to_kegg", lambda f: None)
    monkeypatch.setattr(im, "map_gene_to_reactome", lambda f: f"r-{f}")
    monkeypatch.setattr(im, "map_protein_to_reactome", lambda f: f"rp-{f}")

    res_kegg = im.batch_map_features_to_pathway_ids({"A"}, "gene", pathway_source="KEGG")
    assert res_kegg["A"] == "kegg-A"

    res_reactome = im.batch_map_features_to_pathway_ids({"B"}, "gene", pathway_source="Reactome")
    assert res_reactome["B"] == "r-B"

    res_unknown = im.batch_map_features_to_pathway_ids({"C"}, "lipid", pathway_source="KEGG")
    assert res_unknown["C"] is None


def test_cache_stats_and_clear(monkeypatch):
    im._id_mapping_cache[("k", "a", "")] = "x"
    im._id_mapping_cache[("k", "b", "")] = None
    stats = im.get_cache_stats()
    assert stats["total_cached"] == 2
    assert stats["successful_mappings"] == 1
    im.clear_id_mapping_cache()
    assert im.get_cache_stats()["total_cached"] == 0

