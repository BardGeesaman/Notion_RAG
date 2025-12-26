from __future__ import annotations

import pytest


KGML_SAMPLE = """<?xml version="1.0"?>
<pathway name="path:hsa00010" org="hsa" number="00010" title="Glycolysis / Gluconeogenesis">
  <entry id="1" name="hsa:7157" type="gene">
    <graphics name="TP53" x="200" y="100" type="rectangle"/>
  </entry>
  <entry id="2" name="cpd:C00031" type="compound">
    <graphics name="Glucose" x="100" y="50" type="circle"/>
  </entry>
  <entry id="3" name="path:hsa00020" type="map">
    <graphics name="Citrate cycle" x="50" y="100" type="roundrectangle"/>
  </entry>
  <entry id="4" name="br:br08901" type="brite">
    <graphics name="BRITE" x="10" y="10" type="rectangle"/>
  </entry>
  <relation entry1="1" entry2="2" type="PPrel">
    <subtype name="activation" value="--&gt;"/>
  </relation>
  <relation entry1="2" entry2="3" type="PCrel">
    <subtype name="inhibition" value="--|"/>
  </relation>
  <reaction id="r1" name="rn:R00001" type="reversible">
    <substrate id="2"/>
    <product id="1"/>
  </reaction>
</pathway>
"""


def test_fetch_kegg_pathway_structure_parses_kgml(monkeypatch):
    from amprenta_rag.analysis.pathway import structure as ps

    monkeypatch.setattr(ps, "_rate_limit_kegg", lambda: None)

    class Resp:
        status_code = 200
        text = KGML_SAMPLE

    monkeypatch.setattr(ps.requests, "get", lambda url, timeout=20: Resp())
    st = ps.fetch_kegg_pathway_structure("hsa00010")
    assert st.pathway_id == "hsa00010"
    assert "Glycolysis" in st.name
    # brite entry should be skipped
    assert {n.type for n in st.nodes} >= {"gene", "compound", "sub-pathway"}
    assert len(st.edges) >= 3


def test_node_type_mapping():
    from amprenta_rag.analysis.pathway.structure import _node_type

    assert _node_type("gene") == "gene"
    assert _node_type("ortholog") == "gene"
    assert _node_type("enzyme") == "gene"
    assert _node_type("compound") == "compound"
    assert _node_type("map") == "sub-pathway"
    assert _node_type("reaction") is None
    assert _node_type("brite") is None


def test_coordinate_normalization(monkeypatch):
    from amprenta_rag.analysis.pathway import structure as ps

    monkeypatch.setattr(ps, "_rate_limit_kegg", lambda: None)

    class Resp:
        status_code = 200
        text = KGML_SAMPLE

    monkeypatch.setattr(ps.requests, "get", lambda url, timeout=20: Resp())
    st = ps.fetch_kegg_pathway_structure("hsa00010")
    # Max x in sample is 200, max y is 100 -> normalized should be in [0,1]
    for n in st.nodes:
        assert 0.0 <= n.x <= 1.0
        assert 0.0 <= n.y <= 1.0
    tp53 = next(n for n in st.nodes if n.name == "TP53")
    assert tp53.x == pytest.approx(1.0)
    assert tp53.y == pytest.approx(1.0)


def test_edge_style_activation():
    from amprenta_rag.analysis.pathway.structure import PathwayEdge, edge_style

    s = edge_style(PathwayEdge(source="1", target="2", type="PPrel", subtype="activation"))
    assert s["color"] == "#2E8B57"
    assert s["style"] == "solid"


def test_edge_style_inhibition():
    from amprenta_rag.analysis.pathway.structure import PathwayEdge, edge_style

    s = edge_style(PathwayEdge(source="1", target="2", type="PPrel", subtype="inhibition"))
    assert s["color"] == "#C73E1D"
    assert s["style"] == "dashed"


def test_cache_and_load_structure(tmp_path, monkeypatch):
    from amprenta_rag.analysis.pathway import structure as ps

    monkeypatch.setattr(ps, "CACHE_DIR", tmp_path)
    st = ps.PathwayStructure(
        pathway_id="hsa00010",
        name="x",
        nodes=[ps.PathwayNode(id="1", name="TP53", type="gene", x=0.1, y=0.2, kegg_ids=["hsa:7157"])],
        edges=[ps.PathwayEdge(source="1", target="2", type="PPrel", subtype="activation")],
        organism="hsa",
    )
    ps.cache_pathway_structure("hsa00010", st)
    loaded = ps.load_cached_structure("hsa00010")
    assert loaded is not None
    assert loaded.pathway_id == "hsa00010"
    assert loaded.nodes[0].name == "TP53"


def test_get_pathway_structure_uses_cache(tmp_path, monkeypatch):
    from amprenta_rag.analysis.pathway import structure as ps

    monkeypatch.setattr(ps, "CACHE_DIR", tmp_path)
    st = ps.PathwayStructure(pathway_id="hsa00010", name="x", nodes=[], edges=[], organism="hsa")
    ps.cache_pathway_structure("hsa00010", st)

    called = {"n": 0}

    def fake_fetch(pid):  # noqa: ANN001
        called["n"] += 1
        return ps._empty_structure(pid)

    monkeypatch.setattr(ps, "fetch_kegg_pathway_structure", fake_fetch)
    out = ps.get_pathway_structure("hsa00010")
    assert out.pathway_id == "hsa00010"
    assert called["n"] == 0


