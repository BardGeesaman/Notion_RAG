from __future__ import annotations

from amprenta_rag.graph.analytics import compute_graph_analytics


def test_degree_centrality_computed():
    nodes = [{"id": "a"}, {"id": "b"}, {"id": "c"}]
    edges = [{"source": "a", "target": "b"}, {"source": "b", "target": "c"}]
    out = compute_graph_analytics(nodes, edges, ["degree_centrality"])
    dc = out.get("degree_centrality") or {}
    assert "a" in dc and "b" in dc and "c" in dc
    assert dc["b"] > dc["a"] and dc["b"] > dc["c"]


def test_communities_detected():
    nodes = [{"id": "a"}, {"id": "b"}, {"id": "c"}, {"id": "d"}]
    # two disconnected components -> should yield at least 2 communities
    edges = [{"source": "a", "target": "b"}, {"source": "c", "target": "d"}]
    out = compute_graph_analytics(nodes, edges, ["communities"])
    comms = out.get("communities") or {}
    assert set(comms.keys()) == {"a", "b", "c", "d"}
    assert comms["a"] == comms["b"]
    assert comms["c"] == comms["d"]
    assert comms["a"] != comms["c"]
    assert isinstance(out.get("modularity"), float)


def test_empty_graph_handled():
    out = compute_graph_analytics([], [], ["degree_centrality", "communities"])
    assert out.get("degree_centrality") == {}
    assert out.get("communities") == {}
    assert out.get("modularity") == 0.0


