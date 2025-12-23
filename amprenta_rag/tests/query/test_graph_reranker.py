from __future__ import annotations

from uuid import uuid4

from amprenta_rag.query.graph_reranker import _clear_path_cache, boost_by_graph
from amprenta_rag.query.rag.models import MatchSummary


def _m(score: float, source_type: str, source_id: str) -> MatchSummary:
    return MatchSummary(
        id="c1",
        score=score,
        source="RAG",
        title="t",
        snippet="s",
        tags=[],
        metadata={"source_type": source_type, "source_id": source_id},
    )


def test_boost_increases_score_for_connected_entities(monkeypatch):
    _clear_path_cache()
    qid = uuid4()
    cid = uuid4()

    # Force query entity extraction
    import amprenta_rag.query.graph_reranker as gr

    monkeypatch.setattr(gr, "_extract_query_entities", lambda q: [("compound", qid)])

    # Fake shortest_path with 1 hop + high confidence
    def _fake_shortest_path(*args, **kwargs):
        return {
            "nodes": [("compound", qid), ("dataset", cid)],
            "edges": [{"confidence": 1.0}],
        }

    monkeypatch.setattr(gr, "UUID", gr.UUID)  # no-op, keeps type checker happy
    monkeypatch.setattr("amprenta_rag.graph.traversal.shortest_path", _fake_shortest_path, raising=False)

    m = _m(1.0, "dataset", str(cid))
    out = boost_by_graph("q", [m], alpha=0.3, max_hops=2)
    assert out[0].score > 1.0


def test_no_boost_for_disconnected_entities(monkeypatch):
    _clear_path_cache()
    qid = uuid4()
    cid = uuid4()

    import amprenta_rag.query.graph_reranker as gr

    monkeypatch.setattr(gr, "_extract_query_entities", lambda q: [("compound", qid)])

    def _fake_shortest_path(*args, **kwargs):
        return None

    monkeypatch.setattr("amprenta_rag.graph.traversal.shortest_path", _fake_shortest_path, raising=False)

    m = _m(0.9, "dataset", str(cid))
    out = boost_by_graph("q", [m], alpha=0.3, max_hops=2)
    assert out[0].score == 0.9


def test_decay_by_path_length(monkeypatch):
    _clear_path_cache()
    qid = uuid4()
    cid1 = uuid4()
    cid2 = uuid4()

    import amprenta_rag.query.graph_reranker as gr

    monkeypatch.setattr(gr, "_extract_query_entities", lambda q: [("compound", qid)])

    # Return length-1 path for cid1, length-2 path for cid2
    def _fake_shortest_path(src_type, src_id, tgt_type, tgt_id, **kwargs):
        if str(tgt_id) == str(cid1):
            return {"nodes": [("compound", qid), ("dataset", cid1)], "edges": [{"confidence": 1.0}]}
        return {
            "nodes": [("compound", qid), ("feature", uuid4()), ("dataset", cid2)],
            "edges": [{"confidence": 1.0}, {"confidence": 1.0}],
        }

    monkeypatch.setattr("amprenta_rag.graph.traversal.shortest_path", _fake_shortest_path, raising=False)

    m1 = _m(1.0, "dataset", str(cid1))
    m2 = _m(1.0, "dataset", str(cid2))
    out = boost_by_graph("q", [m1, m2], alpha=1.0, max_hops=2)

    # 1-hop boost (decay=1.0) should outrank 2-hop boost (decay=0.5)
    assert out[0].metadata.get("source_id") == str(cid1)
    assert out[0].score > out[1].score


