from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

import pytest


def test_calculate_pic50_converts_correctly():
    from amprenta_rag.graph.compound_target_network import calculate_pic50

    assert calculate_pic50(1.0) == pytest.approx(9.0, abs=1e-9)
    assert calculate_pic50(10.0) == pytest.approx(8.0, abs=1e-9)
    assert calculate_pic50(1000.0) == pytest.approx(6.0, abs=1e-9)


def test_edge_width_from_pic50_maps_range():
    from amprenta_rag.graph.compound_target_network import edge_width_from_pic50

    assert edge_width_from_pic50(5.0) == pytest.approx(1.0, abs=1e-9)
    assert edge_width_from_pic50(9.0) == pytest.approx(10.0, abs=1e-9)
    # Clamp
    assert edge_width_from_pic50(1.0) == pytest.approx(1.0, abs=1e-9)
    assert edge_width_from_pic50(20.0) == pytest.approx(10.0, abs=1e-9)


class _FakeQuery:
    def __init__(self, items):
        self._items = list(items or [])

    def filter(self, *args, **kwargs):  # noqa: ANN001
        return self

    def filter_by(self, **kwargs):  # noqa: ANN001
        return self

    def order_by(self, *args, **kwargs):  # noqa: ANN001
        return self

    def limit(self, *args, **kwargs):  # noqa: ANN001
        return self

    def all(self):
        return list(self._items)

    def first(self):
        rows = self.all()
        return rows[0] if rows else None


class _FakeDB:
    def __init__(self, *, edges=None, expand_from_compound=None, expand_from_target=None):
        self._edges = list(edges or [])
        self._expand_compound = list(expand_from_compound or [])
        self._expand_target = list(expand_from_target or [])

    def query(self, *args, **kwargs):  # noqa: ANN001
        # Import here so the test module can run without SQLAlchemy in strange environments.
        from amprenta_rag.database.models import Compound, Feature, GraphEdge

        if len(args) == 1 and args[0] is GraphEdge:
            return _FakeQuery(self._edges)
        if len(args) == 1 and args[0] is GraphEdge.target_entity_id:
            return _FakeQuery([(x,) for x in self._expand_compound])
        if len(args) == 1 and args[0] is GraphEdge.source_entity_id:
            return _FakeQuery([(x,) for x in self._expand_target])
        # Label queries (Compound/Feature): return empty; network should still render with fallback.
        if len(args) == 2 and args[0] is Compound.id and args[1] is Compound.compound_id:
            return _FakeQuery([])
        if len(args) == 2 and args[0] is Feature.id and args[1] is Feature.name:
            return _FakeQuery([])
        return _FakeQuery([])


class _DBCtx:
    def __init__(self, db: _FakeDB):
        self._db = db

    def __enter__(self) -> _FakeDB:
        return self._db

    def __exit__(self, exc_type, exc, tb) -> bool:  # noqa: ANN001
        return False


def _edge(*, src, tgt, prov):
    return SimpleNamespace(
        id=uuid4(),
        source_entity_type="compound",
        source_entity_id=src,
        target_entity_type="feature",
        target_entity_id=tgt,
        relationship_type="activity_against",
        confidence=0.7,
        evidence_source="activity_results",
        provenance=prov,
        updated_at=None,
    )


def test_get_compound_target_network_returns_nodes_and_edges():
    from amprenta_rag.graph.compound_target_network import CompoundTargetNetworkService

    c1, c2 = uuid4(), uuid4()
    t1, t2 = uuid4(), uuid4()
    edges = [
        _edge(src=c1, tgt=t1, prov={"best_ic50_nm": 10.0, "assays_count": 2, "activity_type": "IC50"}),
        _edge(src=c2, tgt=t2, prov={"best_ic50_nm": 100.0, "assays_count": 1, "activity_type": "IC50"}),
    ]
    shared = _FakeDB(edges=edges)
    svc = CompoundTargetNetworkService(session_factory=lambda: _DBCtx(shared))
    out = svc.get_compound_target_network(compound_ids=None, target_ids=None, filters=None, max_nodes=500)

    assert set(out.keys()) >= {"nodes", "edges", "meta"}
    assert len(out["edges"]) == 2
    assert len(out["nodes"]) == 4
    assert out["meta"]["node_count"] == 4


def test_expand_from_compound_returns_target_ids():
    from amprenta_rag.graph.compound_target_network import CompoundTargetNetworkService

    c1 = uuid4()
    t_ids = [uuid4(), uuid4()]
    shared = _FakeDB(edges=[], expand_from_compound=t_ids)
    svc = CompoundTargetNetworkService(session_factory=lambda: _DBCtx(shared))
    out = svc.expand_from_compound(c1)
    assert out == t_ids


def test_expand_from_target_returns_compound_ids():
    from amprenta_rag.graph.compound_target_network import CompoundTargetNetworkService

    t1 = uuid4()
    c_ids = [uuid4(), uuid4()]
    shared = _FakeDB(edges=[], expand_from_target=c_ids)
    svc = CompoundTargetNetworkService(session_factory=lambda: _DBCtx(shared))
    out = svc.expand_from_target(t1)
    assert out == c_ids


def test_filters_ic50_range_applied():
    from amprenta_rag.graph.compound_target_network import CompoundTargetNetworkService

    c1, c2 = uuid4(), uuid4()
    t1, t2 = uuid4(), uuid4()
    edges = [
        _edge(src=c1, tgt=t1, prov={"best_ic50_nm": 10.0, "assays_count": 1, "activity_type": "IC50"}),
        _edge(src=c2, tgt=t2, prov={"best_ic50_nm": 10000.0, "assays_count": 1, "activity_type": "IC50"}),
    ]
    shared = _FakeDB(edges=edges)
    svc = CompoundTargetNetworkService(session_factory=lambda: _DBCtx(shared))
    out = svc.get_compound_target_network(
        filters={"ic50_range": {"min_nm": 1.0, "max_nm": 100.0}},
        max_nodes=500,
    )
    assert len(out["edges"]) == 1
    assert out["edges"][0]["data"]["best_ic50_nm"] == 10.0


def test_max_nodes_limit_enforced():
    from amprenta_rag.graph.compound_target_network import CompoundTargetNetworkService

    edges = []
    for _ in range(10):
        edges.append(_edge(src=uuid4(), tgt=uuid4(), prov={"best_ic50_nm": 10.0, "assays_count": 1}))
    shared = _FakeDB(edges=edges)
    svc = CompoundTargetNetworkService(session_factory=lambda: _DBCtx(shared))
    out = svc.get_compound_target_network(max_nodes=3)
    assert out["meta"]["node_count"] <= 3
    assert len(out["nodes"]) <= 3


