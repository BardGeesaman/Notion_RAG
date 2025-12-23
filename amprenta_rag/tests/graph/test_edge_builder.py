from __future__ import annotations

from typing import Any, Dict, List, Optional, Type
from uuid import uuid4

from amprenta_rag.graph.edge_builder import EdgeBuilder


class _FakeQuery:
    def __init__(self, items: List[Any], filters: Optional[Dict[str, Any]] = None):
        self._items = items
        self._filters = filters or {}

    def filter_by(self, **kwargs):
        merged = dict(self._filters)
        merged.update(kwargs)
        return _FakeQuery(self._items, merged)

    def all(self):
        out = []
        for it in self._items:
            ok = True
            for k, v in self._filters.items():
                if getattr(it, k) != v:
                    ok = False
                    break
            if ok:
                out.append(it)
        return out

    def first(self):
        rows = self.all()
        return rows[0] if rows else None


class _FakeDB:
    def __init__(self):
        self.edges: List[Any] = []

    def query(self, _model: Type):
        return _FakeQuery(self.edges)

    def add(self, edge):
        # EdgeBuilder instantiates a GraphEdge, but we only need matching attrs.
        self.edges.append(edge)

    def commit(self):
        return None

    def refresh(self, _edge):
        return None


class _DBCtx:
    def __init__(self, db: _FakeDB):
        self._db = db

    def __enter__(self) -> _FakeDB:
        return self._db

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False


def test_create_edge_idempotent():
    # Patch session factory used by EdgeBuilder
    shared = _FakeDB()
    b = EdgeBuilder(session_factory=lambda: _DBCtx(shared))
    src = uuid4()
    tgt = uuid4()

    b.create_edge(
        source_entity_type="compound",
        source_entity_id=src,
        target_entity_type="feature",
        target_entity_id=tgt,
        relationship_type="activity_against",
        confidence=0.4,
        evidence_source="unit_test",
        provenance={"a": 1},
    )
    e2 = b.create_edge(
        source_entity_type="compound",
        source_entity_id=src,
        target_entity_type="feature",
        target_entity_id=tgt,
        relationship_type="activity_against",
        confidence=0.9,
        evidence_source="unit_test",
        provenance={"b": 2},
    )

    assert len(shared.edges) == 1
    assert e2.confidence == 0.9
    assert e2.provenance["a"] == 1
    assert e2.provenance["b"] == 2


def test_get_neighbors():
    shared = _FakeDB()
    b = EdgeBuilder(session_factory=lambda: _DBCtx(shared))
    a = uuid4()
    b_id = uuid4()
    c = uuid4()

    b.create_edge(
        source_entity_type="compound",
        source_entity_id=a,
        target_entity_type="feature",
        target_entity_id=b_id,
        relationship_type="activity_against",
        confidence=0.8,
        evidence_source="unit_test",
    )
    b.create_edge(
        source_entity_type="feature",
        source_entity_id=b_id,
        target_entity_type="pathway",
        target_entity_id=c,
        relationship_type="in_pathway",
        confidence=0.7,
        evidence_source="unit_test",
    )

    neigh = b.get_neighbors(entity_type="feature", entity_id=b_id, direction="both")
    assert len(neigh) == 2
    rels = sorted([n.relationship_type for n in neigh])
    assert rels == ["activity_against", "in_pathway"]


