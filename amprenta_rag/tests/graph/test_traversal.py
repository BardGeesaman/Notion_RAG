from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID, uuid4

from amprenta_rag.graph import traversal as tr


class _N:
    def __init__(
        self,
        source_entity_type: str,
        source_entity_id: UUID,
        target_entity_type: str,
        target_entity_id: UUID,
        relationship_type: str,
    ):
        self.source_entity_type = source_entity_type
        self.source_entity_id = source_entity_id
        self.target_entity_type = target_entity_type
        self.target_entity_id = target_entity_id
        self.relationship_type = relationship_type
        self.confidence = 0.5
        self.evidence_source = "unit_test"
        self.provenance: Optional[Dict[str, Any]] = None


class _FakeBuilder:
    def __init__(self, adj: Dict[tuple[str, UUID], List[_N]]):
        self._adj = adj

    def get_neighbors(self, *, entity_type: str, entity_id: UUID, direction: str, relationship_types=None, limit=2000):
        neighs = self._adj.get((entity_type, entity_id), [])
        if relationship_types:
            relset = set(relationship_types)
            neighs = [n for n in neighs if n.relationship_type in relset]
        return neighs


def test_k_hop_returns_correct_depth(monkeypatch):
    a = uuid4()
    b = uuid4()
    c = uuid4()
    d = uuid4()

    # a - b - c - d (line)
    n_ab = _N("compound", a, "feature", b, "r1")
    n_bc = _N("feature", b, "feature", c, "r1")
    n_cd = _N("feature", c, "pathway", d, "r1")

    adj = {
        ("compound", a): [n_ab],
        ("feature", b): [n_ab, n_bc],
        ("feature", c): [n_bc, n_cd],
        ("pathway", d): [n_cd],
    }

    monkeypatch.setattr(tr, "EdgeBuilder", lambda: _FakeBuilder(adj))
    out = tr.k_hop_subgraph("compound", a, depth=2)
    node_set = set(out["nodes"])
    assert ("compound", a) in node_set
    assert ("feature", b) in node_set
    assert ("feature", c) in node_set
    assert ("pathway", d) not in node_set  # depth=2 stops before d


def test_shortest_path_finds_path(monkeypatch):
    a = uuid4()
    b = uuid4()
    c = uuid4()

    n_ab = _N("compound", a, "feature", b, "r1")
    n_bc = _N("feature", b, "pathway", c, "r1")

    adj = {
        ("compound", a): [n_ab],
        ("feature", b): [n_ab, n_bc],
        ("pathway", c): [n_bc],
    }

    monkeypatch.setattr(tr, "EdgeBuilder", lambda: _FakeBuilder(adj))
    out = tr.shortest_path("compound", a, "pathway", c, timeout_s=1.0)
    assert out is not None
    assert out["nodes"][0] == ("compound", a)
    assert out["nodes"][-1] == ("pathway", c)


def test_no_path_returns_none(monkeypatch):
    a = uuid4()
    b = uuid4()
    c = uuid4()

    n_ab = _N("compound", a, "feature", b, "r1")
    adj = {
        ("compound", a): [n_ab],
        ("feature", b): [n_ab],
        # c disconnected
    }

    monkeypatch.setattr(tr, "EdgeBuilder", lambda: _FakeBuilder(adj))
    out = tr.shortest_path("compound", a, "pathway", c, timeout_s=0.2)
    assert out is None


