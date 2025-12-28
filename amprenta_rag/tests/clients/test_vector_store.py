from __future__ import annotations

from contextlib import contextmanager

import pytest

import amprenta_rag.clients.vector_store as vs


class DummyCfg:
    def __init__(self, backend: str):
        self.vector_backend = backend


def test_get_vector_store_returns_pgvector(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(vs, "get_config", lambda: DummyCfg("pgvector"))
    store = vs.get_vector_store()
    assert isinstance(store, vs.PgVectorStore)


def test_get_vector_store_pinecone_raises_error(monkeypatch: pytest.MonkeyPatch):
    """Test that attempting to use deprecated pinecone backend raises error."""
    monkeypatch.setattr(vs, "get_config", lambda: DummyCfg("pinecone"))
    with pytest.raises(RuntimeError, match="Pinecone backend is deprecated"):
        vs.get_vector_store()


def test_pgvector_store_query_executes_sql(monkeypatch: pytest.MonkeyPatch):
    executed = {}

    class FakeResult:
        def mappings(self):
            return self

        def all(self):
            return [
                {"chunk_id": "c1", "chunk_metadata": {"k": "v"}, "score": 0.5},
                {"chunk_id": "c2", "chunk_metadata": {}, "score": 0.1},
            ]

    class FakeDB:
        def execute(self, stmt, params):
            executed["sql"] = str(stmt)
            executed["params"] = params
            return FakeResult()

        def commit(self):
            return None

    @contextmanager
    def fake_db_session():
        yield FakeDB()

    store = vs.PgVectorStore(session_factory=fake_db_session)
    out = store.query([0.0, 1.0], top_k=2, filter={"zotero_item_key": {"$eq": "abc"}})
    assert "FROM rag_chunks" in executed["sql"]
    assert executed["params"]["top_k"] == 2
    assert out[0]["id"] == "c1"
    assert out[0]["metadata"] == {"k": "v"}


def test_pgvector_store_query_includes_namespace_filter():
    executed = {}

    class FakeResult:
        def mappings(self):
            return self

        def all(self):
            return []

    class FakeDB:
        def execute(self, stmt, params):
            executed["sql"] = str(stmt)
            executed["params"] = params
            return FakeResult()

        def commit(self):
            return None

    @contextmanager
    def fake_db_session():
        yield FakeDB()

    store = vs.PgVectorStore(session_factory=fake_db_session)
    store.query([0.0, 1.0], top_k=1, namespace="ns1")
    assert "namespace" in executed["params"]
    assert executed["params"]["namespace"] == "ns1"


