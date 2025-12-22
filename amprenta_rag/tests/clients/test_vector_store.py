from __future__ import annotations

from contextlib import contextmanager

import pytest

import amprenta_rag.clients.vector_store as vs


class DummyCfg:
    def __init__(self, backend: str):
        self.vector_backend = backend

        class _P:
            namespace = "default"

        self.pinecone = _P()


def test_get_vector_store_returns_pinecone(monkeypatch: pytest.MonkeyPatch):
    class FakeIndex:
        def query(self, **kwargs):
            return {"matches": []}

        def upsert(self, vectors, namespace):
            return None

        def delete(self, **kwargs):
            return None

    monkeypatch.setattr(vs, "get_config", lambda: DummyCfg("pinecone"))
    # Ensure PineconeStore init doesn't hit real pinecone client
    monkeypatch.setattr(vs.PineconeStore, "__init__", lambda self, index=None: setattr(self, "_index", FakeIndex()))

    store = vs.get_vector_store()
    assert isinstance(store, vs.PineconeStore)


def test_get_vector_store_returns_pgvector(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(vs, "get_config", lambda: DummyCfg("pgvector"))
    store = vs.get_vector_store()
    assert isinstance(store, vs.PgVectorStore)


def test_pinecone_store_query_parses_matches(monkeypatch: pytest.MonkeyPatch):
    class FakeIndex:
        def query(self, **kwargs):
            return {
                "matches": [
                    {"id": "c1", "score": 0.9, "metadata": {"a": 1}},
                    {"id": "c2", "score": 0.1, "metadata": {}},
                ]
            }

    monkeypatch.setattr(vs, "get_config", lambda: DummyCfg("pinecone"))
    store = vs.PineconeStore(index=FakeIndex())
    out = store.query([0.0, 1.0], top_k=2)
    assert out[0]["id"] == "c1"
    assert out[0]["score"] == 0.9
    assert out[0]["metadata"] == {"a": 1}


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


