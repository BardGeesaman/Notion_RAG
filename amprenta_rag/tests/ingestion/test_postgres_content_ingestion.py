from __future__ import annotations

import pytest
from types import SimpleNamespace
from unittest.mock import MagicMock
from amprenta_rag.ingestion import postgres_content_ingestion as pci

# --- Mocks ---

class FakeIndex:
    def __init__(self):
        self.deleted = []
        self.upserted = []

    def delete(self, filter, namespace):
        self.deleted.append((filter, namespace))

    def upsert(self, vectors, namespace):
        self.upserted.extend(vectors)

# --- Tests ---

def test_calculate_content_hash():
    h1 = pci._calculate_content_hash("hello")
    h2 = pci._calculate_content_hash("hello")
    h3 = pci._calculate_content_hash("world")
    assert h1 == h2
    assert h1 != h3
    assert len(h1) == 64

def test_content_already_ingested_force(monkeypatch):
    res, _ = pci._content_already_ingested("id", force=True)
    assert res is False

def test_content_already_ingested_not_found(monkeypatch):
    monkeypatch.setattr(pci, "_query_pinecone_by_filter", lambda f: [])
    res, _ = pci._content_already_ingested("id")
    assert res is False

def test_content_already_ingested_found_no_hash_check(monkeypatch):
    match = MagicMock()
    match.metadata = {"content_hash": "abc"}
    monkeypatch.setattr(pci, "_query_pinecone_by_filter", lambda f: [match])
    
    res, h = pci._content_already_ingested("id")
    assert res is True
    assert h == "abc"

def test_content_already_ingested_found_hash_match(monkeypatch):
    match = MagicMock()
    match.metadata = {"content_hash": "abc"}
    monkeypatch.setattr(pci, "_query_pinecone_by_filter", lambda f: [match])
    
    res, h = pci._content_already_ingested("id", current_hash="abc")
    assert res is True
    assert h == "abc"

def test_content_already_ingested_found_hash_mismatch(monkeypatch):
    match = MagicMock()
    match.metadata = {"content_hash": "old"}
    monkeypatch.setattr(pci, "_query_pinecone_by_filter", lambda f: [match])
    
    res, h = pci._content_already_ingested("id", current_hash="new")
    assert res is False
    assert h == "old"

def test_delete_content_vectors(monkeypatch):
    idx = FakeIndex()
    monkeypatch.setattr(pci, "get_pinecone_index", lambda: idx)
    monkeypatch.setattr(pci, "get_config", lambda: SimpleNamespace(pinecone=SimpleNamespace(namespace="ns")))
    
    pci._delete_content_vectors("id1")
    assert len(idx.deleted) == 1
    assert idx.deleted[0][0]["content_id"]["$eq"] == "id1"

def test_ingest_content_direct_to_pinecone_success(monkeypatch):
    idx = FakeIndex()
    monkeypatch.setattr(pci, "get_pinecone_index", lambda: idx)
    monkeypatch.setattr(pci, "get_config", lambda: SimpleNamespace(pinecone=SimpleNamespace(namespace="ns")))
    
    # Mock helpers
    monkeypatch.setattr(pci, "_content_already_ingested", lambda *a, **k: (False, None))
    monkeypatch.setattr(pci, "_delete_content_vectors", lambda cid: None)
    monkeypatch.setattr(pci, "chunk_text", lambda t: ["chunk1", "chunk2"])
    monkeypatch.setattr(pci, "embed_texts", lambda texts: [[0.1]*1536, [0.2]*1536])
    monkeypatch.setattr(pci, "extract_features_from_text", lambda t: {"diseases": ["d1"]})
    monkeypatch.setattr(pci, "sanitize_metadata", lambda m: m)
    
    content = "some content " * 10  # ensure >= 50 chars
    res = pci.ingest_content_direct_to_pinecone(
        content=content,
        source_type="email",
        title="Title",
        metadata={"extra": "val", "source_url": "http://url"},
        content_id="cid",
    )
    
    assert res == ["cid_chunk_000", "cid_chunk_001"]
    assert len(idx.upserted) == 2
    # Check metadata structure of first chunk
    meta = idx.upserted[0]["metadata"]
    assert meta["content_id"] == "cid"
    assert meta["title"] == "Title"
    assert meta["source_type"] == "email"
    assert meta["chunk_index"] == 0
    assert meta["source_url"] == "http://url"

def test_ingest_content_direct_to_pinecone_skipped(monkeypatch):
    monkeypatch.setattr(pci, "_content_already_ingested", lambda *a, **k: (True, "hash"))
    content = "content " * 10
    res = pci.ingest_content_direct_to_pinecone(content, source_type="email", title="T", content_id="cid")
    assert res == []  # idempotent skip returns empty list

def test_ingest_content_direct_to_pinecone_error(monkeypatch):
    monkeypatch.setattr(pci, "_content_already_ingested", lambda *a, **k: (False, None))
    # Chunking fails
    monkeypatch.setattr(pci, "chunk_text", lambda t: (_ for _ in range(1)).throw(ValueError("boom")))
    content = "content " * 10
    with pytest.raises(ValueError, match="boom"):
        pci.ingest_content_direct_to_pinecone(content, source_type="email", title="T", content_id="cid")

