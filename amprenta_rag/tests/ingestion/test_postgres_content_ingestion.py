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
    """Force flag causes re-ingestion."""
    res, _ = pci._content_already_ingested("id", force=True)
    assert res is False

def test_content_already_ingested_stub(monkeypatch):
    """Test that content_already_ingested is now a stub that always returns False.
    
    Pinecone deduplication removed - function now always returns (False, None).
    NOTE: pgvector-based deduplication tracked in ROADMAP.
    """
    res, h = pci._content_already_ingested("id")
    assert res is False
    assert h is None

def test_delete_content_vectors(monkeypatch):
    """Test _delete_content_vectors still works with vector store."""
    idx = FakeIndex()
    monkeypatch.setattr(pci, "get_vector_store", lambda: idx)
    
    pci._delete_content_vectors("id1")
    assert len(idx.deleted) == 1
    assert idx.deleted[0][0]["content_id"]["$eq"] == "id1"

