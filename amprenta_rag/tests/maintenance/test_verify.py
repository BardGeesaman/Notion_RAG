from __future__ import annotations

import pytest
from types import SimpleNamespace
from amprenta_rag.maintenance import verify as mv

def test_query_db_no_op():
    assert mv._query_db("id", {}) == []

def test_count_rag_chunks_for_lit_page(monkeypatch):
    # Case 1: No rag_db_id
    monkeypatch.setattr(mv, "get_config", lambda: SimpleNamespace(notion=SimpleNamespace(rag_db_id=None)))
    assert mv._count_rag_chunks_for_lit_page("p1") == 0

    # Case 2: Has rag_db_id (but _query_db is no-op)
    monkeypatch.setattr(mv, "get_config", lambda: SimpleNamespace(notion=SimpleNamespace(rag_db_id="db1")))
    # _query_db returns [] by default
    assert mv._count_rag_chunks_for_lit_page("p1") == 0

def test_count_rag_chunks_for_email_page(monkeypatch):
    # Case 1: No rag_db_id
    monkeypatch.setattr(mv, "get_config", lambda: SimpleNamespace(notion=SimpleNamespace(rag_db_id=None)))
    assert mv._count_rag_chunks_for_email_page("p1") == 0

    # Case 2: Has rag_db_id
    monkeypatch.setattr(mv, "get_config", lambda: SimpleNamespace(notion=SimpleNamespace(rag_db_id="db1")))
    assert mv._count_rag_chunks_for_email_page("p1") == 0

def test_verify_rag_metadata(monkeypatch, capsys):
    # Mock config
    cfg = SimpleNamespace(notion=SimpleNamespace(lit_db_id="lit_db", email_db_id="email_db", rag_db_id="rag_db"))
    monkeypatch.setattr(mv, "get_config", lambda: cfg)

    # Mock _query_db to return some fake pages
    def fake_query(db_id, payload):
        if db_id == "lit_db":
            return [{
                "id": "p1",
                "properties": {
                    "Zotero Item Key": {"rich_text": [{"plain_text": "KEY1"}]},
                    "Title": {"title": [{"plain_text": "Paper 1"}]}
                }
            }]
        if db_id == "email_db":
            return [{
                "id": "e1",
                "properties": {
                    "Title": {"title": [{"plain_text": "Email 1"}]},
                    "Type": {"select": {"name": "Email"}}
                }
            }]
        if db_id == "rag_db":
            # Return fake chunks to simulate count > 0
            return [{"id": "chunk1"}]
        return []

    monkeypatch.setattr(mv, "_query_db", fake_query)

    mv.verify_rag_metadata()
    
    captured = capsys.readouterr()
    assert "Verifying RAG metadata" in captured.out
    assert "[KEY1] Paper 1 -> 1 RAG chunks" in captured.out
    assert "[Email] Email 1 -> 1 RAG chunks" in captured.out

