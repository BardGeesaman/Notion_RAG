from __future__ import annotations
import pytest
from types import SimpleNamespace
from amprenta_rag.maintenance import verify

def test_query_db_returns_empty():
    assert verify._query_db("db_id", {}) == []

def test_count_rag_chunks_for_lit_page_no_rag_db(monkeypatch):
    # Mock config to return None for rag_db_id
    mock_config = SimpleNamespace(notion=SimpleNamespace(rag_db_id=None))
    monkeypatch.setattr(verify, "get_config", lambda: mock_config)
    
    assert verify._count_rag_chunks_for_lit_page("page_id") == 0

def test_count_rag_chunks_for_lit_page_with_rag_db(monkeypatch):
    # Mock config
    mock_config = SimpleNamespace(notion=SimpleNamespace(rag_db_id="rag_db"))
    monkeypatch.setattr(verify, "get_config", lambda: mock_config)
    
    # _query_db returns empty list by default in the module
    assert verify._count_rag_chunks_for_lit_page("page_id") == 0

def test_count_rag_chunks_for_email_page_no_rag_db(monkeypatch):
    mock_config = SimpleNamespace(notion=SimpleNamespace(rag_db_id=None))
    monkeypatch.setattr(verify, "get_config", lambda: mock_config)
    assert verify._count_rag_chunks_for_email_page("page_id") == 0

def test_count_rag_chunks_for_email_page_with_rag_db(monkeypatch):
    mock_config = SimpleNamespace(notion=SimpleNamespace(rag_db_id="rag_db"))
    monkeypatch.setattr(verify, "get_config", lambda: mock_config)
    assert verify._count_rag_chunks_for_email_page("page_id") == 0

def test_verify_rag_metadata_runs(monkeypatch, capsys):
    # Mock config with DB IDs
    mock_config = SimpleNamespace(notion=SimpleNamespace(
        lit_db_id="lit_db",
        email_db_id="email_db",
        rag_db_id="rag_db"
    ))
    monkeypatch.setattr(verify, "get_config", lambda: mock_config)
    
    # Since _query_db returns empty list, the loops over pages won't run, 
    # but the function should complete and print start/end messages.
    verify.verify_rag_metadata()
    
    captured = capsys.readouterr()
    assert "Verifying RAG metadata consistency" in captured.out
    assert "Verification complete" in captured.out

def test_verify_rag_metadata_with_mocked_pages(monkeypatch, capsys):
    # Mock config
    mock_config = SimpleNamespace(notion=SimpleNamespace(
        lit_db_id="lit_db",
        email_db_id="email_db",
        rag_db_id="rag_db"
    ))
    monkeypatch.setattr(verify, "get_config", lambda: mock_config)
    
    # Mock _query_db to return some pages for specific calls
    def fake_query_db(db_id, payload):
        if db_id == "lit_db":
            return [{
                "id": "lit_page_1",
                "properties": {
                    "Zotero Item Key": {"rich_text": [{"plain_text": "KEY1"}]},
                    "Title": {"title": [{"plain_text": "Paper Title"}]}
                }
            }]
        elif db_id == "email_db":
            return [{
                "id": "email_page_1",
                "properties": {
                    "Title": {"title": [{"plain_text": "Email Subject"}]},
                    "Type": {"select": {"name": "Email"}}
                }
            }]
        return [] # RAG DB queries return empty chunks
        
    monkeypatch.setattr(verify, "_query_db", fake_query_db)
    
    verify.verify_rag_metadata()
    
    captured = capsys.readouterr()
    assert "Paper Title" in captured.out
    assert "Email Subject" in captured.out
