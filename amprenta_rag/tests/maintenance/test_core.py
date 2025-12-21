from __future__ import annotations
import pytest
from types import SimpleNamespace
from amprenta_rag.maintenance import core

def test_delete_all_pinecone_vectors(monkeypatch, capsys):
    mock_config = SimpleNamespace(pinecone=SimpleNamespace(namespace="test-ns"))
    
    class FakeIndex:
        def delete(self, delete_all=False, namespace=""):
            assert delete_all is True
            assert namespace == "test-ns"
            
    monkeypatch.setattr(core, "get_config", lambda: mock_config)
    monkeypatch.setattr(core, "get_pinecone_index", lambda: FakeIndex())
    
    core.delete_all_pinecone_vectors()
    captured = capsys.readouterr()
    assert "Pinecone vectors deleted" in captured.out

def test_archive_all_pages_in_db(capsys):
    core.archive_all_pages_in_db("db123", "Label")
    captured = capsys.readouterr()
    assert "is a no-op" in captured.out

def test_clear_rag_engine_db_success(monkeypatch, capsys):
    mock_config = SimpleNamespace(notion=SimpleNamespace(rag_db_id="rag123"))
    monkeypatch.setattr(core, "get_config", lambda: mock_config)
    monkeypatch.setattr(core, "archive_all_pages_in_db", lambda db_id, label: print(f"Archived {db_id} {label}"))
    
    core.clear_rag_engine_db()
    captured = capsys.readouterr()
    assert "Archived rag123 RAG Engine" in captured.out

def test_clear_rag_engine_db_missing_config(monkeypatch):
    mock_config = SimpleNamespace(notion=SimpleNamespace(rag_db_id=None))
    monkeypatch.setattr(core, "get_config", lambda: mock_config)
    
    with pytest.raises(RuntimeError, match="NOTION_RAG_DB_ID is not configured"):
        core.clear_rag_engine_db()

def test_reset_all(monkeypatch, capsys):
    calls = []
    
    mock_config = SimpleNamespace(notion=SimpleNamespace(
        rag_db_id="rag_id",
        lit_db_id="lit_id"
    ), pinecone=SimpleNamespace(namespace="ns"))
    
    monkeypatch.setattr(core, "get_config", lambda: mock_config)
    monkeypatch.setattr(core, "delete_all_pinecone_vectors", lambda: calls.append("pinecone"))
    monkeypatch.setattr(core, "archive_all_pages_in_db", lambda db_id, label: calls.append(f"archive {label}"))
    
    core.reset_all()
    
    assert "pinecone" in calls
    assert "archive RAG Engine" in calls
    assert "archive Literature" in calls
    
    captured = capsys.readouterr()
    assert "FULL RAG RESET STARTING" in captured.out
    assert "reset complete" in captured.out

def test_reset_all_missing_ids(monkeypatch, capsys):
    calls = []
    
    mock_config = SimpleNamespace(notion=SimpleNamespace(
        rag_db_id=None,
        lit_db_id=None
    ), pinecone=SimpleNamespace(namespace="ns"))
    
    monkeypatch.setattr(core, "get_config", lambda: mock_config)
    monkeypatch.setattr(core, "delete_all_pinecone_vectors", lambda: calls.append("pinecone"))
    
    core.reset_all()
    
    captured = capsys.readouterr()
    assert "NOTION_RAG_DB_ID not set" in captured.out
    assert "NOTION_LIT_DB_ID not set" in captured.out

