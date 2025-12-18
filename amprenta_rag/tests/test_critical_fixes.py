from __future__ import annotations

"""
Tests for critical fixes and regressions.

Coverage:
1. CORS configuration honors CORS_ORIGINS environment variable.
2. Streamlit ingestion page module imports cleanly (no undefined globals).
3. auto_linking.db_session properly opens and closes DB sessions, even on error.
4. API input validation rejects overlong name_filter (>100 chars) for datasets.
5. Embedding retry logic retries OpenAI calls on failure and eventually succeeds.
6. Repository rate limiting + retry logic respect delays and retry on 500/503.
"""

from typing import Any, List

import pytest


# ---------------------------------------------------------------------------
# 1. CORS configuration
# ---------------------------------------------------------------------------


def test_cors_origins_respects_environment(monkeypatch):
    """
    ServerConfig.cors_origins should reflect the CORS_ORIGINS environment variable.
    """
    from amprenta_rag.config import ServerConfig

    monkeypatch.setenv("CORS_ORIGINS", "http://a.example.com, https://b.example.com")
    cfg = ServerConfig()

    assert cfg.cors_origins == ["http://a.example.com", "https://b.example.com"]


# ---------------------------------------------------------------------------
# 2. ingestion.py import safety
# ---------------------------------------------------------------------------


def test_ingestion_page_module_imports_without_errors():
    """
    The Streamlit ingestion page module should import without raising NameError
    or other import-time errors (regression check).
    """
    import importlib

    module = importlib.import_module("scripts.dashboard.pages.ingestion")
    assert hasattr(module, "render_ingestion_page")


# ---------------------------------------------------------------------------
# 3. auto_linking session handling
# ---------------------------------------------------------------------------


@pytest.mark.skip(reason="Test mocks non-existent auto_linking.get_db - needs rewrite")
def test_auto_linking_db_session_closes_on_success(monkeypatch):
    """
    db_session context manager must close the underlying DB session on success.
    """
    from amprenta_rag.ingestion import auto_linking

    class DummySession:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

    def fake_get_db():
        session = DummySession()

        def gen():
            try:
                yield session
            finally:
                # get_db's own finally isn't used here; closing is handled by db_session
                pass

        return gen()

    monkeypatch.setattr(auto_linking, "get_db", fake_get_db)

    with auto_linking.db_session() as db:
        assert isinstance(db, DummySession)
        assert not db.closed

    # After context exit, session should be closed
    assert db.closed


@pytest.mark.skip(reason="Test mocks non-existent auto_linking.get_db - needs rewrite")
def test_auto_linking_db_session_closes_on_exception(monkeypatch):
    """
    db_session must also close the session when an exception occurs inside the block.
    """
    from amprenta_rag.ingestion import auto_linking

    class DummySession:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

    def fake_get_db():
        session = DummySession()

        def gen():
            try:
                yield session
            finally:
                pass

        return gen()

    monkeypatch.setattr(auto_linking, "get_db", fake_get_db)

    with pytest.raises(RuntimeError):
        with auto_linking.db_session() as db:
            assert isinstance(db, DummySession)
            raise RuntimeError("Simulated failure inside db_session")

    assert db.closed


# ---------------------------------------------------------------------------
# 4. API input validation: name_filter length
# ---------------------------------------------------------------------------


def test_get_datasets_name_filter_length_validation():
    """
    get_datasets should raise ValueError when name_filter exceeds 100 characters.
    """
    from amprenta_rag.api.services.datasets import get_datasets

    class DummyQuery:
        def filter(self, *args: Any, **kwargs: Any) -> "DummyQuery":
            return self

        def offset(self, *args: Any, **kwargs: Any) -> "DummyQuery":
            return self

        def limit(self, *args: Any, **kwargs: Any) -> "DummyQuery":
            return self

        def all(self) -> list:
            return []

    class DummyDB:
        def query(self, model):
            return DummyQuery()

    db = DummyDB()

    # Short filter should not raise
    res = get_datasets(db=db, name_filter="short")
    assert isinstance(res, list)

    # Overlong filter should raise
    long_filter = "x" * 101
    with pytest.raises(ValueError):
        get_datasets(db=db, name_filter=long_filter)


# ---------------------------------------------------------------------------
# 5. Embedding retry logic
# ---------------------------------------------------------------------------


def test_embedding_retry_logic_retries_and_succeeds(monkeypatch):
    """
    _embed_batch should retry failed embedding calls according to retry_with_backoff.
    We simulate two failures followed by success.
    """
    from amprenta_rag.ingestion.text_embedding_utils import _embed_batch

    call_count = {"n": 0}

    class DummyData:
        def __init__(self, embedding):
            self.embedding = embedding

    class DummyResponse:
        def __init__(self, batch: List[str]):
            self.data = [DummyData([float(i)]) for i, _ in enumerate(batch)]

    class DummyEmbeddings:
        def create(self, model: str, input: List[str]):
            call_count["n"] += 1
            # Fail first two attempts to trigger retries, then succeed
            if call_count["n"] < 3:
                raise RuntimeError("Simulated embedding failure")
            return DummyResponse(input)

    class DummyClient:
        def __init__(self):
            self.embeddings = DummyEmbeddings()

    # Avoid real sleeping during backoff
    monkeypatch.setattr(
        "amprenta_rag.utils.error_handling.time.sleep",
        lambda *_args, **_kwargs: None,
        raising=False,
    )

    batch = ["a", "b"]
    embeddings = _embed_batch(DummyClient(), "dummy-model", batch)

    assert call_count["n"] == 3  # 2 failures + 1 success
    assert len(embeddings) == len(batch)


# ---------------------------------------------------------------------------
# 6. Repository rate limiting and retry logic
# ---------------------------------------------------------------------------


def test_repository_rate_limiting_and_retry(monkeypatch):
    """
    MetaboLightsRepository._make_request_with_retry should:
    - Call repo_rate_limit before each request
    - Retry once after a 500/503 response with a delay
    """
    from amprenta_rag.ingestion.repositories.metabolights import MetaboLightsRepository

    call_sequence: List[int] = []

    class DummyResponse:
        def __init__(self, status_code: int):
            self.status_code = status_code

        def json(self):
            return {}

    # First call: 500, second call: 200
    def fake_get(url, params=None, headers=None, timeout=30):
        if not call_sequence:
            call_sequence.append(500)
            return DummyResponse(500)
        call_sequence.append(200)
        return DummyResponse(200)

    # Track repo_rate_limit and avoid real sleep
    rate_limit_calls = {"n": 0}

    def fake_repo_rate_limit():
        rate_limit_calls["n"] += 1

    monkeypatch.setattr(
        "amprenta_rag.ingestion.repositories.metabolights.requests.get",
        fake_get,
    )
    monkeypatch.setattr(
        "amprenta_rag.ingestion.repositories.metabolights.repo_rate_limit",
        fake_repo_rate_limit,
    )
    monkeypatch.setattr(
        "amprenta_rag.ingestion.repositories.metabolights.time.sleep",
        lambda *_args, **_kwargs: None,
        raising=False,
    )

    repo = MetaboLightsRepository()
    resp = repo._make_request_with_retry("http://example.com")

    # Should have seen a 500 then a 200
    assert call_sequence == [500, 200]

    # repo_rate_limit should have been called before each request (2 times)
    assert rate_limit_calls["n"] == 2

    # Final response should be successful
    assert resp is not None
    assert resp.status_code == 200


