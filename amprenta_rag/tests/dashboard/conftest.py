from __future__ import annotations

from contextlib import contextmanager
from typing import Dict

import pytest


@pytest.fixture
def test_program() -> Dict[str, str]:
    return {"id": "prog-uuid-1", "name": "ALS Program", "disease": ["ALS"]}


@pytest.fixture
def test_dataset() -> Dict[str, str]:
    return {"id": "ds-uuid-1", "name": "ALS CSF Lipidomics", "omics_type": "lipidomics"}


@pytest.fixture
def test_signature() -> Dict[str, str]:
    return {"id": "sig-uuid-1", "name": "ALS-Core-6Ceramides"}


@pytest.fixture
def mock_db_session(monkeypatch):
    """
    Mock dashboard db_session only when explicitly requested by a test.

    This avoids interfering with module imports during collection.
    """

    @contextmanager
    def _dummy_db_session():
        yield None

    monkeypatch.setattr(
        "scripts.dashboard.db_session.db_session",
        _dummy_db_session,
        raising=False,
    )
    yield


# Marker for UI tests
pytestmark = pytest.mark.ui

