from __future__ import annotations

from contextlib import contextmanager
from typing import Dict, Any

import os
import subprocess
import time

import pytest
import requests


@pytest.fixture
def test_program() -> Dict[str, Any]:
    return {"id": "prog-uuid-1", "name": "ALS Program", "disease": ["ALS"]}


@pytest.fixture
def test_dataset() -> Dict[str, Any]:
    return {"id": "ds-uuid-1", "name": "ALS CSF Lipidomics", "omics_type": "lipidomics"}


@pytest.fixture
def test_signature() -> Dict[str, Any]:
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


@pytest.fixture(scope="session")
def streamlit_server():
    """Start Streamlit server for E2E tests, stop when done."""
    port = 8502
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

    env = os.environ.copy()
    env["PYTHONPATH"] = project_root

    proc = subprocess.Popen(
        ["streamlit", "run", "scripts/dashboard/app.py", f"--server.port={port}", "--server.headless=true"],
        cwd=project_root,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Wait for server to be ready
    max_wait = 30
    for _ in range(max_wait):
        try:
            resp = requests.get(f"http://localhost:{port}/_stcore/health", timeout=1)
            if resp.status_code == 200:
                break
        except requests.exceptions.RequestException:
            pass
        time.sleep(1)
    else:
        proc.terminate()
        raise RuntimeError(f"Streamlit server failed to start on port {port}")

    try:
        yield f"http://localhost:{port}"
    finally:
        # Cleanup
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()

