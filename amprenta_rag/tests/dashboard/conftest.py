from __future__ import annotations

from contextlib import contextmanager
from typing import Dict, Any

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest
import requests


def _kill_port(port: int) -> None:
    """Best-effort kill any process listening on a TCP port (avoids E2E port conflicts)."""
    try:
        res = subprocess.run(
            ["lsof", "-ti", f":{port}"],
            capture_output=True,
            text=True,
            check=False,
        )
        pids = [p.strip() for p in (res.stdout or "").splitlines() if p.strip()]
        for pid in pids:
            subprocess.run(["kill", "-9", pid], capture_output=True, text=True, check=False)
        if pids:
            time.sleep(0.5)
    except Exception:
        # Best-effort; don't fail the test suite if we can't kill the port.
        return


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
def fastapi_server():
    """Start FastAPI server for E2E tests."""
    port = 8000
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

    env = os.environ.copy()
    env["PYTHONPATH"] = project_root
    env["DISABLE_AUTH"] = "1"
    env.setdefault("AMPRENTA_DEBUG_NAV", "1")

    log_path = Path("/tmp/fastapi_e2e.log")
    log_f = log_path.open("w", encoding="utf-8")

    _kill_port(port)
    proc = subprocess.Popen(
        [sys.executable, "-m", "uvicorn", "amprenta_rag.api.main:app", f"--port={port}", "--host=0.0.0.0"],
        cwd=project_root,
        env=env,
        stdout=log_f,
        stderr=log_f,
    )

    # Wait for server to be ready
    max_wait = 30
    for _ in range(max_wait):
        try:
            resp = requests.get(f"http://localhost:{port}/health", timeout=1)
            if resp.status_code == 200:
                break
        except requests.exceptions.RequestException:
            pass
        time.sleep(1)
    else:
        try:
            log_f.flush()
        except Exception:
            pass
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
        try:
            log_f.close()
        except Exception:
            pass
        tail = ""
        try:
            tail = log_path.read_text(encoding="utf-8")[-4000:]
        except Exception:
            pass
        raise RuntimeError(f"FastAPI server failed to start on port {port}\n{tail}")

    try:
        yield f"http://localhost:{port}"
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
        try:
            log_f.close()
        except Exception:
            pass


@pytest.fixture(scope="session")
def streamlit_server(fastapi_server):
    """Start Streamlit server for E2E tests, stop when done."""
    port = 8502
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

    env = os.environ.copy()
    env["PYTHONPATH"] = project_root
    env["DISABLE_AUTH"] = "1"
    env["AMPRENTA_DEBUG_NAV"] = "1"
    env["API_URL"] = fastapi_server

    log_path = Path("/tmp/streamlit_e2e.log")
    log_f = log_path.open("w", encoding="utf-8")

    _kill_port(port)
    proc = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", "scripts/dashboard/app.py", f"--server.port={port}", "--server.headless=true"],
        cwd=project_root,
        env=env,
        stdout=log_f,
        stderr=log_f,
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
        try:
            log_f.flush()
        except Exception:
            pass
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
        try:
            log_f.close()
        except Exception:
            pass
        tail = ""
        try:
            tail = log_path.read_text(encoding="utf-8")[-4000:]
        except Exception:
            pass
        raise RuntimeError(f"Streamlit server failed to start on port {port}\n{tail}")

    try:
        yield f"http://localhost:{port}"
    finally:
        # Cleanup
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
        try:
            log_f.close()
        except Exception:
            pass


@pytest.fixture(scope="session")
def browser_type_launch_args(browser_type_launch_args):
    """
    Optionally force Playwright to run headed for E2E debugging.

    Default behavior should remain whatever pytest-playwright chooses (typically
    headless in CI). To force headed, set `E2E_HEADED=1`.
    """
    if os.getenv("E2E_HEADED") == "1":
        return {**browser_type_launch_args, "headless": False}
    return browser_type_launch_args

