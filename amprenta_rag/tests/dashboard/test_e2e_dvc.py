"""
E2E tests for DVC auto-versioning UI + backend imports.
"""

from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path
from uuid import uuid4

import pytest
import requests
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


@pytest.fixture(scope="session")
def dvc_fastapi_server():
    """Start a tiny in-memory FastAPI server for Dataset Details E2E."""
    import threading

    import uvicorn
    from fastapi import FastAPI, HTTPException

    port = 8005
    _kill_port(port)

    app = FastAPI()
    store = {}

    @app.get("/health")
    def _health():
        return {"ok": True}

    @app.post("/api/v1/datasets/")
    def _create_dataset(payload: dict):
        ds_id = str(uuid4())
        now = "2025-01-01T00:00:00Z"
        store[ds_id] = {
            "id": ds_id,
            "name": payload.get("name") or "DVC Test Dataset",
            "omics_type": payload.get("omics_type") or "transcriptomics",
            "description": payload.get("description") or "Test dataset for E2E",
            "file_paths": payload.get("file_paths") or [],
            "file_urls": payload.get("file_urls") or [],
            "organism": payload.get("organism") or [],
            "sample_type": payload.get("sample_type") or [],
            "disease": payload.get("disease") or [],
            "external_ids": payload.get("external_ids"),
            "created_at": now,
            "updated_at": now,
            "version": 1,
            "ingestion_status": "pending",
            "feature_ids": {},
            "dvc_version": None,
            "dvc_metadata": {"versions": []},
            "dvc_pushed": False,
        }
        return store[ds_id]

    @app.get("/api/v1/datasets/{dataset_id}")
    def _get_dataset(dataset_id: str):
        ds = store.get(dataset_id)
        if not ds:
            raise HTTPException(status_code=404, detail="Dataset not found")
        return ds

    @app.patch("/api/v1/datasets/{dataset_id}")
    def _patch_dataset(dataset_id: str, payload: dict):
        ds = store.get(dataset_id)
        if not ds:
            raise HTTPException(status_code=404, detail="Dataset not found")
        ds.update(payload)
        ds["updated_at"] = "2025-01-01T00:00:01Z"
        ds["version"] = int(ds.get("version") or 1) + 1
        return ds

    @app.delete("/api/v1/datasets/{dataset_id}")
    def _delete_dataset(dataset_id: str):
        store.pop(dataset_id, None)
        return {"ok": True}

    config = uvicorn.Config(app, host="127.0.0.1", port=port, log_level="warning")
    server = uvicorn.Server(config)
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()

    # wait for server
    for _ in range(30):
        try:
            r = requests.get(f"http://127.0.0.1:{port}/health", timeout=1)
            if r.status_code == 200:
                break
        except Exception:
            pass
        time.sleep(0.25)
    else:
        raise RuntimeError("Failed to start dvc_fastapi_server")

    yield f"http://127.0.0.1:{port}"

    server.should_exit = True
    try:
        thread.join(timeout=3)
    except Exception:
        pass


@pytest.fixture
def test_dataset(dvc_fastapi_server: str):
    """Create a test dataset for DVC tests."""
    import httpx

    with httpx.Client(timeout=15) as client:
        resp = client.post(
            f"{dvc_fastapi_server}/api/v1/datasets/",
            json={
                "name": "DVC Test Dataset",
                "omics_type": "transcriptomics",
                "description": "Test dataset for E2E",
            },
        )
        resp.raise_for_status()
        ds_id = resp.json().get("id")
        assert ds_id

    try:
        yield str(ds_id)
    finally:
        try:
            with httpx.Client(timeout=15) as client:
                client.delete(f"{dvc_fastapi_server}/api/v1/datasets/{ds_id}")
        except Exception:
            pass


def _kill_port(port: int) -> None:
    """Best-effort kill any process listening on a TCP port."""
    try:
        res = subprocess.run(["lsof", "-ti", f":{port}"], capture_output=True, text=True, check=False)
        pids = [p.strip() for p in (res.stdout or "").splitlines() if p.strip()]
        for pid in pids:
            subprocess.run(["kill", "-9", pid], capture_output=True, text=True, check=False)
        if pids:
            time.sleep(0.5)
    except Exception:
        return


def _main_container(page: Page):
    main = page.locator('[data-testid="stMainBlockContainer"]')
    if main.count() > 0:
        return main.first
    return page.locator('[data-testid="stAppViewContainer"]').first


@pytest.fixture(scope="session")
def dataset_details_streamlit_server(dvc_fastapi_server: str):
    """Start a Streamlit server specifically for dataset_details.py page."""
    port = 8503
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

    env = os.environ.copy()
    env["PYTHONPATH"] = project_root
    env["DISABLE_AUTH"] = "1"
    env["API_URL"] = dvc_fastapi_server

    log_path = Path("/tmp/streamlit_dataset_details_e2e.log")
    log_f = log_path.open("w", encoding="utf-8")

    _kill_port(port)
    proc = subprocess.Popen(
        ["streamlit", "run", "scripts/dashboard/pages/dataset_details.py", f"--server.port={port}", "--server.headless=true"],
        cwd=project_root,
        env=env,
        stdout=log_f,
        stderr=log_f,
    )

    # Wait for server to be ready
    for _ in range(30):
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
        tail = ""
        try:
            tail = log_path.read_text(encoding="utf-8")[-4000:]
        except Exception:
            pass
        raise RuntimeError(f"Streamlit dataset_details server failed to start on port {port}\n{tail}")

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


class TestDVCE2E:
    def test_dataset_detail_has_version_history_tab(
        self, page: Page, dataset_details_streamlit_server: str, test_dataset: str
    ):
        ds_id = test_dataset

        page.goto(dataset_details_streamlit_server)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(1500)

        # Load dataset by ID
        page.locator('input[aria-label="Dataset ID"]').fill(ds_id)
        page.keyboard.press("Tab")
        page.get_by_role("button", name="Load Dataset").click()
        page.wait_for_timeout(2000)

        main = _main_container(page)
        expect(main.get_by_role("tab", name="Version History")).to_be_visible(timeout=20000)

    def test_dvc_manager_imports(self):
        from amprenta_rag.ingestion import dvc_manager  # noqa: F401


