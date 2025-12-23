"""
E2E tests for DVC auto-versioning UI + backend imports.
"""

from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

import pytest
import requests
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


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
def dataset_details_streamlit_server(fastapi_server: str):
    """Start a Streamlit server specifically for dataset_details.py page."""
    port = 8503
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

    env = os.environ.copy()
    env["PYTHONPATH"] = project_root
    env["DISABLE_AUTH"] = "1"
    env["API_URL"] = fastapi_server

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


def _get_any_dataset_id(api_base: str) -> str | None:
    try:
        resp = requests.get(f"{api_base}/api/v1/datasets?limit=1", timeout=10)
        resp.raise_for_status()
        data = resp.json()
        if isinstance(data, list) and data:
            return str(data[0].get("id"))
    except Exception:
        return None
    return None


class TestDVCE2E:
    def test_dataset_detail_has_version_history_tab(
        self, page: Page, fastapi_server: str, dataset_details_streamlit_server: str
    ):
        ds_id = _get_any_dataset_id(fastapi_server)
        if not ds_id:
            pytest.skip("No datasets available to load Dataset Details page.")

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


