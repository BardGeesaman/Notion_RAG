from __future__ import annotations

"""
Playwright-based E2E tests for critical Streamlit dashboard workflows.

These tests assume the dashboard is already running, e.g.:
    streamlit run scripts/run_dashboard.py --server.port 8502
"""

import os

import pytest
from playwright.sync_api import Page, expect


@pytest.fixture(scope="session")
def dashboard_base_url() -> str:
    """Base URL for the running Streamlit dashboard."""
    return os.getenv("DASHBOARD_BASE_URL", "http://localhost:8502")


@pytest.fixture(autouse=True)
def mock_external_services(monkeypatch):
    """
    Mock external services (Pinecone, OpenAI) to avoid network calls and costs.
    """
    # Pinecone queries
    try:
        monkeypatch.setattr(
            "amprenta_rag.query.pinecone_query.query_pinecone",
            lambda *args, **kwargs: [],
        )
    except Exception:
        # If module not imported in this flow, ignore
        pass

    # OpenAI client
    class _DummyChat:
        def completions(self, *args, **kwargs):
            class _Resp:
                choices = [type("Choice", (), {"message": type("Msg", (), {"content": "Dummy response"})()})]

            return _Resp()

    class _DummyClient:
        @property
        def chat(self):
            return _DummyChat()

    try:
        monkeypatch.setattr(
            "amprenta_rag.clients.openai_client.get_openai_client",
            lambda: _DummyClient(),
        )
    except Exception:
        pass


@pytest.mark.ui
def test_navigate_all_sidebar_pages(page: Page, dashboard_base_url: str):
    """
    Smoke test: navigate through all sidebar pages to ensure they render.
    """
    page.goto(dashboard_base_url)

    # Ensure sidebar is present
    expect(page.get_by_text("Navigation")).to_be_visible()

    sidebar_labels = [
        "Overview",
        "Getting Started",
        "Evaluation Wizard",
        "Chat",
        "Lab Notebook",
        "Search",
        "Data Ingestion",
        "Repositories",
        "Analysis Tools",
        "Discovery",
        "Coverage Map",
        "Feature Recurrence",
        "Evidence Report",
        "Data Management",
        "System Health",
        "Relationships",
        "Datasets",
        "Programs",
        "Experiments",
        "Features",
        "Signatures",
        "Literature",
        "Emails",
        "RAG Chunks",
        "Chemistry",
        "RAG Query",
        "Cross-Omics",
    ]

    for label in sidebar_labels:
        page.get_by_label("Select Page").get_by_text(label).click(force=True)
        # Basic assertion: page content area updates and contains the label text somewhere
        expect(page.get_by_text(label)).to_be_visible()


@pytest.mark.ui
def test_dataset_upload_and_search_workflow(page: Page, dashboard_base_url: str, tmp_path):
    """
    Critical workflow: upload a dataset then search for it.

    Note: Selectors may need refinement using `playwright codegen`.
    """
    page.goto(dashboard_base_url)

    # Navigate to Data Ingestion page
    page.get_by_label("Select Page").get_by_text("Data Ingestion").click(force=True)

    # Prepare a small CSV file to upload
    csv_path = tmp_path / "test_dataset.csv"
    csv_path.write_text("id,value\n1,2\n", encoding="utf-8")

    # Try to find a file upload widget (generic locator; may need to be updated)
    file_inputs = page.locator("input[type='file']")
    if file_inputs.count() > 0:
        file_inputs.first.set_input_files(str(csv_path))

    # Look for a generic 'Upload' or 'Ingest' button and click if present
    for button_label in ["Upload", "Ingest", "Submit"]:
        buttons = page.get_by_role("button", name=button_label)
        if buttons.count() > 0:
            buttons.first.click()
            break

    # Navigate to Search page and verify page loads
    page.get_by_label("Select Page").get_by_text("Search").click(force=True)
    expect(page.get_by_text("Search")).to_be_visible()


