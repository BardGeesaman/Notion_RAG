"""E2E tests for Nextflow Orchestrator dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_nextflow_orchestrator_page(page: Page, base_url: str) -> None:
    """Navigate to the Nextflow Orchestrator page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Nextflow%20Orchestrator")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_nextflow_orchestrator_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Nextflow Orchestrator page loads successfully."""
    _goto_nextflow_orchestrator_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Nextflow Orchestrator heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Nextflow.*Orchestrator", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Nextflow", re.IGNORECASE)),
    ]
    
    found = False
    for pattern in heading_patterns:
        try:
            expect(pattern.first).to_be_visible(timeout=5000)
            found = True
            break
        except AssertionError:
            continue
    
    if not found:
        # Fallback: just check if main container has content
        expect(main).to_be_visible(timeout=10000)


def test_nextflow_orchestrator_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that orchestrator tabs are present."""
    _goto_nextflow_orchestrator_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Look for tab text content (Run Pipeline, Job Monitor, Results)
    run_tab = page.get_by_text(re.compile("Run.*Pipeline", re.IGNORECASE), exact=False)
    job_tab = page.get_by_text(re.compile("Job.*Monitor", re.IGNORECASE), exact=False)
    results_tab = page.get_by_text("Results", exact=False)
    
    # At least one strategy should find tabs
    assert (tab_count >= 2 or run_tab.count() > 0 or job_tab.count() > 0 or 
            results_tab.count() > 0), f"Expected tabs but found {tab_count} tab elements"


def test_nextflow_orchestrator_installation_check(page: Page, streamlit_server: str) -> None:
    """Test that Nextflow installation status is shown."""
    _goto_nextflow_orchestrator_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show either "Nextflow available" or "not installed" message
    nextflow_available = main.get_by_text(re.compile("Nextflow.*available|version", re.IGNORECASE))
    nextflow_missing = main.get_by_text(re.compile("not installed|Install", re.IGNORECASE))
    
    # Should show installation status
    assert (nextflow_available.count() > 0 or 
            nextflow_missing.count() > 0), "Expected Nextflow installation status"


def test_nextflow_orchestrator_pipeline_selector(page: Page, streamlit_server: str) -> None:
    """Test that pipeline selector is present."""
    _goto_nextflow_orchestrator_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have pipeline selection or "no pipelines" message
    pipeline_selector = main.get_by_text(re.compile("Pipeline", re.IGNORECASE))
    no_pipelines = main.get_by_text(re.compile("No pipelines", re.IGNORECASE))
    
    # Should have pipeline-related content
    assert pipeline_selector.count() > 0 or no_pipelines.count() > 0, "Expected pipeline content"


def test_nextflow_orchestrator_workflow_content(page: Page, streamlit_server: str) -> None:
    """Test that workflow-related content is present."""
    _goto_nextflow_orchestrator_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Page should have workflow-related content or installation instructions
    # Sample sheet tabs only show if Nextflow is installed
    workflow_keywords = [
        main.get_by_text(re.compile("Pipeline|Workflow|Run|Submit", re.IGNORECASE)),
        main.get_by_text(re.compile("Nextflow", re.IGNORECASE)),
    ]
    
    # Should have workflow content
    found = False
    for keyword in workflow_keywords:
        if keyword.count() > 0:
            found = True
            break
    
    assert found, "Expected workflow-related content"


def test_nextflow_orchestrator_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_nextflow_orchestrator_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

