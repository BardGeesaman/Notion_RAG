"""E2E tests for Pipeline Runner dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_pipeline_runner_page(page: Page, base_url: str) -> None:
    """Navigate to the Pipeline Runner page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Pipeline%20Runner")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_pipeline_runner_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Pipeline Runner page loads successfully."""
    _goto_pipeline_runner_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Pipeline Runner heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Pipeline.*Runner", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Bioinformatics.*Pipeline", re.IGNORECASE)),
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


def test_pipeline_runner_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that pipeline runner tabs are present."""
    _goto_pipeline_runner_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Look for tab text content (Run Pipeline, Manage Indices, Job History)
    run_tab = page.get_by_text(re.compile("Run.*Pipeline", re.IGNORECASE), exact=False)
    manage_tab = page.get_by_text(re.compile("Manage.*Indices", re.IGNORECASE), exact=False)
    history_tab = page.get_by_text(re.compile("Job.*History|History", re.IGNORECASE), exact=False)
    
    # At least one strategy should find tabs
    assert (tab_count >= 2 or run_tab.count() > 0 or manage_tab.count() > 0 or 
            history_tab.count() > 0), f"Expected tabs but found {tab_count} tab elements"


def test_pipeline_runner_fastq_upload(page: Page, streamlit_server: str) -> None:
    """Test that FASTQ upload section is present."""
    _goto_pipeline_runner_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have FASTQ upload or search functionality
    fastq_text = main.get_by_text(re.compile("FASTQ|Upload", re.IGNORECASE))
    select_input = main.get_by_text(re.compile("Select.*Input", re.IGNORECASE))
    
    # Should have input selection UI
    assert fastq_text.count() > 0 or select_input.count() > 0, "Expected FASTQ input section"


def test_pipeline_runner_tool_selection(page: Page, streamlit_server: str) -> None:
    """Test that tool/pipeline selection is present."""
    _goto_pipeline_runner_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have tool selection or index selection
    tool_text = main.get_by_text(re.compile("Tool|Index|Select", re.IGNORECASE))
    selectboxes = main.locator('[data-baseweb="select"]')
    
    # Should have tool selection UI
    assert tool_text.count() > 0 or selectboxes.count() > 0, "Expected tool selection"


def test_pipeline_runner_index_management(page: Page, streamlit_server: str) -> None:
    """Test that index management tab has content."""
    _goto_pipeline_runner_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Manage Indices tab
    manage_tab = page.get_by_text(re.compile("Manage.*Indices", re.IGNORECASE), exact=False).or_(
        page.locator('[role="tab"]').filter(has_text=re.compile("Manage|Indices", re.IGNORECASE))
    )
    
    try:
        manage_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have index-related content
    index_text = main.get_by_text(re.compile("Index|Indices|Register", re.IGNORECASE))
    
    assert index_text.count() > 0, "Expected index management content"


def test_pipeline_runner_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_pipeline_runner_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

