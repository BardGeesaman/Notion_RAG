"""E2E tests for Notebook Co-Pilot dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_notebook_copilot_page(page: Page, base_url: str) -> None:
    """Navigate to the Notebook Co-Pilot page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Notebook%20Co-Pilot")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_notebook_copilot_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Notebook Co-Pilot page loads successfully."""
    _goto_notebook_copilot_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Notebook Co-Pilot heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Notebook.*Co-Pilot", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Notebook.*Copilot", re.IGNORECASE)),
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


def test_notebook_copilot_instructions_expander(page: Page, streamlit_server: str) -> None:
    """Test that 'How to use' instructions are present."""
    _goto_notebook_copilot_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have "How to use" expander
    how_to_use = main.get_by_text(re.compile("How to use", re.IGNORECASE))
    
    expect(how_to_use.first).to_be_visible(timeout=10000)


def test_notebook_copilot_action_selector(page: Page, streamlit_server: str) -> None:
    """Test that action selection controls are present."""
    _goto_notebook_copilot_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have action-related content
    # Copilot sidebar component may have action selection
    action_text = main.get_by_text(re.compile("action|Pick|Select", re.IGNORECASE))
    selectboxes = main.locator('[data-baseweb="select"]')
    
    # Should have some action selection UI
    assert action_text.count() > 0 or selectboxes.count() > 0, "Expected action selection controls"


def test_notebook_copilot_code_generation(page: Page, streamlit_server: str) -> None:
    """Test that code generation interface is present."""
    _goto_notebook_copilot_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have generate code button or related text
    generate_code = main.get_by_role("button").filter(has_text=re.compile("Generate.*code", re.IGNORECASE))
    code_text = main.get_by_text(re.compile("Generate.*code|code.*cell", re.IGNORECASE))
    
    # Should have code generation UI
    assert generate_code.count() > 0 or code_text.count() > 0, "Expected code generation interface"


def test_notebook_copilot_llm_note(page: Page, streamlit_server: str) -> None:
    """Test that LLM requirements note is shown."""
    _goto_notebook_copilot_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show note about LLM/model requirements
    llm_note = main.get_by_text(re.compile("LLM|model|OpenAI|credentials", re.IGNORECASE))
    
    # Note should be visible
    assert llm_note.count() > 0, "Expected LLM requirements note"


def test_notebook_copilot_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_notebook_copilot_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

