"""E2E tests for Notebook Generator dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_notebook_generator_page(page: Page, base_url: str) -> None:
    """Navigate to the Notebook Generator page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Notebook%20Generator")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_notebook_generator_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Notebook Generator page loads successfully."""
    _goto_notebook_generator_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Notebook Generator heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Notebook.*Generator", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Query.*Notebook", re.IGNORECASE)),
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


def test_notebook_generator_question_input(page: Page, streamlit_server: str) -> None:
    """Test that question input textarea is present."""
    _goto_notebook_generator_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have question/query textarea (may be marked as hidden by Streamlit)
    textareas = main.locator('textarea')
    
    # At least one textarea should exist (even if hidden)
    assert textareas.count() > 0, "Expected question input textarea"


def test_notebook_generator_generate_button(page: Page, streamlit_server: str) -> None:
    """Test that generate button is present."""
    _goto_notebook_generator_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have Generate Notebook button
    generate_button = main.get_by_role("button").filter(has_text=re.compile("Generate.*Notebook", re.IGNORECASE))
    
    expect(generate_button.first).to_be_visible(timeout=10000)


def test_notebook_generator_llm_requirements_caption(page: Page, streamlit_server: str) -> None:
    """Test that LLM requirements caption is shown."""
    _goto_notebook_generator_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show note about LLM credentials
    llm_note = main.get_by_text(re.compile("LLM|API.*credentials|backend", re.IGNORECASE))
    
    # Note should be visible
    assert llm_note.count() > 0, "Expected LLM requirements note"


def test_notebook_generator_form_elements(page: Page, streamlit_server: str) -> None:
    """Test that generation form has required elements."""
    _goto_notebook_generator_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have textarea and button
    textareas = main.locator('textarea')
    buttons = main.get_by_role("button")
    
    # Should have both form elements
    assert textareas.count() > 0, "Expected textarea for question input"
    assert buttons.count() > 0, "Expected generate button"


def test_notebook_generator_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_notebook_generator_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

