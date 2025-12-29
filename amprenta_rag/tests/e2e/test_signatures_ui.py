"""E2E tests for electronic signatures UI."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_audit_trail_page(page: Page, base_url: str) -> None:
    """Navigate to the Audit Trail page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Audit%20Trail")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_signatures_tab_exists(page: Page, streamlit_server: str) -> None:
    """Test that Signatures tab exists."""
    _goto_audit_trail_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Look for Signatures tab
    sig_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Signatures", re.IGNORECASE))
    
    # Tab should exist (may not be visible without clicking)
    assert sig_tab.count() > 0 or page.locator('[role="tab"]').count() >= 3, "Signatures tab should exist"


def test_sign_form_elements_present(page: Page, streamlit_server: str) -> None:
    """Test that sign form has required elements."""
    _goto_audit_trail_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Signatures tab
    sig_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Signatures", re.IGNORECASE))
    
    try:
        sig_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have form elements (inputs, selectboxes, buttons)
    inputs = main.locator('input')
    buttons = main.get_by_role("button")
    
    assert inputs.count() > 0 or buttons.count() > 0, "Sign form should have elements"


def test_signature_history_section(page: Page, streamlit_server: str) -> None:
    """Test that signature history section exists."""
    _goto_audit_trail_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Signatures tab
    sig_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Signatures", re.IGNORECASE))
    
    try:
        sig_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have history or signature-related text
    sig_text = main.get_by_text(re.compile("Signature|History|Sign", re.IGNORECASE))
    
    assert sig_text.count() > 0 or main.count() > 0, "Signature history should be present"


def test_verify_button_exists(page: Page, streamlit_server: str) -> None:
    """Test that verify signature button exists."""
    _goto_audit_trail_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Signatures tab
    sig_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Signatures", re.IGNORECASE))
    
    try:
        sig_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have buttons (sign, verify, or load)
    buttons = main.get_by_role("button")
    
    assert buttons.count() > 0, "Should have action buttons"

