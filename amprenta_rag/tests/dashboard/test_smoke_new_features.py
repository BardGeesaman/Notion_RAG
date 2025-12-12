"""Smoke tests for new dashboard features."""
import pytest
from playwright.sync_api import Page, expect


@pytest.mark.parametrize("page_name,expected_heading", [
    ("Compare", "Compare"),
    ("Timeline", "Activity Timeline"),
    ("Data Quality", "Data Quality"),
])
def test_page_loads(page: Page, base_url: str, page_name: str, expected_heading: str) -> None:
    """
    Test that new feature pages load correctly.
    
    Args:
        page: Playwright page fixture
        base_url: Base URL fixture from pytest-base-url
        page_name: Name of the page to navigate to
        expected_heading: Expected heading text on the page
    """
    # Navigate to dashboard
    page.goto(base_url)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(3000)  # Wait for Streamlit to fully load
    
    # Expand "Other Pages" expander
    page.locator("text=Other Pages").first.click()
    page.wait_for_timeout(500)  # Wait for expander to open
    
    # Click the page button
    page.get_by_role("button", name=page_name, exact=True).click()
    page.wait_for_timeout(2000)  # Wait for page to load
    
    # Assert main heading is present
    expect(page.locator(f"text={expected_heading}").first).to_be_attached(timeout=10000)
    
    # Check for error alerts - allow info alerts but flag if contains "Error"
    alerts = page.locator("div[role='alert']").all()
    for alert in alerts:
        alert_text = alert.text_content() or ""
        if "error" in alert_text.lower() and "❌" in alert_text:
            pytest.fail(f"Error alert found on {page_name} page: {alert_text}")
    
    # Also check for Streamlit error messages
    error_elements = page.locator("text=/error/i").all()
    for error_elem in error_elements:
        # Skip if it's just part of normal text
        parent = error_elem.locator("..")
        if parent.get_attribute("role") == "alert":
            pytest.fail(f"Error message found on {page_name} page: {error_elem.text_content()}")
