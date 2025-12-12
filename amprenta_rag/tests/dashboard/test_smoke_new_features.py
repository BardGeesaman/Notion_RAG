"""Smoke tests for new dashboard features."""
import pytest
from playwright.sync_api import Page, expect


@pytest.mark.parametrize("page_name,expected_text,expander", [
    ("Compare", "Compare", "📚 Other Pages"),
    ("Timeline", "Timeline", "📚 Other Pages"),
    ("Data Quality", "Quality", "📚 Other Pages"),
    ("Workflows", "Workflow", "⚙️ Admin"),
    ("Literature Analysis", "Literature", "📊 Analysis"),
    ("Candidate Selection", "Candidate", "📊 Analysis"),
    ("Email Settings", "Email", "📚 Other Pages"),
    ("Data Lineage", "Lineage", "📚 Other Pages"),
    ("Experiments", "Experiment", "🔍 Discovery"),
    ("Chemistry", "Chemistry", "📊 Analysis"),
    ("RAG Query", "RAG", "📚 Other Pages"),
])
def test_page_loads(page: Page, base_url: str, page_name: str, expected_text: str, expander: str) -> None:
    """
    Test that new feature pages load correctly.
    
    Args:
        page: Playwright page fixture
        base_url: Base URL fixture from pytest-base-url
        page_name: Name of the page to navigate to
        expected_text: Expected text on the page (heading or key text)
        expander: Name of the expander that contains this page
    """
    # Navigate to dashboard
    page.goto(base_url)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(3000)  # Wait for Streamlit to fully load
    
    # Click the correct expander to reveal the page button
    try:
        expander_element = page.locator(f"text={expander}").first
        if expander_element.is_visible():
            expander_element.click()
            page.wait_for_timeout(500)  # Wait for expander to open
    except Exception as e:
        pytest.fail(f"Could not find or click expander '{expander}': {e}")
    
    # Click the page button
    try:
        page.get_by_role("button", name=page_name, exact=True).click()
    except Exception:
        # Try without exact match if exact match fails
        button = page.get_by_role("button", name=page_name).first
        if button.is_visible():
            button.click()
        else:
            pytest.fail(f"Could not find button for page: {page_name}")
    
    page.wait_for_timeout(2000)  # Wait for page to load
    
    # Assert expected text is present (more flexible than exact heading match)
    expect(page.locator(f"text={expected_text}").first).to_be_attached(timeout=10000)
    
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
