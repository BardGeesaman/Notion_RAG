"""Playwright test utilities with error detection."""

def check_page_errors(page):
    """Raise exception if page shows Streamlit error."""
    if page.locator("text=Error rendering").count() > 0:
        error_text = page.locator(".stException").first.text_content() or "Unknown error"
        raise Exception(f"Page error: {error_text[:200]}")
