"""E2E tests for Job Queue Management page."""

import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def job_queue_page(page: Page, streamlit_server: str) -> Page:
    """Navigate to job queue page and wait for load."""
    page.goto(f"{streamlit_server}/?page=Job%20Queue")
    page.wait_for_load_state("domcontentloaded")
    
    # Wait for Streamlit to finish loading
    page.wait_for_selector('[data-testid="stMainBlockContainer"]', timeout=10000)
    return page


class TestJobQueuePage:
    """Test job queue management page functionality."""
    
    def test_page_loads(self, job_queue_page: Page) -> None:
        """Test that page renders with correct title."""
        # Look for the main title
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        expect(main.get_by_text("Job Queue Management")).to_be_visible(timeout=10000)

    def test_tabs_present(self, job_queue_page: Page) -> None:
        """Test that all 3 tabs are visible."""
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        
        # Check for tab labels
        expect(main.get_by_text("🔄 Active Jobs")).to_be_visible()
        expect(main.get_by_text("📜 Job History")).to_be_visible()
        expect(main.get_by_text("📊 Statistics")).to_be_visible()

    def test_filter_dropdowns_present(self, job_queue_page: Page) -> None:
        """Test that type and status filter dropdowns exist in Active Jobs tab."""
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        
        # Should be on Active Jobs tab by default, look for filter dropdowns
        # Streamlit selectboxes have specific data-testid patterns
        job_type_selectbox = main.locator('[data-testid="stSelectbox"]').first
        expect(job_type_selectbox).to_be_visible()
        
        # Look for multiple selectboxes (job type and status filters)
        selectboxes = main.locator('[data-testid="stSelectbox"]')
        expect(selectboxes).to_have_count_greater_than_or_equal(2)

    def test_refresh_button_present(self, job_queue_page: Page) -> None:
        """Test that refresh button is clickable."""
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        
        # Look for refresh button (may have emoji or text)
        refresh_btn = main.get_by_role("button").filter(has_text="Refresh")
        expect(refresh_btn).to_be_visible()

    def test_flower_link_present(self, job_queue_page: Page) -> None:
        """Test that Statistics tab has Flower dashboard link."""
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        
        # Click on Statistics tab
        statistics_tab = main.get_by_text("📊 Statistics")
        statistics_tab.click()
        
        # Wait for tab content to load and look for Flower link
        expect(main.get_by_text("Flower Dashboard")).to_be_visible(timeout=5000)

    def test_auto_refresh_checkbox_present(self, job_queue_page: Page) -> None:
        """Test that auto-refresh checkbox is present in Active Jobs tab."""
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        
        # Look for auto-refresh checkbox
        auto_refresh_checkbox = main.locator('[data-testid="stCheckbox"]').filter(has_text="Auto-refresh")
        expect(auto_refresh_checkbox).to_be_visible()

    def test_job_history_pagination_present(self, job_queue_page: Page) -> None:
        """Test that Job History tab has pagination controls."""
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        
        # Click on Job History tab
        history_tab = main.get_by_text("📜 Job History")
        history_tab.click()
        
        # Look for pagination controls (page size selectbox and page number input)
        expect(main.get_by_text("Page Size")).to_be_visible(timeout=5000)
        expect(main.get_by_text("Page")).to_be_visible()

    def test_statistics_metrics_present(self, job_queue_page: Page) -> None:
        """Test that Statistics tab shows metrics."""
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        
        # Click on Statistics tab
        statistics_tab = main.get_by_text("📊 Statistics")
        statistics_tab.click()
        
        # Look for metric labels (these are shown even when no data is available)
        expect(main.get_by_text("Total Jobs")).to_be_visible(timeout=5000)
        expect(main.get_by_text("Running")).to_be_visible()
        expect(main.get_by_text("Failed Today")).to_be_visible()

    def test_export_csv_button_in_history(self, job_queue_page: Page) -> None:
        """Test that Job History tab has export CSV button."""
        main = job_queue_page.locator('[data-testid="stMainBlockContainer"]').first
        
        # Click on Job History tab
        history_tab = main.get_by_text("📜 Job History")
        history_tab.click()
        
        # Look for export CSV button
        export_btn = main.get_by_role("button").filter(has_text="Export CSV")
        expect(export_btn).to_be_visible(timeout=5000)
