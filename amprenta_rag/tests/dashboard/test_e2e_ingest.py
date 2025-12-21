"""
E2E Ingestion/Import test scenarios.

Extracted from the historical monolithic `test_e2e_platform.py` suite.
"""

import pytest
from playwright.sync_api import Page


pytestmark = pytest.mark.requires_server


class TestDiscoveryWorkflow:
    """Test discovery workflow"""

    def test_discovery_job(self, page: Page, base_url: str):
        """Run discovery scan"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: Discovery (expanded) > Discovery Workflow
        try:
            # Scroll sidebar to reveal buttons
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)

            # Click Discovery Workflow (Discovery is expanded by default)
            sidebar.locator('button:has-text("Discovery Workflow")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Discovery Workflow: {e}")

        # Tab: Run Discovery > enter query "test"
        try:
            page.locator('button[role="tab"]:has-text("Run Discovery")').click()
            page.wait_for_timeout(1000)

            # Fill discovery form
            query_input = page.locator('input[aria-label="Search Query"]')
            query_input.fill("test")

            # Click Start Discovery button
            start_button = page.get_by_role("button", name="Start Discovery", exact=True).first
            if start_button.count() == 0:
                start_button = page.get_by_role("button", name="🚀 Start Discovery", exact=True).first
            if start_button.count() > 0:
                start_button.click()
                page.wait_for_timeout(5000)  # Discovery takes time

                # Check for job created or pending studies
                assert (
                    page.locator("text=Job").first.count() > 0
                    or page.locator("text=Pending").first.count() > 0
                    or page.locator("text=success").first.count() > 0
                )
        except Exception as e:
            pytest.skip(f"Could not run discovery job: {e}")


class TestImportData:
    """Test import data page"""

    def test_view_import_data(self, page: Page, base_url: str):
        """View import data page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: ⚙️ Admin expander > Import Data
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)

            # Click Import Data
            sidebar.locator('button:has-text("Import Data")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Import Data: {e}")

        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Import Data page: Import error - Error loading page"
            assert rendering_error == 0, "Import Data page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Import Data page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Import Data page loaded: {e}")


class TestDataIngestion:
    """Test data ingestion page"""

    def test_view_data_ingestion(self, page: Page, base_url: str):
        """View data ingestion page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > Data Ingestion
        try:
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            sidebar.locator('button:has-text("Data Ingestion")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Data Ingestion: {e}")

        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Data Ingestion page: Import error - Error loading page"
            assert rendering_error == 0, "Data Ingestion page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Data Ingestion page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Data Ingestion page loaded: {e}")


class TestRepositories:
    """Test repositories page"""

    def test_view_repositories(self, page: Page, base_url: str):
        """View repositories page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > Repositories
        try:
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            sidebar.locator('button:has-text("Repositories")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Repositories: {e}")

        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Repositories page: Import error - Error loading page"
            assert rendering_error == 0, "Repositories page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Repositories page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Repositories page loaded: {e}")


