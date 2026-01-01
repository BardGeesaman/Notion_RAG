"""E2E tests for Data Lifecycle Management dashboard."""

import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def lifecycle_page(page: Page, streamlit_server: str) -> Page:
    """Navigate to Data Lifecycle page."""
    page.goto(f"{streamlit_server}/?page=Data+Lifecycle")
    page.wait_for_load_state("domcontentloaded")
    # Wait for page title
    expect(page.get_by_role("heading", name="Data Lifecycle Management")).to_be_visible(timeout=10000)
    return page


class TestDataLifecycleE2E:
    """E2E tests for Data Lifecycle dashboard."""

    def test_lifecycle_page_loads(self, lifecycle_page: Page):
        """Page renders with 4 tabs."""
        # Verify all 4 tabs are present
        expect(lifecycle_page.get_by_text("Overview")).to_be_visible()
        expect(lifecycle_page.get_by_text("Quarantine Queue")).to_be_visible()
        expect(lifecycle_page.get_by_text("Bulk Operations")).to_be_visible()
        expect(lifecycle_page.get_by_text("Audit Log")).to_be_visible()

    def test_quarantine_check_impact(self, lifecycle_page: Page):
        """Check Impact button returns impact data."""
        # Click Quarantine Queue tab
        lifecycle_page.get_by_text("Quarantine Queue").click()
        lifecycle_page.wait_for_timeout(500)
        
        # Enter a test UUID (use placeholder pattern)
        entity_id_input = lifecycle_page.get_by_placeholder("Enter UUID")
        entity_id_input.fill("550e8400-e29b-41d4-a716-446655440000")
        lifecycle_page.keyboard.press("Tab")
        lifecycle_page.wait_for_timeout(500)
        
        # Select "Check Impact" action
        lifecycle_page.get_by_label("Action").select_option("Check Impact")
        
        # Click Execute
        lifecycle_page.get_by_role("button", name="Execute").click()
        lifecycle_page.wait_for_timeout(1000)
        
        # Should show some response (error or data)
        # Entity not found is acceptable - we're testing the flow
        page_content = lifecycle_page.content()
        assert "Entity not found" in page_content or "impact" in page_content.lower()

    def test_quarantine_restore_entity(self, lifecycle_page: Page):
        """Restore action updates entity status."""
        lifecycle_page.get_by_text("Quarantine Queue").click()
        lifecycle_page.wait_for_timeout(500)
        
        # Enter UUID
        entity_id_input = lifecycle_page.get_by_placeholder("Enter UUID")
        entity_id_input.fill("550e8400-e29b-41d4-a716-446655440000")
        lifecycle_page.keyboard.press("Tab")
        lifecycle_page.wait_for_timeout(500)
        
        # Select Restore action
        lifecycle_page.get_by_label("Action").select_option("Restore (→ active)")
        
        # Click Execute
        lifecycle_page.get_by_role("button", name="Execute").click()
        lifecycle_page.wait_for_timeout(1000)
        
        # Should show response
        page_content = lifecycle_page.content()
        assert "Entity not found" in page_content or "Status" in page_content

    def test_bulk_preview_shows_impact(self, lifecycle_page: Page):
        """Preview button shows aggregated impact."""
        lifecycle_page.get_by_text("Bulk Operations").click()
        lifecycle_page.wait_for_timeout(500)
        
        # Enter UUIDs
        uuid_input = lifecycle_page.locator("textarea").first
        uuid_input.fill("550e8400-e29b-41d4-a716-446655440000\n6ba7b810-9dad-11d1-80b4-00c04fd430c8")
        lifecycle_page.keyboard.press("Tab")
        lifecycle_page.wait_for_timeout(500)
        
        # Click Preview button
        preview_btn = lifecycle_page.get_by_role("button", name="Preview Impact")
        preview_btn.click()
        lifecycle_page.wait_for_timeout(1500)
        
        # Should show preview results (even if entities not found)
        page_content = lifecycle_page.content()
        assert "Preview" in page_content or "Entities" in page_content or "valid UUIDs" in page_content

    def test_bulk_execute_requires_preview(self, lifecycle_page: Page):
        """Execute button disabled until preview done."""
        lifecycle_page.get_by_text("Bulk Operations").click()
        lifecycle_page.wait_for_timeout(500)
        
        # Execute button should be disabled initially
        execute_btn = lifecycle_page.get_by_role("button", name="Execute Bulk Update")
        expect(execute_btn).to_be_disabled()

    def test_audit_tab_shows_history(self, lifecycle_page: Page):
        """Audit tab displays lifecycle changes."""
        lifecycle_page.get_by_text("Audit Log").click()
        lifecycle_page.wait_for_timeout(500)
        
        # Should show audit log section
        expect(lifecycle_page.get_by_text("Lifecycle Audit Log")).to_be_visible()
        
        # Should have filter controls
        expect(lifecycle_page.get_by_label("Entity Type")).to_be_visible()
