"""End-to-end tests for Review Threads functionality in Review Queue."""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


class TestReviewThreadsE2E:
    """Test Review Threads functionality in Review Queue page."""
    
    @pytest.fixture(autouse=True)
    def setup_page(self, page: Page, streamlit_server: str):
        """Navigate to Review Queue page before each test."""
        page.goto(f"{streamlit_server}/?page=Review+Queue")
        page.wait_for_load_state("domcontentloaded")
        
        # Wait for page to load
        page.wait_for_selector("[data-testid='stMainBlockContainer']", timeout=10000)
    
    def test_review_queue_page_loads(self, page: Page):
        """Test that Review Queue page renders with tabs in expanders."""
        # Check header
        expect(page.get_by_text("Review Queue")).to_be_visible()
        expect(page.get_by_text("Approve/reject notebook dashboards")).to_be_visible()
        
        # Look for pending reviews text (may show "0" or actual count)
        expect(page.get_by_text("Pending reviews:")).to_be_visible()
        
        # If there are review expanders, check for tabs
        # Note: This test may show "No pending reviews" which is fine
        expanders = page.locator("[data-testid='stExpander']")
        if expanders.count() > 0:
            # Click first expander to open it
            expanders.first.click()
            page.wait_for_timeout(1000)
            
            # Check for tabs within the expander
            expect(page.get_by_role("tab", name="📝 Review")).to_be_visible()
            expect(page.get_by_role("tab", name="💬 Discussion")).to_be_visible()
            expect(page.get_by_role("tab", name="📋 Diff")).to_be_visible()
    
    def test_discussion_tab_visible(self, page: Page):
        """Test that Discussion tab loads thread view components."""
        # Look for review expanders
        expanders = page.locator("[data-testid='stExpander']")
        
        if expanders.count() == 0:
            # If no pending reviews, we can't test the tabs
            expect(page.get_by_text("No pending reviews")).to_be_visible()
            pytest.skip("No pending reviews to test Discussion tab")
        
        # Click first expander to open it
        expanders.first.click()
        page.wait_for_timeout(1000)
        
        # Click Discussion tab
        discussion_tab = page.get_by_role("tab", name="💬 Discussion")
        expect(discussion_tab).to_be_visible()
        discussion_tab.click()
        page.wait_for_timeout(1000)
        
        # Check for thread view components
        discussion_content = page.get_by_label("💬 Discussion")
        expect(discussion_content.get_by_text("💬 Discussion Threads")).to_be_visible()
        expect(discussion_content.get_by_text("📝 Start New Discussion")).to_be_visible()
        
        # Check for thread creation form
        expect(discussion_content.get_by_placeholder("e.g., 'Question about methodology in cell 3'")).to_be_visible()
        expect(discussion_content.get_by_role("button", name="💬 Create Thread")).to_be_visible()
    
    def test_create_new_thread(self, page: Page):
        """Test creating a new discussion thread."""
        # Look for review expanders
        expanders = page.locator("[data-testid='stExpander']")
        
        if expanders.count() == 0:
            pytest.skip("No pending reviews to test thread creation")
        
        # Click first expander to open it
        expanders.first.click()
        page.wait_for_timeout(1000)
        
        # Click Discussion tab
        discussion_tab = page.get_by_role("tab", name="💬 Discussion")
        discussion_tab.click()
        page.wait_for_timeout(1000)
        
        # Scope to discussion tab content
        discussion_content = page.get_by_label("💬 Discussion")
        
        # Fill thread title
        thread_title_input = discussion_content.get_by_placeholder("e.g., 'Question about methodology in cell 3'")
        thread_title_input.fill("Test thread for E2E")
        page.keyboard.press("Tab")  # Trigger Streamlit rerun
        page.wait_for_timeout(1000)
        
        # Fill optional cell number
        cell_input = discussion_content.get_by_label("Cell # (optional)")
        cell_input.fill("5")
        page.keyboard.press("Tab")
        page.wait_for_timeout(1000)
        
        # Click create thread button
        create_button = discussion_content.get_by_role("button", name="💬 Create Thread")
        expect(create_button).to_be_enabled()
        create_button.click()
        page.wait_for_timeout(2000)  # Wait for API call and rerun
        
        # Check for success message or API unavailable warning
        # Since API may not be running, either response is acceptable
        success_visible = discussion_content.get_by_text("✅ Thread created successfully!").is_visible()
        api_warning_visible = discussion_content.get_by_text("⚠️ API server unavailable").is_visible()
        
        assert success_visible or api_warning_visible, "Expected either success message or API unavailable warning"
    
    def test_add_comment_to_thread(self, page: Page):
        """Test adding a comment to an existing thread."""
        # Look for review expanders
        expanders = page.locator("[data-testid='stExpander']")
        
        if expanders.count() == 0:
            pytest.skip("No pending reviews to test comment addition")
        
        # Click first expander to open it
        expanders.first.click()
        page.wait_for_timeout(1000)
        
        # Click Discussion tab
        discussion_tab = page.get_by_role("tab", name="💬 Discussion")
        discussion_tab.click()
        page.wait_for_timeout(1000)
        
        # Scope to discussion tab content
        discussion_content = page.get_by_label("💬 Discussion")
        
        # Look for existing threads (demo data should be shown if API unavailable)
        thread_expanders = discussion_content.locator("[data-testid='stExpander']")
        
        if thread_expanders.count() > 0:
            # Click first thread expander
            thread_expanders.first.click()
            page.wait_for_timeout(1000)
            
            # Look for comment textarea
            comment_textareas = discussion_content.get_by_placeholder("Share your thoughts...")
            if comment_textareas.count() > 0:
                # Fill comment
                comment_textareas.first.fill("This is a test comment for E2E testing")
                page.keyboard.press("Tab")
                page.wait_for_timeout(1000)
                
                # Click Add Comment button
                add_comment_buttons = discussion_content.get_by_role("button", name="💬 Add Comment")
                if add_comment_buttons.count() > 0:
                    add_comment_buttons.first.click()
                    page.wait_for_timeout(2000)
                    
                    # Check for success or API warning
                    success_visible = discussion_content.get_by_text("✅ Comment added successfully!").is_visible()
                    api_warning_visible = discussion_content.get_by_text("⚠️ API server unavailable").is_visible()
                    
                    assert success_visible or api_warning_visible, "Expected either success message or API unavailable warning"
        else:
            # No threads visible, check for "No discussions yet" message
            expect(discussion_content.get_by_text("💭 No discussions yet")).to_be_visible()
    
    def test_diff_tab_visible(self, page: Page):
        """Test that Diff tab shows diff viewer components."""
        # Look for review expanders
        expanders = page.locator("[data-testid='stExpander']")
        
        if expanders.count() == 0:
            pytest.skip("No pending reviews to test Diff tab")
        
        # Click first expander to open it
        expanders.first.click()
        page.wait_for_timeout(1000)
        
        # Click Diff tab
        diff_tab = page.get_by_role("tab", name="📋 Diff")
        expect(diff_tab).to_be_visible()
        diff_tab.click()
        page.wait_for_timeout(1000)
        
        # Check for diff view components
        diff_content = page.get_by_label("📋 Diff")
        expect(diff_content.get_by_text("📋 Notebook Changes")).to_be_visible()
        
        # The diff content will either show metrics or "No changes detected" or API warning
        # Any of these is acceptable for the E2E test
        has_metrics = diff_content.get_by_text("Added").is_visible()
        has_no_changes = diff_content.get_by_text("✨ No changes detected").is_visible()
        has_api_warning = diff_content.get_by_text("⚠️ API server unavailable").is_visible()
        
        assert has_metrics or has_no_changes or has_api_warning, "Expected diff content, no changes message, or API warning"
    
    def test_diff_shows_summary(self, page: Page):
        """Test that diff viewer shows summary metrics when changes exist."""
        # Look for review expanders
        expanders = page.locator("[data-testid='stExpander']")
        
        if expanders.count() == 0:
            pytest.skip("No pending reviews to test Diff summary")
        
        # Click first expander to open it
        expanders.first.click()
        page.wait_for_timeout(1000)
        
        # Click Diff tab
        diff_tab = page.get_by_role("tab", name="📋 Diff")
        diff_tab.click()
        page.wait_for_timeout(1000)
        
        # Scope to diff tab content
        diff_content = page.get_by_label("📋 Diff")
        
        # Check for either demo data metrics, no changes message, or API warning
        # Demo data should show metrics like "Added", "Removed", "Modified"
        if diff_content.get_by_text("⚠️ API server unavailable").is_visible():
            # API unavailable, demo data should be shown
            expect(diff_content.get_by_text("Added")).to_be_visible()
            expect(diff_content.get_by_text("Removed")).to_be_visible()
            expect(diff_content.get_by_text("Modified")).to_be_visible()
            expect(diff_content.get_by_text("Unchanged")).to_be_visible()
        elif diff_content.get_by_text("✨ No changes detected").is_visible():
            # No changes case is also valid
            expect(diff_content.get_by_text("The notebook matches the reviewed snapshot")).to_be_visible()
        else:
            # Real API response with changes
            # Should show at least one of the metric labels
            metrics_visible = (
                diff_content.get_by_text("Added").is_visible() or
                diff_content.get_by_text("Removed").is_visible() or
                diff_content.get_by_text("Modified").is_visible() or
                diff_content.get_by_text("Unchanged").is_visible()
            )
            assert metrics_visible, "Expected to see at least one diff metric"
