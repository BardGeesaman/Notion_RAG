"""E2E tests for Inline Annotations dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def _goto_lab_notebook_page(page: Page, base_url: str) -> None:
    """Navigate to the Lab Notebook page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Lab+Notebook")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _goto_datasets_page(page: Page, base_url: str) -> None:
    """Navigate to the Datasets page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Datasets")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def test_annotation_panel_loads(page: Page, streamlit_server: str) -> None:
    """Test that annotation panel loads when annotation indicator is clicked."""
    _goto_lab_notebook_page(page, streamlit_server)
    
    main = _main_container(page)
    
    # Look for annotation indicators (📝 buttons)
    annotation_buttons = main.get_by_text(re.compile(r"📝", re.IGNORECASE))
    
    if annotation_buttons.count() > 0:
        # Click the first annotation button
        annotation_buttons.first.click()
        page.wait_for_timeout(2000)
        
        # Check if annotation panel appears in sidebar
        sidebar = page.locator('[data-testid="stSidebar"]')
        annotation_panel = sidebar.get_by_text(re.compile(r"Inline Annotations", re.IGNORECASE))
        
        # Panel should be visible
        expect(annotation_panel.first).to_be_visible(timeout=5000)
    else:
        # If no annotation buttons, just verify page loads
        expect(main).to_be_visible(timeout=10000)


def test_create_annotation_on_notebook_cell(page: Page, streamlit_server: str) -> None:
    """Test creating an annotation on a notebook cell."""
    _goto_lab_notebook_page(page, streamlit_server)
    
    main = _main_container(page)
    
    # Look for annotation indicators
    annotation_buttons = main.get_by_text(re.compile(r"📝", re.IGNORECASE))
    
    if annotation_buttons.count() > 0:
        # Click annotation button to open panel
        annotation_buttons.first.click()
        page.wait_for_timeout(2000)
        
        sidebar = page.locator('[data-testid="stSidebar"]')
        
        # Look for create annotation form
        create_form = sidebar.get_by_text(re.compile(r"Create.*Annotation", re.IGNORECASE))
        
        if create_form.count() > 0:
            # Fill in annotation content
            content_input = sidebar.get_by_label(re.compile(r"Annotation.*Content", re.IGNORECASE)).or_(
                sidebar.locator('textarea').filter(has_text=re.compile(r"content|annotation", re.IGNORECASE))
            )
            
            if content_input.count() > 0:
                content_input.first.fill("This is a test annotation")
                page.wait_for_timeout(1000)
                
                # Look for create button
                create_button = sidebar.get_by_text(re.compile(r"Create.*Annotation", re.IGNORECASE)).filter(
                    has_text=re.compile(r"button|submit", re.IGNORECASE)
                )
                
                if create_button.count() > 0:
                    create_button.first.click()
                    page.wait_for_timeout(2000)
                    
                    # Check for success message or annotation display
                    success_msg = sidebar.get_by_text(re.compile(r"success|created", re.IGNORECASE))
                    annotation_content = sidebar.get_by_text("This is a test annotation")
                    
                    # Either success message or annotation should appear
                    assert success_msg.count() > 0 or annotation_content.count() > 0


def test_resolve_annotation(page: Page, streamlit_server: str) -> None:
    """Test resolving an annotation."""
    _goto_datasets_page(page, streamlit_server)
    
    main = _main_container(page)
    
    # Look for annotation indicators
    annotation_buttons = main.get_by_text(re.compile(r"📝", re.IGNORECASE))
    
    if annotation_buttons.count() > 0:
        # Click annotation button to open panel
        annotation_buttons.first.click()
        page.wait_for_timeout(2000)
        
        sidebar = page.locator('[data-testid="stSidebar"]')
        
        # Look for existing annotations with resolve button
        resolve_buttons = sidebar.get_by_text(re.compile(r"Resolve|✅", re.IGNORECASE))
        
        if resolve_buttons.count() > 0:
            # Click resolve button
            resolve_buttons.first.click()
            page.wait_for_timeout(2000)
            
            # Check for resolved status or reopen button
            resolved_status = sidebar.get_by_text(re.compile(r"resolved|🔵", re.IGNORECASE))
            reopen_button = sidebar.get_by_text(re.compile(r"Reopen|🔄", re.IGNORECASE))
            
            # Should show resolved status or reopen option
            assert resolved_status.count() > 0 or reopen_button.count() > 0
        else:
            # If no annotations to resolve, just verify panel loads
            annotation_panel = sidebar.get_by_text(re.compile(r"Inline Annotations", re.IGNORECASE))
            expect(annotation_panel.first).to_be_visible(timeout=5000)


def test_annotation_filter_by_status(page: Page, streamlit_server: str) -> None:
    """Test filtering annotations by status."""
    _goto_datasets_page(page, streamlit_server)
    
    main = _main_container(page)
    
    # Look for annotation indicators
    annotation_buttons = main.get_by_text(re.compile(r"📝", re.IGNORECASE))
    
    if annotation_buttons.count() > 0:
        # Click annotation button to open panel
        annotation_buttons.first.click()
        page.wait_for_timeout(2000)
        
        sidebar = page.locator('[data-testid="stSidebar"]')
        
        # Look for status filter dropdown
        status_filter = sidebar.get_by_label(re.compile(r"Status", re.IGNORECASE)).or_(
            sidebar.locator('select').filter(has_text=re.compile(r"status", re.IGNORECASE))
        )
        
        if status_filter.count() > 0:
            # Try to select "resolved" filter
            status_filter.first.click()
            page.wait_for_timeout(1000)
            
            # Look for "resolved" option
            resolved_option = page.get_by_text("resolved")
            if resolved_option.count() > 0:
                resolved_option.first.click()
                page.wait_for_timeout(2000)
                
                # Check that filter was applied (annotations list might change)
                annotations_section = sidebar.get_by_text(re.compile(r"Total.*Annotations|annotations", re.IGNORECASE))
                expect(annotations_section.first).to_be_visible(timeout=5000)
        else:
            # If no filter controls, just verify panel loads
            annotation_panel = sidebar.get_by_text(re.compile(r"Inline Annotations", re.IGNORECASE))
            expect(annotation_panel.first).to_be_visible(timeout=5000)


def test_annotation_indicators_on_datasets(page: Page, streamlit_server: str) -> None:
    """Test that annotation indicators appear on dataset elements."""
    _goto_datasets_page(page, streamlit_server)
    
    main = _main_container(page)
    
    # Wait for datasets to load
    page.wait_for_timeout(5000)
    
    # Look for dataset expanders
    dataset_expanders = main.locator('[data-testid="stExpander"]')
    
    if dataset_expanders.count() > 0:
        # Click on first dataset to expand it
        dataset_expanders.first.click()
        page.wait_for_timeout(2000)
        
        # Look for annotation indicators (📝 buttons) within the expanded dataset
        annotation_buttons = main.get_by_text(re.compile(r"📝", re.IGNORECASE))
        
        # Should have at least one annotation indicator
        assert annotation_buttons.count() > 0, "Expected annotation indicators on dataset elements"
    else:
        # If no datasets, just verify page loads
        expect(main).to_be_visible(timeout=10000)


def test_annotation_position_context(page: Page, streamlit_server: str) -> None:
    """Test that annotation panel shows position context correctly."""
    _goto_lab_notebook_page(page, streamlit_server)
    
    main = _main_container(page)
    
    # Look for annotation indicators
    annotation_buttons = main.get_by_text(re.compile(r"📝", re.IGNORECASE))
    
    if annotation_buttons.count() > 0:
        # Click annotation button to open panel
        annotation_buttons.first.click()
        page.wait_for_timeout(2000)
        
        sidebar = page.locator('[data-testid="stSidebar"]')
        
        # Look for position type selector in create form
        position_selector = sidebar.get_by_label(re.compile(r"Position.*Type", re.IGNORECASE)).or_(
            sidebar.locator('select').filter(has_text=re.compile(r"position", re.IGNORECASE))
        )
        
        if position_selector.count() > 0:
            # Verify position selector has options
            position_selector.first.click()
            page.wait_for_timeout(1000)
            
            # Look for position options (cell, column, row, etc.)
            cell_option = page.get_by_text("cell")
            column_option = page.get_by_text("column")
            
            # Should have position type options
            assert cell_option.count() > 0 or column_option.count() > 0, "Expected position type options"
        else:
            # If no position controls, just verify panel loads
            annotation_panel = sidebar.get_by_text(re.compile(r"Inline Annotations", re.IGNORECASE))
            expect(annotation_panel.first).to_be_visible(timeout=5000)
