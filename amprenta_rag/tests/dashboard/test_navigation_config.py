"""Tests for dashboard navigation configuration."""

from __future__ import annotations

import pytest

from scripts.dashboard.core.config import (
    PAGE_REGISTRY,
    PAGE_GROUPS,
    ALL_PAGES,
    ADMIN_PAGES,
)


class TestNavigationConfig:
    """Test navigation configuration consistency."""
    
    def test_no_duplicate_pages(self):
        """Test that no page appears in multiple groups."""
        all_grouped_pages = []
        
        for group_name, pages in PAGE_GROUPS.items():
            all_grouped_pages.extend(pages)
        
        # Check for duplicates
        seen = set()
        duplicates = []
        
        for page in all_grouped_pages:
            if page in seen:
                duplicates.append(page)
            else:
                seen.add(page)
        
        assert not duplicates, f"Duplicate pages found in PAGE_GROUPS: {duplicates}"
    
    def test_all_registry_pages_grouped(self):
        """Test that every PAGE_REGISTRY key is in some PAGE_GROUP."""
        all_grouped_pages = []
        
        for group_name, pages in PAGE_GROUPS.items():
            all_grouped_pages.extend(pages)
        
        grouped_pages_set = set(all_grouped_pages)
        registry_pages_set = set(PAGE_REGISTRY.keys())
        
        ungrouped_pages = registry_pages_set - grouped_pages_set
        
        assert not ungrouped_pages, f"Pages in PAGE_REGISTRY but not in any PAGE_GROUP: {sorted(ungrouped_pages)}"
        
        # Also check for pages in groups but not in registry
        unregistered_pages = grouped_pages_set - registry_pages_set
        assert not unregistered_pages, f"Pages in PAGE_GROUPS but not in PAGE_REGISTRY: {sorted(unregistered_pages)}"
    
    def test_admin_pages_derived(self):
        """Test that ADMIN_PAGES matches PAGE_GROUPS['Admin']."""
        expected_admin_pages = PAGE_GROUPS.get("Admin", [])
        
        assert ADMIN_PAGES == expected_admin_pages, (
            f"ADMIN_PAGES mismatch. Expected: {expected_admin_pages}, "
            f"Got: {ADMIN_PAGES}"
        )
    
    def test_all_pages_derived(self):
        """Test that ALL_PAGES matches PAGE_REGISTRY.keys()."""
        expected_all_pages = list(PAGE_REGISTRY.keys())
        
        # Sort both for comparison since order doesn't matter
        assert sorted(ALL_PAGES) == sorted(expected_all_pages), (
            f"ALL_PAGES mismatch. Expected: {sorted(expected_all_pages)}, "
            f"Got: {sorted(ALL_PAGES)}"
        )
    
    def test_page_groups_structure(self):
        """Test that PAGE_GROUPS has expected structure."""
        assert isinstance(PAGE_GROUPS, dict), "PAGE_GROUPS must be a dictionary"
        
        for group_name, pages in PAGE_GROUPS.items():
            assert isinstance(group_name, str), f"Group name must be string: {group_name}"
            assert isinstance(pages, list), f"Group pages must be list: {group_name}"
            assert all(isinstance(page, str) for page in pages), f"All pages must be strings in group: {group_name}"
    
    def test_page_registry_structure(self):
        """Test that PAGE_REGISTRY has expected structure."""
        assert isinstance(PAGE_REGISTRY, dict), "PAGE_REGISTRY must be a dictionary"
        
        for page_name, (module_path, function_name) in PAGE_REGISTRY.items():
            assert isinstance(page_name, str), f"Page name must be string: {page_name}"
            assert isinstance(module_path, str), f"Module path must be string for page: {page_name}"
            assert isinstance(function_name, str), f"Function name must be string for page: {page_name}"
            
            # Check module path format
            assert module_path.startswith("scripts.dashboard.pages"), (
                f"Module path should start with 'scripts.dashboard.pages' for page: {page_name}"
            )
    
    def test_essential_groups_exist(self):
        """Test that essential page groups exist."""
        essential_groups = ["Home", "Admin", "Discovery", "Chemistry"]
        
        for group in essential_groups:
            assert group in PAGE_GROUPS, f"Essential group '{group}' missing from PAGE_GROUPS"
            assert len(PAGE_GROUPS[group]) > 0, f"Essential group '{group}' is empty"
    
    def test_essential_pages_exist(self):
        """Test that essential pages exist in registry."""
        essential_pages = ["Overview", "Chemistry", "Company Settings"]
        
        for page in essential_pages:
            assert page in PAGE_REGISTRY, f"Essential page '{page}' missing from PAGE_REGISTRY"