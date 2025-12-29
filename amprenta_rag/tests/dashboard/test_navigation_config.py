"""Tests for navigation configuration."""

from __future__ import annotations

import pytest

from scripts.dashboard.core.config import (
    PAGE_GROUPS,
    GROUP_ORDER,
    GROUP_ICONS,
    ALL_PAGES,
)


class TestPageGroups:
    """Tests for PAGE_GROUPS structure."""

    def test_all_groups_in_order(self):
        """Test that all groups are in GROUP_ORDER."""
        for group in PAGE_GROUPS.keys():
            assert group in GROUP_ORDER, f"Group {group} not in GROUP_ORDER"

    def test_all_groups_have_icons(self):
        """Test that all groups have icons."""
        for group in PAGE_GROUPS.keys():
            assert group in GROUP_ICONS, f"Group {group} missing icon"

    def test_no_duplicate_pages_across_groups(self):
        """Test that pages don't appear in multiple groups."""
        all_grouped_pages = []
        for pages in PAGE_GROUPS.values():
            all_grouped_pages.extend(pages)
        
        # Check for duplicates
        assert len(all_grouped_pages) == len(set(all_grouped_pages)), "Pages appear in multiple groups"

    def test_groups_contain_valid_pages(self):
        """Test that grouped pages exist in ALL_PAGES."""
        for group, pages in PAGE_GROUPS.items():
            for page in pages:
                assert page in ALL_PAGES, f"Page '{page}' in group '{group}' not in ALL_PAGES"


class TestGroupOrder:
    """Tests for GROUP_ORDER."""

    def test_group_order_length(self):
        """Test GROUP_ORDER has all groups."""
        assert len(GROUP_ORDER) == len(PAGE_GROUPS)

    def test_home_first(self):
        """Test Home group is first."""
        assert GROUP_ORDER[0] == "Home"

    def test_admin_near_end(self):
        """Test Admin is near end."""
        assert "Admin" in GROUP_ORDER[-3:], "Admin should be near end"


class TestGroupIcons:
    """Tests for GROUP_ICONS."""

    def test_all_icons_are_emoji(self):
        """Test all icons are emoji characters."""
        for icon in GROUP_ICONS.values():
            assert len(icon) >= 1, "Icons should be non-empty"
            # Basic check: emoji are typically 1-2 chars
            assert len(icon) <= 3, "Icons should be short (emoji)"

