"""Tests for Phase 1.1b: Package __init__.py Notion import removals."""

import inspect
from amprenta_rag import chemistry, migration


class TestChemistryPackageImports:
    """Test chemistry package imports after Notion removal."""

    def test_package_imports_successfully(self):
        """Test that chemistry package imports without errors."""
        # This import is implicit by importing from amprenta_rag above,
        # but explicit import ensures we are testing what we think we are.
        assert True  # If we get here, import succeeded

    def test_removed_notion_functions_not_available(self):
        """Test that removed Notion functions are not available."""
        removed_functions = [
            'create_compound_feature_page',
            'find_or_create_compound_page',
            'create_hts_campaign_page',
            'find_or_create_hts_campaign_page',
            'create_biochemical_hit_page',
        ]

        for func_name in removed_functions:
            assert not hasattr(chemistry, func_name), \
                f"{func_name} should be removed from chemistry package"

    def test_no_notion_imports_in_file(self):
        """Verify chemistry/__init__.py doesn't contain Notion imports."""
        chem_source = inspect.getsource(chemistry)
        assert 'from amprenta_rag.chemistry.notion_integration' not in chem_source, \
            "Should not import notion_integration"


class TestMigrationPackageImports:
    """Test migration package imports after Notion removal."""

    def test_package_imports_successfully(self):
        """Test that migration package imports without errors."""
        assert True  # If we get here, import succeeded

    def test_no_notion_imports_in_file(self):
        """Verify migration/__init__.py doesn't contain Notion imports."""
        mig_source = inspect.getsource(migration)
        assert 'from amprenta_rag.migration.dual_write' not in mig_source, \
            "Should not import dual_write"

