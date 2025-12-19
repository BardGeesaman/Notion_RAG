"""Tests for Phase 2.7: Cross-Omics Helpers (Postgres-only)."""


from amprenta_rag.query.cross_omics.helpers import (
    extract_relation_ids,
    extract_select_values,
    extract_text_property,
    get_chunk_text,
    group_chunks_by_omics_type,
)


class TestModuleImport:
    """Test module imports successfully."""

    def test_helpers_module_imports(self):
        """Test helpers module imports without errors."""
        from amprenta_rag.query.cross_omics import helpers
        assert helpers is not None


class TestUtilityFunctions:
    """Test utility functions still work."""

    def test_extract_relation_ids(self):
        """Test extract_relation_ids with mock page dict."""
        mock_page = {
            "properties": {
                "Related Items": {
                    "type": "relation",
                    "relation": [
                        {"id": "id-1"},
                        {"id": "id-2"},
                    ]
                }
            }
        }
        result = extract_relation_ids(mock_page, "Related Items")
        assert result == ["id-1", "id-2"]

    def test_extract_relation_ids_empty(self):
        """Test extract_relation_ids with missing property."""
        mock_page = {"properties": {}}
        result = extract_relation_ids(mock_page, "NonExistent")
        assert result == []

    def test_extract_select_values_single(self):
        """Test extract_select_values with single select."""
        mock_page = {
            "properties": {
                "Status": {
                    "type": "select",
                    "select": {"name": "Active"}
                }
            }
        }
        result = extract_select_values(mock_page, "Status")
        assert result == ["Active"]

    def test_extract_select_values_multi(self):
        """Test extract_select_values with multi-select."""
        mock_page = {
            "properties": {
                "Tags": {
                    "type": "multi_select",
                    "multi_select": [
                        {"name": "tag1"},
                        {"name": "tag2"},
                    ]
                }
            }
        }
        result = extract_select_values(mock_page, "Tags")
        assert result == ["tag1", "tag2"]

    def test_extract_text_property_title(self):
        """Test extract_text_property with title."""
        mock_page = {
            "properties": {
                "Name": {
                    "type": "title",
                    "title": [{"plain_text": "Test Name"}]
                }
            }
        }
        result = extract_text_property(mock_page, "Name")
        assert result == "Test Name"

    def test_extract_text_property_rich_text(self):
        """Test extract_text_property with rich_text."""
        mock_page = {
            "properties": {
                "Description": {
                    "type": "rich_text",
                    "rich_text": [
                        {"plain_text": "Part 1"},
                        {"plain_text": " Part 2"},
                    ]
                }
            }
        }
        result = extract_text_property(mock_page, "Description")
        assert result == "Part 1 Part 2"

    def test_get_chunk_text_returns_snippet(self):
        """Test get_chunk_text returns snippet from metadata."""
        mock_chunk = {
            "metadata": {
                "snippet": "This is a test snippet."
            }
        }
        result = get_chunk_text(mock_chunk)
        assert result == "This is a test snippet."

    def test_get_chunk_text_no_snippet(self):
        """Test get_chunk_text returns None when no snippet."""
        mock_chunk = {"metadata": {}}
        result = get_chunk_text(mock_chunk)
        assert result is None

    def test_group_chunks_by_omics_type(self):
        """Test group_chunks_by_omics_type groups correctly."""
        mock_chunks = [
            {"metadata": {"omics_type": "Proteomics"}},
            {"metadata": {"omics_type": "Lipidomics"}},
            {"metadata": {"omics_type": "Proteomics"}},
            {"metadata": {"omics_type": "Metabolomics"}},
            {"metadata": {}},  # No omics_type -> Other
        ]
        result = group_chunks_by_omics_type(mock_chunks)
        assert "Proteomics" in result
        assert "Lipidomics" in result
        assert "Metabolomics" in result
        assert len(result["Proteomics"]) == 2
        assert len(result["Lipidomics"]) == 1
        assert len(result["Metabolomics"]) == 1
        assert len(result["Other"]) == 1

