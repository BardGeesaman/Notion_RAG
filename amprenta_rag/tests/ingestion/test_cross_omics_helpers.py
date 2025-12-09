"""
Tests for cross_omics helpers module.
"""


from amprenta_rag.query.cross_omics.helpers import (
    extract_relation_ids,
    extract_select_values,
    extract_text_property,
    group_chunks_by_omics_type,
)


def test_extract_relation_ids():
    """Test extracting relation IDs from Notion properties."""
    # Function expects a full page dict with properties
    page = {
        "properties": {
            "Related Datasets": {
                "type": "relation",
                "relation": [
                    {"id": "abc-123-def"},
                    {"id": "xyz-789-ghi"},
                ],
            },
        },
    }
    result = extract_relation_ids(page, "Related Datasets")
    assert len(result) == 2
    # Function returns IDs with dashes (not removed)
    assert "abc-123-def" in result
    assert "xyz-789-ghi" in result


def test_extract_select_values():
    """Test extracting select values from Notion properties."""
    # Test single select - function expects full page dict
    page = {
        "properties": {
            "Status": {
                "type": "select",
                "select": {"name": "Active"},
            },
        },
    }
    result = extract_select_values(page, "Status")
    assert result == ["Active"]

    # Test multi-select
    page = {
        "properties": {
            "Modalities": {
                "type": "multi_select",
                "multi_select": [
                    {"name": "Gene"},
                    {"name": "Protein"},
                ],
            },
        },
    }
    result = extract_select_values(page, "Modalities")
    assert len(result) == 2
    assert "Gene" in result
    assert "Protein" in result


def test_extract_text_property():
    """Test extracting text from title or rich_text properties."""
    # Test title - function expects full page dict
    page = {
        "properties": {
            "Name": {
                "type": "title",
                "title": [{"plain_text": "Test Dataset"}],
            },
        },
    }
    result = extract_text_property(page, "Name")
    assert result == "Test Dataset"

    # Test rich_text
    page = {
        "properties": {
            "Description": {
                "type": "rich_text",
                "rich_text": [
                    {"plain_text": "First part"},
                    {"plain_text": " second part"},
                ],
            },
        },
    }
    result = extract_text_property(page, "Description")
    assert result == "First part second part"


def test_group_chunks_by_omics_type():
    """Test grouping chunks by omics type."""
    chunks = [
        {"metadata": {"omics_type": "Lipidomics"}},
        {"metadata": {"omics_type": "Metabolomics"}},
        {"metadata": {"omics_type": "Proteomics"}},
        {"metadata": {"omics_type": "Transcriptomics"}},
        {"metadata": {"omics_type": "Unknown"}},
        {"metadata": {}},  # No omics_type
    ]
    result = group_chunks_by_omics_type(chunks)
    assert len(result["Lipidomics"]) == 1
    assert len(result["Metabolomics"]) == 1
    assert len(result["Proteomics"]) == 1
    assert len(result["Transcriptomics"]) == 1
    assert len(result["Other"]) == 2  # Unknown + no type

