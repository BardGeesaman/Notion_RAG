"""
Tests for metadata extraction helpers.
"""

from amprenta_rag.ingestion.metadata.helpers import (
    get_multi_names,
    get_relation_ids,
    get_select_name,
)


def test_get_select_name():
    """Test extracting select values."""
    props = {
        "Status": {
            "select": {"name": "Active"},
        },
    }
    result = get_select_name(props, "Status")
    assert result == "Active"

    # Missing property
    result = get_select_name(props, "Missing")
    assert result is None


def test_get_multi_names():
    """Test extracting multi-select values."""
    props = {
        "Disease": {
            "multi_select": [
                {"name": "ALS"},
                {"name": "AD"},
            ],
        },
    }
    result = get_multi_names(props, "Disease")
    assert len(result) == 2
    assert "ALS" in result
    assert "AD" in result

    # Empty multi-select
    props = {
        "Disease": {
            "multi_select": [],
        },
    }
    result = get_multi_names(props, "Disease")
    assert result == []


def test_get_relation_ids():
    """Test extracting relation IDs."""
    props = {
        "Related Datasets": {
            "relation": [
                {"id": "abc-123-def"},
                {"id": "xyz-789-ghi"},
            ],
        },
    }
    result = get_relation_ids(props, "Related Datasets")
    assert len(result) == 2
    # IDs should have dashes removed
    assert "abc123def" in result or "abc-123-def" in result

