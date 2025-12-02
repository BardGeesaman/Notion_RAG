# tests/test_sanitize_metadata.py

from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata


def test_sanitize_metadata_removes_none_and_converts_types():
    meta = {
        "ok_str": "value",
        "ok_int": 123,
        "none_val": None,
        "list_mixed": [1, "two", None, 3.5],
        "nested_dict": {"a": 1, "b": "x"},
    }

    cleaned = sanitize_metadata(meta)

    # None should be removed entirely
    assert "none_val" not in cleaned

    # Primitive values should pass through unchanged
    assert cleaned["ok_str"] == "value"
    assert cleaned["ok_int"] == 123

    # Lists should become list of strings (no None)
    assert cleaned["list_mixed"] == ["1", "two", "3.5"]

    # Dicts and other non-primitive values should be stringified
    assert isinstance(cleaned["nested_dict"], str)
    assert "a" in cleaned["nested_dict"] and "b" in cleaned["nested_dict"]