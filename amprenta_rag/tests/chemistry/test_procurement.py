from __future__ import annotations

from amprenta_rag.chemistry import procurement


def test_get_vendor_info_contains_ids():
    info = procurement.get_vendor_info()
    assert isinstance(info, list)
    ids = {v["id"] for v in info}
    assert "sigma_aldrich" in ids
    assert "fisher" in ids


def test_search_vendors_returns_mock_results():
    results = procurement.search_vendors("aspirin")
    assert len(results) >= 3
    first = results[0]
    assert {"vendor", "catalog_id", "price", "availability", "url"}.issubset(first.keys())

