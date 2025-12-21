from __future__ import annotations

from amprenta_rag.ingestion import design_integration as di


def test_merge_designs_prefers_higher_confidence():
    existing = {"design_type": "observational", "confidence": 0.5, "sample_groups": {"all": ["s1"]}}
    new = {"design_type": "case_control", "confidence": 0.8, "sample_groups": {"case": ["s1"], "control": ["s2"]}}
    merged = di.merge_designs(existing, new)
    assert merged["design_type"] == "case_control"
    assert merged["sample_groups"]["case"] == ["s1"]


def test_merge_designs_fills_missing_metadata():
    existing = {"design_type": "time_course", "confidence": 0.6, "sample_groups": {}}
    new = {"design_type": "time_course", "confidence": 0.5, "design_metadata": {"extra": True}}
    merged = di.merge_designs(existing, new)
    assert merged["design_metadata"]["extra"] is True

