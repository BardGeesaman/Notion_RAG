from __future__ import annotations

from amprenta_rag.ingestion import design_extraction as de


def test_detect_design_type_patterns():
    design, score = de.detect_design_type(
        sample_names=["case1", "control1"],
        sample_attributes={"attr": ["treated", "placebo"]},
        study_description="after treatment week1",
    )
    assert design in {"case_control", "intervention", "time_course", "dose_response", "multi_factorial"}
    assert 0.5 <= score <= 1.0


def test_extract_sample_groups_case_control():
    groups = de.extract_sample_groups(
        ["CaseA", "ControlB", "UnknownC"],
        {"status": ["case", "control", "idk"]},
        design_type="case_control",
    )
    assert "control" in groups and "case" in groups and "unknown" in groups


def test_extract_sample_groups_time_course():
    groups = de.extract_sample_groups(
        ["T0_sample", "day3_sample", "other"],
        {},
        design_type="time_course",
    )
    assert "timepoints" in groups and any(groups["timepoints"].values())


def test_extract_geo_design_basic():
    geo_meta = {
        "Samples": [
            {"Title": "Sample1", "Characteristics": ["status: case"]},
            {"Title": "Sample2", "Characteristics": ["status: control"]},
        ],
        "Summary": "case control study",
    }
    res = de.extract_geo_design(geo_meta)
    assert res["design_type"]
    assert res["sample_groups"]


def test_extract_mw_design_uses_factors():
    mw_meta = {
        "study_id": "S1",
        "factors": [
            {"Factor": "time", "FactorValues": ["T0", "T1"]},
        ],
    }
    res = de.extract_mw_design(mw_meta)
    assert res["design_type"]
    # sample_groups may be empty depending on factors; ensure no exception and structure present
    assert "sample_groups" in res

