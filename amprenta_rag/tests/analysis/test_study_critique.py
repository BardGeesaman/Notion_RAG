from __future__ import annotations

import uuid

from amprenta_rag.analysis import study_critique as sc


class FakeExperiment:
    def __init__(
        self,
        sample_groups=None,
        design_type=None,
        organism=None,
        description="",
        datasets=None,
    ):
        self.sample_groups = sample_groups
        self.design_type = design_type
        self.organism = organism
        self.description = description
        self.datasets = datasets or []


class _Query:
    def __init__(self, obj=None):
        self.obj = obj

    def filter(self, *_args, **_kwargs):
        return self

    def first(self):
        return self.obj


class FakeDB:
    def __init__(self, experiment=None):
        self.experiment = experiment

    def query(self, *_args, **_kwargs):
        return _Query(self.experiment)


def test_assess_study_quality_not_found() -> None:
    db = FakeDB(None)
    result = sc.assess_study_quality(uuid.uuid4(), db)
    assert result["quality_score"] == 0
    assert "Experiment not found" in result["summary"]


def test_assess_study_quality_flaws_and_gaps() -> None:
    exp = FakeExperiment(
        sample_groups={"case": {"samples": ["s1"]}, "treated": {"samples": []}},
        design_type=None,
        organism=[],
        description="",
        datasets=[],
    )
    db = FakeDB(exp)

    flaws = sc.detect_design_flaws(exp, db)
    assert any("Missing control group" in f["flaw"] for f in flaws)
    assert any("Design type not specified" in f["flaw"] for f in flaws)

    gaps = sc.identify_data_gaps(exp, db)
    assert "Missing organism field" in gaps
    assert "No linked datasets" in gaps

    result = sc.assess_study_quality(uuid.uuid4(), db)
    assert result["quality_score"] < 100
    assert result["issues"]

