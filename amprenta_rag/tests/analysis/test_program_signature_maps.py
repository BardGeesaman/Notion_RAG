from __future__ import annotations

import sys
from types import ModuleType

from amprenta_rag.analysis import program_signature_maps as psm


class FakeSig:
    def __init__(self, sig_id: str, name: str):
        self.id = sig_id
        self.name = name


class FakeCompMatch:
    def __init__(self, feature_type: str, match_score: float, matched: bool = True):
        self.feature_type = feature_type
        self.match_score = match_score
        self.matched = matched


class FakeScoreResult:
    def __init__(self, total_score: float, component_matches):
        self.total_score = total_score
        self.component_matches = component_matches


def test_compute_program_signature_scores(monkeypatch) -> None:
    monkeypatch.setattr(psm, "_get_program_name", lambda pid: f"Program-{pid}")
    monkeypatch.setattr(psm, "get_program_datasets", lambda pid: ["d1", "d2"])
    monkeypatch.setattr(psm, "load_signatures_from_postgres", lambda: [FakeSig("s1", "Sig1"), FakeSig("s2", "Sig2")])
    monkeypatch.setattr(psm, "extract_dataset_features_by_type", lambda *args, **kwargs: {"gene": {"A", "B"}})

    def _score(signature, dataset_features_by_type):
        return FakeScoreResult(
            total_score=1.0,
            component_matches=[FakeCompMatch("gene", 0.5, True)],
        )

    fake_mos = ModuleType("amprenta_rag.ingestion.multi_omics_scoring")
    fake_mos.score_multi_omics_signature_against_dataset = _score
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.multi_omics_scoring", fake_mos)

    scores = psm.compute_program_signature_scores("prog1", use_cache=False)

    assert len(scores) == 2
    first = scores[0]
    assert first.program_name == "Program-prog1"
    assert first.coverage_fraction == 1.0
    assert first.score_by_omics["gene"] == 0.5
    assert set(first.matching_datasets) == {"d1", "d2"}


def test_compute_program_omics_coverage(monkeypatch) -> None:
    monkeypatch.setattr(psm, "_get_program_name", lambda pid: "Program-X")
    monkeypatch.setattr(psm, "get_program_datasets", lambda pid: ["d1", "d2"])

    def _extract(dataset_id, use_cache=True):
        if dataset_id == "d1":
            return {"gene": {"A"}, "protein": {"P1", "P2"}}
        return {"gene": {"B"}, "protein": set()}

    monkeypatch.setattr(psm, "extract_dataset_features_by_type", _extract)

    cov = psm.compute_program_omics_coverage("prog1", use_cache=False)

    assert cov.total_datasets == 2
    assert cov.datasets_by_omics["gene"] == 2
    assert cov.datasets_by_omics["protein"] == 1
    assert cov.features_by_omics["gene"] == 2
    assert "omics data" not in cov.coverage_summary.lower()


def test_convergence_indicators() -> None:
    scores = [
        psm.ProgramSignatureScore(
            program_id="p",
            signature_id="s1",
            program_name="P",
            signature_name="S1",
            overall_score=1.0,
            score_by_omics={"gene": 0.1, "protein": 0.2},
        ),
        psm.ProgramSignatureScore(
            program_id="p",
            signature_id="s2",
            program_name="P",
            signature_name="S2",
            overall_score=0.5,
            score_by_omics={"gene": 0.1},
        ),
    ]

    metrics = psm.compute_convergence_indicators(scores)

    assert metrics["multi_omics_signature_count"] == 1
    assert metrics["convergence_fraction"] == 0.5
    assert metrics["avg_omics_per_signature"] == 1.5


def test_generate_program_signature_map_and_report(monkeypatch) -> None:
    fake_scores = [
        psm.ProgramSignatureScore(
            program_id="p",
            signature_id="s1",
            program_name="P",
            signature_name="Sig1",
            overall_score=0.9,
            score_by_omics={"gene": 0.9},
            matching_datasets=["d1", "d2"],
            coverage_fraction=1.0,
        )
    ]
    fake_cov = psm.ProgramOmicsCoverage(
        program_id="p",
        program_name="Prog Name",
        total_datasets=2,
        datasets_by_omics={"gene": 2},
        features_by_omics={"gene": 3},
        coverage_summary="Gene: 2 dataset(s), 3 unique features",
    )

    monkeypatch.setattr(psm, "compute_program_signature_scores", lambda pid, use_cache=True: fake_scores)
    monkeypatch.setattr(psm, "compute_program_omics_coverage", lambda pid, use_cache=True: fake_cov)
    monkeypatch.setattr(psm, "compute_convergence_indicators", lambda scores: {"convergence_fraction": 1.0, "avg_omics_per_signature": 1.0, "multi_omics_signature_count": 1})

    prog_map = psm.generate_program_signature_map("p", top_n=1, use_cache=False)

    assert prog_map.top_signatures == fake_scores
    assert prog_map.convergence_indicators["convergence_fraction"] == 1.0

    report = psm.generate_program_map_report(prog_map, include_all_signatures=True)
    assert "Program-Signature Map" in report
    assert "Omics Coverage" in report
    assert "Sig1" in report


def test_get_program_datasets_and_name_fallbacks(monkeypatch) -> None:
    class _DummyQuery:
        def __init__(self, result=None):
            self.result = result

        def filter(self, *_args, **_kwargs):
            return self

        def first(self):
            return self.result

    class _DummyDB:
        def __init__(self, result=None):
            self.result = result

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def query(self, *_args, **_kwargs):
            return _DummyQuery(self.result)

    # No program found path
    monkeypatch.setattr(psm, "db_session", lambda: _DummyDB())
    assert psm.get_program_datasets("not-a-uuid") == []

    # Program found but no name
    class _Prog:
        def __init__(self):
            self.name = None
            self.datasets = []

    monkeypatch.setattr(psm, "db_session", lambda: _DummyDB(_Prog()))
    assert psm._get_program_name("prog1234") == "Program prog1234"

    # Exception path
    class _Boom:
        def __enter__(self):
            return self

        def query(self, *_args, **_kwargs):
            raise RuntimeError("boom")

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(psm, "db_session", lambda: _Boom())
    assert psm.get_program_datasets("prog") == []

