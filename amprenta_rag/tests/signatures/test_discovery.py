import json
from pathlib import Path
from typing import List

import pytest

from amprenta_rag.signatures.discovery import (
    FeatureOccurrence,
    SignatureCandidate,
    compute_direction_consistency,
    compute_feature_co_occurrence,
    discover_signature_candidates,
    export_candidates_to_json,
    export_candidates_to_tsv,
)


def _make_occurrence(
    name: str,
    ftype: str = "gene",
    direction: str | None = "up",
    dataset_id: str = "d1",
) -> FeatureOccurrence:
    return FeatureOccurrence(
        feature_name=name,
        feature_type=ftype,
        direction=direction,
        dataset_id=dataset_id,
    )


def test_compute_feature_co_occurrence_filters_minimum() -> None:
    all_features = {
        "d1": [_make_occurrence("A"), _make_occurrence("B"), _make_occurrence("C")],
        "d2": [_make_occurrence("A"), _make_occurrence("B")],
        "d3": [_make_occurrence("A"), _make_occurrence("C")],
    }

    counts = compute_feature_co_occurrence(all_features, min_co_occurrence=2)

    # (A, B) appears in d1 and d2; (A, C) appears in d1 and d3; (B, C) only once
    assert counts == {("A", "B"): 2, ("A", "C"): 2}


def test_compute_direction_consistency_majority_vote() -> None:
    all_features = {
        "d1": [_make_occurrence("A", direction="up")],
        "d2": [_make_occurrence("A", direction="down")],
        "d3": [_make_occurrence("A", direction="up")],
    }

    score = compute_direction_consistency("A", ["d1", "d2", "d3"], all_features)

    # Two of three datasets agree on "up"
    assert score == pytest.approx(2 / 3)


def test_discover_signature_candidates_happy_path(monkeypatch: pytest.MonkeyPatch) -> None:
    dataset_ids = ["d1", "d2", "d3"]

    all_features = {
        "d1": [
            _make_occurrence("A", "gene", "up", "d1"),
            _make_occurrence("B", "gene", "up", "d1"),
            _make_occurrence("C", "protein", "down", "d1"),
        ],
        "d2": [
            _make_occurrence("A", "gene", "up", "d2"),
            _make_occurrence("B", "gene", "up", "d2"),
            _make_occurrence("C", "protein", "down", "d2"),
        ],
        "d3": [
            _make_occurrence("A", "gene", "up", "d3"),
            _make_occurrence("B", "gene", "up", "d3"),
            _make_occurrence("C", "protein", "down", "d3"),
        ],
    }

    monkeypatch.setattr(
        "amprenta_rag.signatures.discovery.extract_features_from_datasets",
        lambda *_args, **_kwargs: all_features,
    )

    candidates = discover_signature_candidates(
        dataset_ids,
        min_support=3,
        min_features=3,
        max_features=5,
        min_co_occurrence=2,
        min_confidence=0.1,
    )

    assert len(candidates) == 1
    cand = candidates[0]
    assert cand.support_count == 3
    assert set(cand.datasets) == set(dataset_ids)
    assert cand.co_occurrence_score == pytest.approx(1.0)
    assert {"gene", "protein"} == cand.feature_types
    # All directions are consistent within each feature, so confidence should be high
    assert cand.confidence >= 0.5
    assert len(cand.features) == 3
    assert any("Discovered" in cand.name for cand in candidates)


def test_discover_signature_candidates_no_features(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "amprenta_rag.signatures.discovery.extract_features_from_datasets",
        lambda *_args, **_kwargs: {},
    )

    candidates = discover_signature_candidates(["d1", "d2"], min_support=1)

    assert candidates == []


def test_export_candidates_to_tsv(tmp_path: Path) -> None:
    candidates: List[SignatureCandidate] = [
        SignatureCandidate(
            name="Discovered-gene-1",
            features=[_make_occurrence("A", "gene", "up")],
            feature_types={"gene"},
            datasets=["d1", "d2"],
            co_occurrence_score=0.9,
            direction_consistency=1.0,
            support_count=2,
            confidence=0.8,
        )
    ]

    out_file = tmp_path / "candidates.tsv"
    export_candidates_to_tsv(candidates, str(out_file))

    content = out_file.read_text().splitlines()
    assert content[0].startswith("signature_name")
    assert "Discovered-gene-1" in content[1]
    assert "A" in content[1]
    assert "0.800" in content[1]


def test_export_candidates_to_json(tmp_path: Path) -> None:
    candidates = [
        SignatureCandidate(
            name="Discovered-gene-1",
            features=[_make_occurrence("A", "gene", "up")],
            feature_types={"gene"},
            datasets=["d1"],
            co_occurrence_score=1.0,
            direction_consistency=1.0,
            support_count=1,
            confidence=0.7,
        )
    ]

    out_file = tmp_path / "candidates.json"
    export_candidates_to_json(candidates, str(out_file))

    data = json.loads(out_file.read_text())
    assert data["summary"]["total_candidates"] == 1
    assert data["candidates"][0]["name"] == "Discovered-gene-1"
    assert data["candidates"][0]["features"][0]["feature_name"] == "A"


