from __future__ import annotations

import uuid
from typing import Any, Dict, List, Set

from amprenta_rag.analysis import pattern_detection as pd


class _Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):  # type: ignore[override]
        return ("eq", self.name, other)

    def in_(self, items: List[Any]):
        return ("in", self.name, items)


class FakeDataset:
    id = _Field("dataset_id")

    def __init__(self, id_: str, features=None, signatures=None):
        self.id = uuid.UUID(id_)
        self.features = features or []
        self.signatures = signatures or []


class FakeSignature:
    id = _Field("signature_id")
    datasets: List[FakeDataset] = []

    def __init__(self, id_: str, datasets: List[FakeDataset]):
        self.id = id_
        self.datasets = datasets


class FakeFeature:
    id = _Field("feature_id")

    def __init__(self, id_: str, name: str):
        self.id = id_
        self.name = name


class FakeComponent:
    signature_id = _Field("signature_id")
    feature_id = _Field("feature_id")

    def __init__(self, signature_id: str, feature_name: str | None = None, feature_id: str | None = None):
        self.signature_id = signature_id
        self.feature_name = feature_name
        self.feature_id = feature_id


class FakeQuery:
    def __init__(self, cls, db):
        self.cls = cls
        self.db = db
        self.dataset_filter: Set[str] | None = None
        self.sig_filter: str | None = None
        self.feature_filter: str | None = None
        self.dataset_id: str | None = None

    def join(self, *_args, **_kwargs):
        return self

    def filter(self, condition):
        if isinstance(condition, tuple) and len(condition) == 3:
            kind, name, value = condition
            if kind == "in" and name == "dataset_id":
                self.dataset_filter = {str(v) for v in value}
            elif kind == "eq" and name == "signature_id":
                self.sig_filter = value
            elif kind == "eq" and name == "dataset_id":
                self.dataset_id = str(value)
            elif kind == "eq" and name == "feature_id":
                self.feature_filter = value
        return self

    def distinct(self):
        return self

    def all(self):
        if self.cls is pd.Signature:
            if self.dataset_filter is None:
                return self.db.signatures
            return [
                s for s in self.db.signatures if any(str(d.id) in self.dataset_filter for d in s.datasets)
            ]
        if self.cls is pd.SignatureComponent:
            return self.db.components_by_sig.get(self.sig_filter or "", [])
        return []

    def first(self):
        if self.cls is pd.Feature:
            return self.db.features_by_id.get(self.feature_filter or "")
        if self.cls is pd.Dataset:
            return self.db.datasets.get(self.dataset_id or "")
        return None


class FakeDB:
    def __init__(
        self,
        signatures: List[FakeSignature],
        components_by_sig: Dict[str, List[FakeComponent]],
        features_by_id: Dict[str, FakeFeature],
        datasets: Dict[str, FakeDataset],
    ):
        self.signatures = signatures
        self.components_by_sig = components_by_sig
        self.features_by_id = features_by_id
        self.datasets = datasets

    def query(self, cls):
        return FakeQuery(cls, self)


def _uuid(n: int) -> str:
    return f"00000000-0000-0000-0000-00000000000{n}"


def test_find_recurring_features_filters_by_min_occurrence(monkeypatch) -> None:
    d1 = FakeDataset(_uuid(1))
    d2 = FakeDataset(_uuid(2))
    sig = FakeSignature("sig1", datasets=[d1, d2])
    comp1 = FakeComponent("sig1", feature_name="GeneA")
    comp2 = FakeComponent("sig1", feature_name="GeneB")

    db = FakeDB(
        signatures=[sig],
        components_by_sig={"sig1": [comp1, comp2]},
        features_by_id={},
        datasets={_uuid(1): d1, _uuid(2): d2},
    )

    monkeypatch.setattr(pd, "Signature", FakeSignature)
    monkeypatch.setattr(pd, "SignatureComponent", FakeComponent)
    monkeypatch.setattr(pd, "Dataset", FakeDataset)
    monkeypatch.setattr(pd, "Feature", FakeFeature)

    results = pd.find_recurring_features([_uuid(1), _uuid(2)], db, min_occurrence=2)

    assert len(results) == 2
    assert results[0]["feature_name"] in {"GeneA", "GeneB"}
    assert results[0]["occurrence_count"] == 2


def test_calculate_overlap_collects_features_and_components(monkeypatch) -> None:
    d1 = FakeDataset(_uuid(1))
    d2 = FakeDataset(_uuid(2))

    sig1 = FakeSignature("sig1", datasets=[d1])
    sig2 = FakeSignature("sig2", datasets=[d2])

    featA = FakeFeature("fA", "A")
    featB = FakeFeature("fB", "B")

    d1.features = [featA]
    d1.signatures = [sig1]
    d2.features = [featB]
    d2.signatures = [sig2]

    comp1 = FakeComponent("sig1", feature_name=None, feature_id="fB")
    comp2 = FakeComponent("sig2", feature_name="C")

    db = FakeDB(
        signatures=[sig1, sig2],
        components_by_sig={"sig1": [comp1], "sig2": [comp2]},
        features_by_id={"fB": featB},
        datasets={_uuid(1): d1, _uuid(2): d2},
    )

    monkeypatch.setattr(pd, "Signature", FakeSignature)
    monkeypatch.setattr(pd, "SignatureComponent", FakeComponent)
    monkeypatch.setattr(pd, "Dataset", FakeDataset)
    monkeypatch.setattr(pd, "Feature", FakeFeature)

    overlap = pd.calculate_overlap([_uuid(1), _uuid(2)], db)

    assert overlap[_uuid(1)] == {"A", "B"}
    assert overlap[_uuid(2)] == {"B", "C"}


def test_get_cross_dataset_summary_uses_helpers(monkeypatch) -> None:
    monkeypatch.setattr(
        pd,
        "calculate_overlap",
        lambda dataset_ids, db: {
            "d1": {"A", "B"},
            "d2": {"B", "C"},
        },
    )
    monkeypatch.setattr(
        pd,
        "find_recurring_features",
        lambda dataset_ids, db, min_occurrence=2: [
            {"feature_name": "B", "occurrence_count": 2, "datasets": ["d1", "d2"]}
        ],
    )

    summary = pd.get_cross_dataset_summary(["d1", "d2"], db=None)

    assert summary["total_features_per_dataset"] == {"d1": 2, "d2": 2}
    assert summary["overlap_counts"]  # should have pairwise overlap
    assert summary["top_recurring_features"][0]["feature_name"] == "B"
    assert summary["unique_features_per_dataset"]["d1"] == 1
    assert summary["unique_features_per_dataset"]["d2"] == 1


