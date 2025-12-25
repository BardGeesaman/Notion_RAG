from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4

import numpy as np
import pytest


pd = pytest.importorskip("pandas")


def test_load_biomarker_dataset_filters_low_coverage(monkeypatch):
    from amprenta_rag.ml.biomarker import datasets as dsm

    exp_id = uuid4()
    ds_id = uuid4()
    fake_exp = SimpleNamespace(id=exp_id, datasets=[SimpleNamespace(id=ds_id, omics_type="transcriptomics")])

    class FakeQuery:
        def filter(self, *args, **kwargs):  # noqa: ANN001
            return self

        def first(self):
            return fake_exp

    class FakeDB:
        def query(self, *args, **kwargs):  # noqa: ANN001
            return FakeQuery()

    @contextmanager
    def fake_db_session():
        yield FakeDB()

    monkeypatch.setattr(dsm, "db_session", fake_db_session)

    # features x samples with low coverage feature "bad"
    mat = pd.DataFrame(
        {
            "S1": [1.0, np.nan, 3.0],
            "S2": [2.0, np.nan, np.nan],
            "S3": [3.0, 5.0, 6.0],
            "S4": [4.0, np.nan, 7.0],
        },
        index=["good1", "bad", "good2"],
    )
    monkeypatch.setattr(dsm, "extract_feature_matrix", lambda dataset_id, db: mat)

    X, y, features = dsm.load_biomarker_dataset(
        experiment_id=exp_id,
        group1_samples=["S1", "S2"],
        group2_samples=["S3", "S4"],
        omics_type="transcriptomics",
        min_coverage=0.6,
    )

    assert "bad" not in features
    assert X.shape[0] == 4
    assert X.shape[1] == 2  # only good1, good2
    assert not np.isnan(X).any()


def test_load_biomarker_dataset_returns_binary_labels(monkeypatch):
    from amprenta_rag.ml.biomarker import datasets as dsm

    exp_id = uuid4()
    ds_id = uuid4()
    fake_exp = SimpleNamespace(id=exp_id, datasets=[SimpleNamespace(id=ds_id, omics_type="transcriptomics")])

    class FakeQuery:
        def filter(self, *args, **kwargs):  # noqa: ANN001
            return self

        def first(self):
            return fake_exp

    class FakeDB:
        def query(self, *args, **kwargs):  # noqa: ANN001
            return FakeQuery()

    @contextmanager
    def fake_db_session():
        yield FakeDB()

    monkeypatch.setattr(dsm, "db_session", fake_db_session)

    mat = pd.DataFrame({"S1": [1.0], "S2": [2.0], "S3": [3.0]}, index=["f0"])
    monkeypatch.setattr(dsm, "extract_feature_matrix", lambda dataset_id, db: mat)

    X, y, features = dsm.load_biomarker_dataset(
        experiment_id=exp_id,
        group1_samples=["S1", "S2"],
        group2_samples=["S3"],
        omics_type="transcriptomics",
        min_coverage=0.5,
    )

    assert features == ["f0"]
    assert X.shape == (3, 1)
    assert y.tolist() == [0, 0, 1]


