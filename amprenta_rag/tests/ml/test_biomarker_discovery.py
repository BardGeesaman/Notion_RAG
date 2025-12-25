from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4

import numpy as np


def test_consensus_ranking_averages_ranks():
    from amprenta_rag.ml.biomarker.discovery import BiomarkerDiscoveryService

    svc = BiomarkerDiscoveryService()
    consensus = svc._consensus_ranking(  # noqa: SLF001
        {
            "m1": ["A", "B", "C", "D"],
            "m2": ["B", "C", "A", "D"],
        }
    )
    # B: (2+1)/2 = 1.5 should rank ahead of A: (1+3)/2 = 2.0
    assert consensus[0]["feature"] == "B"
    ranks = {r["feature"]: r["avg_rank"] for r in consensus}
    assert ranks["B"] == 1.5
    assert ranks["A"] == 2.0
    assert ranks["D"] == 4.0


def test_discover_returns_expected_structure(monkeypatch):
    from amprenta_rag.ml.biomarker import discovery as dmod
    from amprenta_rag.ml.biomarker.discovery import BiomarkerDiscoveryService

    exp_id = uuid4()

    fake_exp = SimpleNamespace(
        id=exp_id,
        sample_groups={"control": ["S1"], "case": ["S2"]},
        datasets=[SimpleNamespace(id=uuid4(), omics_type="transcriptomics")],
    )

    class FakeQuery:
        def __init__(self, obj):
            self.obj = obj

        def filter(self, *args, **kwargs):  # noqa: ANN001
            return self

        def first(self):
            return self.obj

    class FakeDB:
        def query(self, *args, **kwargs):  # noqa: ANN001
            return FakeQuery(fake_exp)

    @contextmanager
    def fake_db_session():
        yield FakeDB()

    monkeypatch.setattr(dmod, "db_session", fake_db_session)

    # Avoid DB/datasets, provide synthetic data.
    X = np.zeros((10, 4), dtype=float)
    y = np.array([0] * 5 + [1] * 5, dtype=int)
    feat = ["A", "B", "C", "D"]
    monkeypatch.setattr(dmod, "load_biomarker_dataset", lambda **kwargs: (X, y, feat))

    # Mock method implementations to avoid sklearn/scipy.
    monkeypatch.setattr(dmod, "t_test_selection", lambda X, y, feature_names=None: [(f, 0.0, i / 100.0) for i, f in enumerate(feat)])
    monkeypatch.setattr(dmod, "fdr_correction", lambda pvals: ([True] * len(pvals), list(pvals)))

    class FakeStab:
        def __init__(self, n_bootstrap=50, threshold=0.0):  # noqa: ANN001
            self.threshold = threshold

        def fit(self, X, y, feature_names=None, random_state=None):  # noqa: ANN001
            self.feature_names_ = list(feature_names or feat)
            return self

        def get_ranked_features(self):
            return [(self.feature_names_[0], 1.0, 0.5), (self.feature_names_[1], 0.9, 0.2)]

    class FakeImp:
        def __init__(self, n_folds=5, n_estimators=100):  # noqa: ANN001
            pass

        def fit(self, X, y, feature_names=None, random_state=None):  # noqa: ANN001
            self.feature_names_ = list(feature_names or feat)
            return self

        def get_ranked_features(self):
            return [(self.feature_names_[1], 0.5, 0.1), (self.feature_names_[2], 0.4, 0.2)]

    monkeypatch.setattr(dmod, "StabilitySelector", FakeStab)
    monkeypatch.setattr(dmod, "CVFeatureImportance", FakeImp)

    svc = BiomarkerDiscoveryService()
    out = svc.discover(experiment_id=exp_id, group1=["S1"], group2=["S2"], methods=["statistical", "stability", "importance"])

    assert set(out.keys()) >= {"consensus_ranking", "method_results"}
    assert isinstance(out["consensus_ranking"], list)
    assert isinstance(out["method_results"], dict)
    assert "statistical" in out["method_results"]
    assert "stability" in out["method_results"]
    assert "importance" in out["method_results"]


