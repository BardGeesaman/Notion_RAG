from __future__ import annotations

from typing import Any, Dict

import pandas as pd
import pytest

# Graceful skipping in minimal test environments
papermill = pytest.importorskip("papermill")
apscheduler = pytest.importorskip("apscheduler")
sklearn = pytest.importorskip("sklearn")

from amprenta_rag.notebook import automl_helpers as ah


def test_load_dataset_returns_dataframe(monkeypatch, tmp_path):
    p = tmp_path / "d.csv"
    p.write_text("x,y\n1,0\n2,1\n", encoding="utf-8")

    monkeypatch.setattr(ah, "_get_dataset_metadata", lambda dataset_id: {"file_paths": [str(p)]})
    df = ah.load_dataset_as_dataframe("00000000-0000-0000-0000-000000000000")
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["x", "y"]
    assert len(df) == 2


def test_register_model_creates_entry(monkeypatch):
    calls: list[Dict[str, Any]] = []

    class _Reg:
        def register_model(self, **kwargs):
            calls.append(kwargs)
            return {"ok": True}

    monkeypatch.setattr(ah, "get_registry", lambda: _Reg())
    out = ah.register_trained_model(
        model={"m": 1},
        name="demo_model",
        metrics={"auc": 0.9},
        dataset_id="00000000-0000-0000-0000-000000000000",
        model_type="automl_classification",
        framework="xgboost",
    )
    assert out["ok"] is True
    assert calls and calls[0]["name"] == "demo_model"
    assert calls[0]["metrics"]["auc"] == 0.9


def test_classification_report_structure():
    rep = ah.generate_classification_report([0, 1, 1], [0, 1, 0], y_proba=[0.1, 0.9, 0.2])
    assert "accuracy" in rep
    assert "auc" in rep


def test_regression_report_structure():
    rep = ah.generate_regression_report([1.0, 2.0, 3.0], [1.1, 1.9, 2.8])
    assert set(rep.keys()) == {"rmse", "mae", "r2"}


