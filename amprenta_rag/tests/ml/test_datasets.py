from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

xgb = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")
_ = (xgb, sklearn)


def test_load_endpoint_returns_correct_shapes(monkeypatch, tmp_path: Path):
    from amprenta_rag.ml.admet.datasets import ChEMBLDatasetLoader
    from amprenta_rag.ml.admet import predictor as predmod

    data_dir = tmp_path / "admet"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "herg_train.parquet").write_text("placeholder")

    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "canonical_smiles": ["CCO", "CCN", "CCC"],
            "y": [0, 1, 0],
        }
    )

    monkeypatch.setattr(pd, "read_parquet", lambda p: df)
    monkeypatch.setattr(predmod.ADMETPredictor, "_get_features", lambda self, s: np.zeros((2054,), dtype=np.float32))

    loader = ChEMBLDatasetLoader(str(data_dir))
    X, y, smiles = loader.load_endpoint("herg", "train")
    assert X.shape == (3, 2054)
    assert y.shape == (3,)
    assert len(smiles) == 3


def test_load_endpoint_skips_invalid_smiles(monkeypatch, tmp_path: Path):
    from amprenta_rag.ml.admet.datasets import ChEMBLDatasetLoader
    from amprenta_rag.ml.admet import predictor as predmod

    data_dir = tmp_path / "admet"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "herg_train.parquet").write_text("placeholder")

    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "canonical_smiles": ["CCO", "NOT_A_SMILES", ""],
            "y": [0, 1, 0],
        }
    )

    def fake_features(self, smiles: str):  # noqa: ANN001
        if smiles == "NOT_A_SMILES":
            return None
        if not smiles:
            return None
        return np.ones((2054,), dtype=np.float32)

    monkeypatch.setattr(pd, "read_parquet", lambda p: df)
    monkeypatch.setattr(predmod.ADMETPredictor, "_get_features", fake_features)

    loader = ChEMBLDatasetLoader(str(data_dir))
    X, y, smiles = loader.load_endpoint("herg", "train")
    assert X.shape == (1, 2054)
    assert y.tolist() == [0.0]
    assert smiles == ["CCO"]


def test_get_summary_loads_json(tmp_path: Path):
    from amprenta_rag.ml.admet.datasets import ChEMBLDatasetLoader

    data_dir = tmp_path / "admet"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "summary.json").write_text(json.dumps({"hello": "world"}))

    loader = ChEMBLDatasetLoader(str(data_dir))
    out = loader.get_summary()
    assert out["hello"] == "world"


def test_list_endpoints_finds_parquets(tmp_path: Path):
    from amprenta_rag.ml.admet.datasets import ChEMBLDatasetLoader

    data_dir = tmp_path / "admet"
    data_dir.mkdir(parents=True, exist_ok=True)
    for fn in ("herg_train.parquet", "herg_test.parquet", "logs_train.parquet", "weird.txt"):
        (data_dir / fn).write_text("x")

    loader = ChEMBLDatasetLoader(str(data_dir))
    eps = loader.list_endpoints()
    assert eps == ["herg", "logs"]


