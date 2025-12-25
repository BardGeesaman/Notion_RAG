from __future__ import annotations

from typing import Any, Dict

import numpy as np
import pytest

requests = pytest.importorskip("requests")
rdkit = pytest.importorskip("rdkit")
_ = (requests, rdkit)


def test_normalize_to_nm_various_units():
    from amprenta_rag.ml.qsar.datasets import _normalize_to_nm

    assert _normalize_to_nm(1.0, "nM") == 1.0
    assert _normalize_to_nm(1.0, "uM") == 1000.0
    assert _normalize_to_nm(1.0, "μM") == 1000.0
    assert _normalize_to_nm(1.0, "mM") == 1_000_000.0
    assert _normalize_to_nm(1.0, "M") == 1_000_000_000.0
    assert _normalize_to_nm(1.0, "weird") is None


def test_aggregate_ic50_median():
    from amprenta_rag.ml.qsar.datasets import _aggregate_ic50

    assert _aggregate_ic50([1.0, 2.0, 100.0]) == 2.0
    assert _aggregate_ic50([10.0, 20.0]) == 15.0


def test_load_target_chembl(monkeypatch):
    from amprenta_rag.ml.qsar import datasets as ds
    from amprenta_rag.ml.admet import predictor as predmod

    def fake_get(url: str, params: Dict[str, Any], timeout: int = 30):  # noqa: ANN001
        class R:
            def raise_for_status(self):  # noqa: ANN001
                return None

            def json(self):  # noqa: ANN001
                if "activity.json" in url:
                    return {
                        "activities": [
                            {
                                "molecule_chembl_id": "CHEMBL1",
                                "standard_value": 0.5,
                                "standard_units": "uM",
                            },
                            {
                                "molecule_chembl_id": "CHEMBL1",
                                "standard_value": 500,
                                "standard_units": "nM",
                            },
                            {
                                "molecule_chembl_id": "CHEMBL2",
                                "standard_value": 2000,
                                "standard_units": "nM",
                            },
                        ]
                    }
                if "molecule/CHEMBL1.json" in url:
                    return {"molecule_structures": {"canonical_smiles": "CCO"}}
                if "molecule/CHEMBL2.json" in url:
                    return {"molecule_structures": {"canonical_smiles": "CCN"}}
                return {}

        return R()

    monkeypatch.setattr(ds, "_chembl_get_json", lambda url, params: fake_get(url, params).json())
    monkeypatch.setattr(predmod.ADMETPredictor, "_get_features", lambda self, s: np.zeros((2054,), dtype=np.float32))

    X, y, smiles, meta = ds.fetch_chembl_target_data("CHEMBL203", threshold_nm=1000, limit=5000)
    assert X.shape == (2, 2054)
    assert y.tolist() == [1, 0]  # CCO active (median 500 nM), CCN inactive (2000 nM)
    assert smiles == ["CCO", "CCN"]
    assert meta["source"] == "chembl"


def test_load_target_validates_criteria(monkeypatch):
    from amprenta_rag.ml.qsar.datasets import TargetDatasetLoader
    from amprenta_rag.ml.qsar import datasets as ds

    monkeypatch.setattr(
        ds,
        "fetch_chembl_target_data",
        lambda *args, **kwargs: (np.zeros((10, 2054), dtype=np.float32), np.zeros((10,), dtype=np.int64), ["C"] * 10, {}),
    )

    loader = TargetDatasetLoader()
    with pytest.raises(ValueError, match="Insufficient compounds"):
        loader.load_target("CHEMBL203", source="chembl", min_compounds=100)



