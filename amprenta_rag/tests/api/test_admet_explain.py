from __future__ import annotations

from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


client = TestClient(app)


def test_explain_endpoint_returns_shap(monkeypatch):
    import amprenta_rag.api.routers.admet as admet_router

    class _P:
        def predict_with_uncertainty(self, smiles_list, endpoints, include_shap, shap_top_k):  # noqa: ANN001
            assert smiles_list == ["CCO"]
            assert endpoints == ["herg"]
            assert include_shap is True
            assert shap_top_k == 5
            return [
                {
                    "smiles": "CCO",
                    "predictions": {
                        "herg": {
                            "mean": 0.12,
                            "std": 0.01,
                            "ci_low": 0.10,
                            "ci_high": 0.14,
                            "in_domain": True,
                            "similarity": 0.9,
                            "calibrated": False,
                            "shap": {
                                "top_features": [{"name": "MorganBit_0000", "value": 1.0, "rank": 1}],
                                "other_sum": 0.0,
                                "base_value": 0.5,
                            },
                        }
                    },
                    "error": None,
                }
            ]

    monkeypatch.setattr(admet_router, "get_admet_predictor", lambda: _P())

    resp = client.post("/api/admet/explain", json={"smiles": "CCO", "endpoint": "herg", "top_k": 5})
    assert resp.status_code == 200
    data = resp.json()
    assert data["smiles"] == "CCO"
    assert data["endpoint"] == "herg"
    assert data["prediction"]["mean"] == 0.12
    assert data["shap"]["top_features"][0]["name"] == "MorganBit_0000"
    assert data["error"] is None


def test_explain_endpoint_invalid_smiles(monkeypatch):
    import amprenta_rag.api.routers.admet as admet_router

    class _P:
        def predict_with_uncertainty(self, smiles_list, endpoints, include_shap, shap_top_k):  # noqa: ANN001
            return [{"smiles": smiles_list[0], "predictions": {}, "error": "Invalid SMILES"}]

    monkeypatch.setattr(admet_router, "get_admet_predictor", lambda: _P())

    resp = client.post("/api/admet/explain", json={"smiles": "BAD", "endpoint": "herg", "top_k": 10})
    assert resp.status_code == 200
    data = resp.json()
    assert data["error"] == "Invalid SMILES"
    assert data["prediction"] == {}
    assert data["shap"] is None


def test_explain_endpoint_shap_not_installed(monkeypatch):
    import amprenta_rag.api.routers.admet as admet_router

    class _P:
        def predict_with_uncertainty(self, smiles_list, endpoints, include_shap, shap_top_k):  # noqa: ANN001
            return [
                {
                    "smiles": smiles_list[0],
                    "predictions": {
                        "herg": {
                            "mean": 0.5,
                            "std": 0.0,
                            "ci_low": 0.5,
                            "ci_high": 0.5,
                            "in_domain": None,
                            "similarity": None,
                            "calibrated": False,
                            "shap": {"error": "shap not installed"},
                        }
                    },
                    "error": None,
                }
            ]

    monkeypatch.setattr(admet_router, "get_admet_predictor", lambda: _P())

    resp = client.post("/api/admet/explain", json={"smiles": "CCO", "endpoint": "herg", "top_k": 10})
    assert resp.status_code == 200
    data = resp.json()
    assert data["error"] is None
    assert data["shap"]["error"] == "shap not installed"



