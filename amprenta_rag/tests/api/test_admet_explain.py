from __future__ import annotations

import asyncio
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient
from httpx import ASGITransport, AsyncClient

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


class TestAsyncADMETAPI:
    """Test async execution of ADMET API endpoints."""

    @pytest.mark.asyncio
    async def test_predict_async(self):
        """Test ADMET prediction endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.admet._sync_predict_admet') as mock_predict:
            # Mock the prediction function
            from amprenta_rag.api.schemas import ADMETPredictResponse, ADMETCompoundPrediction, ADMETEndpointPrediction
            
            mock_response = ADMETPredictResponse(
                results=[
                    ADMETCompoundPrediction(
                        smiles="CCO",
                        predictions={
                            "herg": ADMETEndpointPrediction(
                                mean=0.25,
                                std=0.05,
                                calibrated=True
                            )
                        },
                        error=None
                    )
                ],
                model_info={
                    "endpoints": ["herg"],
                    "models": {
                        "herg": {
                            "name": "test_model",
                            "status": "active"
                        }
                    }
                }
            )
            mock_predict.return_value = mock_response
            
            transport = ASGITransport(app=app)
            async with AsyncClient(transport=transport, base_url="http://test") as client:
                response = await client.post(
                    "/api/admet/predict",
                    json={
                        "smiles": ["CCO"],
                        "endpoints": ["herg"],
                        "include_uncertainty": True
                    }
                )
                
                assert response.status_code == 200
                data = response.json()
                assert len(data["results"]) == 1
                assert data["results"][0]["smiles"] == "CCO"
                assert data["results"][0]["predictions"]["herg"]["mean"] == 0.25
                assert data["model_info"]["endpoints"] == ["herg"]
                
                # Verify async execution
                mock_predict.assert_called_once()
                call_args = mock_predict.call_args[0]
                assert call_args[0].smiles == ["CCO"]
                assert call_args[0].endpoints == ["herg"]
                assert call_args[0].include_uncertainty == True

    @pytest.mark.asyncio
    async def test_explain_async(self):
        """Test ADMET explanation endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.admet._sync_explain_admet') as mock_explain:
            # Mock the explanation function
            from amprenta_rag.api.schemas import ADMETExplainResponse
            
            mock_response = ADMETExplainResponse(
                smiles="CCO",
                endpoint="herg",
                prediction={
                    "mean": 0.12,
                    "std": 0.01,
                    "calibrated": False
                },
                shap={
                    "top_features": [
                        {"name": "AsyncFeature1", "value": 0.8, "rank": 1},
                        {"name": "AsyncFeature2", "value": 0.6, "rank": 2}
                    ],
                    "base_value": 0.5,
                    "other_sum": 0.1
                },
                error=None
            )
            mock_explain.return_value = mock_response
            
            transport = ASGITransport(app=app)
            async with AsyncClient(transport=transport, base_url="http://test") as client:
                response = await client.post(
                    "/api/admet/explain",
                    json={
                        "smiles": "CCO",
                        "endpoint": "herg",
                        "top_k": 5
                    }
                )
                
                assert response.status_code == 200
                data = response.json()
                assert data["smiles"] == "CCO"
                assert data["endpoint"] == "herg"
                assert data["prediction"]["mean"] == 0.12
                assert len(data["shap"]["top_features"]) == 2
                assert data["shap"]["top_features"][0]["name"] == "AsyncFeature1"
                
                # Verify async execution
                mock_explain.assert_called_once()
                call_args = mock_explain.call_args[0]
                assert call_args[0].smiles == "CCO"
                assert call_args[0].endpoint == "herg"
                assert call_args[0].top_k == 5

    @pytest.mark.asyncio
    async def test_concurrent_admet(self):
        """Test multiple simultaneous ADMET requests."""
        with patch('amprenta_rag.api.routers.admet._sync_predict_admet') as mock_predict:
            with patch('amprenta_rag.api.routers.admet._sync_explain_admet') as mock_explain:
                # Mock prediction function
                from amprenta_rag.api.schemas import (
                    ADMETPredictResponse, ADMETCompoundPrediction, ADMETEndpointPrediction,
                    ADMETExplainResponse
                )
                
                mock_predict.return_value = ADMETPredictResponse(
                    results=[
                        ADMETCompoundPrediction(
                            smiles="CCC",
                            predictions={
                                "cyp3a4": ADMETEndpointPrediction(mean=0.75, calibrated=False)
                            },
                            error=None
                        )
                    ],
                    model_info={"endpoints": ["cyp3a4"], "models": {}}
                )
                
                mock_explain.return_value = ADMETExplainResponse(
                    smiles="CCO",
                    endpoint="herg",
                    prediction={"mean": 0.3},
                    shap={"top_features": [], "base_value": 0.5},
                    error=None
                )
                
                transport = ASGITransport(app=app)
                async with AsyncClient(transport=transport, base_url="http://test") as client:
                    # Define async request functions
                    async def predict_request():
                        return await client.post(
                            "/api/admet/predict",
                            json={
                                "smiles": ["CCC"],
                                "endpoints": ["cyp3a4"],
                                "include_uncertainty": False
                            }
                        )
                    
                    async def explain_request():
                        return await client.post(
                            "/api/admet/explain",
                            json={
                                "smiles": "CCO",
                                "endpoint": "herg",
                                "top_k": 3
                            }
                        )
                    
                    # Make 2 concurrent requests (predict + explain)
                    tasks = [
                        predict_request(),
                        explain_request()
                    ]
                    
                    results = await asyncio.gather(*tasks)
                    
                    # All requests should succeed
                    assert len(results) == 2
                    assert all(r.status_code == 200 for r in results)
                    
                    # Verify responses
                    predict_data = results[0].json()
                    explain_data = results[1].json()
                    
                    assert predict_data["results"][0]["smiles"] == "CCC"
                    assert predict_data["results"][0]["predictions"]["cyp3a4"]["mean"] == 0.75
                    
                    assert explain_data["smiles"] == "CCO"
                    assert explain_data["endpoint"] == "herg"
                    assert explain_data["prediction"]["mean"] == 0.3
                    
                    # Verify all calls were made
                    assert mock_predict.call_count == 1
                    assert mock_explain.call_count == 1



