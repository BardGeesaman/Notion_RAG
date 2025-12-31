from __future__ import annotations

import asyncio
from unittest.mock import patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient
from httpx import ASGITransport, AsyncClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_get_methods_returns_list():
    resp = client.get("/api/biomarker/methods")
    assert resp.status_code == 200
    data = resp.json()
    assert isinstance(data, list)
    assert set(data) >= {"statistical", "stability", "importance"}


def test_discover_endpoint_returns_consensus(monkeypatch):
    from amprenta_rag.api.routers import biomarker as biomarker_router

    class FakeSvc:
        def discover(self, experiment_id, group1, group2, methods=None):  # noqa: ANN001
            assert group1 == ["S1"]
            assert group2 == ["S2"]
            return {
                "consensus_ranking": [
                    {"feature": "A", "avg_rank": 1.0, "methods": 1},
                    {"feature": "B", "avg_rank": 2.0, "methods": 1},
                ],
                "method_results": {
                    "statistical": [
                        {"feature": "A", "t_stat": 1.0, "p_value": 0.01, "p_adj": 0.01},
                        {"feature": "B", "t_stat": 0.5, "p_value": 0.5, "p_adj": 0.5},
                    ]
                },
            }

    monkeypatch.setattr(biomarker_router, "BiomarkerDiscoveryService", lambda: FakeSvc())

    exp_id = uuid4()
    resp = client.post(
        "/api/biomarker/discover",
        json={
            "experiment_id": str(exp_id),
            "group1_samples": ["S1"],
            "group2_samples": ["S2"],
            "methods": ["statistical"],
            "fdr_threshold": 0.05,
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert "consensus_ranking" in data
    assert "method_results" in data
    # FDR filtering should keep only feature A
    assert [r["feature"] for r in data["consensus_ranking"]] == ["A"]
    assert [r["feature"] for r in data["method_results"]["statistical"]] == ["A"]


class TestAsyncBiomarkerAPI:
    """Test async execution of biomarker API endpoints."""

    @pytest.mark.asyncio
    async def test_discover_async(self):
        """Test biomarker discovery endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.biomarker._sync_discover_biomarkers') as mock_discover:
            # Mock the discovery function
            from amprenta_rag.api.schemas import BiomarkerDiscoverResponse, BiomarkerFeature
            
            mock_response = BiomarkerDiscoverResponse(
                consensus_ranking=[
                    BiomarkerFeature(
                        feature="AsyncFeature1",
                        avg_rank=1.0,
                        methods=2,
                        score=0.95
                    ),
                    BiomarkerFeature(
                        feature="AsyncFeature2", 
                        avg_rank=2.0,
                        methods=2,
                        score=0.85
                    )
                ],
                method_results={
                    "statistical": [
                        {
                            "feature": "AsyncFeature1",
                            "t_stat": 3.5,
                            "p_value": 0.001,
                            "p_adj": 0.005,
                            "fold_change": 2.1
                        },
                        {
                            "feature": "AsyncFeature2",
                            "t_stat": 2.8,
                            "p_value": 0.01,
                            "p_adj": 0.02,
                            "fold_change": 1.8
                        }
                    ],
                    "stability": [
                        {
                            "feature": "AsyncFeature1",
                            "stability_score": 0.92
                        },
                        {
                            "feature": "AsyncFeature2",
                            "stability_score": 0.88
                        }
                    ]
                }
            )
            mock_discover.return_value = mock_response
            
            transport = ASGITransport(app=app)
            async with AsyncClient(transport=transport, base_url="http://test") as client:
                exp_id = uuid4()
                response = await client.post(
                    "/api/biomarker/discover",
                    json={
                        "experiment_id": str(exp_id),
                        "group1_samples": ["Group1_S1", "Group1_S2"],
                        "group2_samples": ["Group2_S1", "Group2_S2"],
                        "methods": ["statistical", "stability"],
                        "fdr_threshold": 0.05
                    }
                )
                
                assert response.status_code == 200
                data = response.json()
                assert len(data["consensus_ranking"]) == 2
                assert data["consensus_ranking"][0]["feature"] == "AsyncFeature1"
                assert data["consensus_ranking"][0]["avg_rank"] == 1.0
                assert data["consensus_ranking"][1]["feature"] == "AsyncFeature2"
                assert "statistical" in data["method_results"]
                assert "stability" in data["method_results"]
                
                # Verify async execution
                mock_discover.assert_called_once()
                call_args = mock_discover.call_args[0][0]
                assert str(call_args.experiment_id) == str(exp_id)
                assert call_args.group1_samples == ["Group1_S1", "Group1_S2"]
                assert call_args.group2_samples == ["Group2_S1", "Group2_S2"]
                assert call_args.methods == ["statistical", "stability"]
                assert call_args.fdr_threshold == 0.05

    @pytest.mark.asyncio
    async def test_concurrent_discover(self):
        """Test multiple simultaneous biomarker discovery requests."""
        with patch('amprenta_rag.api.routers.biomarker._sync_discover_biomarkers') as mock_discover:
            # Mock function to return different results for different experiments
            from amprenta_rag.api.schemas import BiomarkerDiscoverResponse, BiomarkerFeature
            
            def mock_discover_func(request):
                exp_id_str = str(request.experiment_id)
                feature_name = f"Feature_Exp_{exp_id_str[-4:]}"  # Use last 4 chars of UUID
                
                return BiomarkerDiscoverResponse(
                    consensus_ranking=[
                        BiomarkerFeature(
                            feature=feature_name,
                            avg_rank=1.0,
                            methods=1,
                            score=0.9
                        )
                    ],
                    method_results={
                        "statistical": [
                            {
                                "feature": feature_name,
                                "t_stat": 2.5,
                                "p_value": 0.01,
                                "p_adj": 0.02
                            }
                        ]
                    }
                )
            
            mock_discover.side_effect = mock_discover_func
            
            transport = ASGITransport(app=app)
            async with AsyncClient(transport=transport, base_url="http://test") as client:
                # Define async request functions
                async def discover_for_experiment(exp_suffix: str):
                    exp_id = uuid4()
                    return await client.post(
                        "/api/biomarker/discover",
                        json={
                            "experiment_id": str(exp_id),
                            "group1_samples": [f"Exp{exp_suffix}_G1_S1"],
                            "group2_samples": [f"Exp{exp_suffix}_G2_S1"],
                            "methods": ["statistical"]
                            # fdr_threshold will use default value 0.05
                        }
                    )
                
                # Make 3 concurrent discovery requests
                tasks = [
                    discover_for_experiment("A"),
                    discover_for_experiment("B"),
                    discover_for_experiment("C")
                ]
                
                results = await asyncio.gather(*tasks)
                
                # All requests should succeed
                assert len(results) == 3
                assert all(r.status_code == 200 for r in results)
                
                # Verify responses are different (different experiment IDs lead to different features)
                features = []
                for result in results:
                    data = result.json()
                    assert len(data["consensus_ranking"]) == 1
                    feature_name = data["consensus_ranking"][0]["feature"]
                    features.append(feature_name)
                    assert feature_name.startswith("Feature_Exp_")
                    assert "statistical" in data["method_results"]
                
                # All features should be different
                assert len(set(features)) == 3
                
                # Verify all calls were made
                assert mock_discover.call_count == 3


