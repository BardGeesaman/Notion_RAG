"""Tests for Network Hub API endpoints (PPI, expression overlay)."""

import pytest
from uuid import uuid4
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient
from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user


class TestPPIEndpoints:
    @pytest.fixture
    def client(self):
        mock_user = MagicMock()
        mock_user.id = uuid4()
        app.dependency_overrides[get_current_user] = lambda: mock_user
        try:
            yield TestClient(app)
        finally:
            app.dependency_overrides.clear()

    @patch("amprenta_rag.api.routers.graph.get_interaction_partners")
    @patch("amprenta_rag.api.routers.graph.interactions_to_cytoscape")
    def test_get_ppi_single_gene(self, mock_cyto, mock_partners, client):
        mock_partners.return_value = []
        mock_cyto.return_value = ([], [])
        
        response = client.get("/api/graph/ppi/TP53")
        assert response.status_code == 200
        data = response.json()
        assert "nodes" in data
        assert "edges" in data

    @patch("amprenta_rag.api.routers.graph.get_interactions")
    @patch("amprenta_rag.api.routers.graph.interactions_to_cytoscape")
    def test_post_ppi_network(self, mock_cyto, mock_interactions, client):
        mock_interactions.return_value = []
        mock_cyto.return_value = ([], [])
        
        response = client.post("/api/graph/ppi/network", json={
            "proteins": ["TP53", "MDM2"],
            "species": 9606,
            "min_score": 400
        })
        assert response.status_code == 200


class TestExpressionOverlayEndpoint:
    @pytest.fixture
    def client(self):
        mock_user = MagicMock()
        mock_user.id = uuid4()
        app.dependency_overrides[get_current_user] = lambda: mock_user
        try:
            yield TestClient(app)
        finally:
            app.dependency_overrides.clear()

    @patch("amprenta_rag.services.expression_overlay.get_expression_overlay")
    def test_expression_overlay(self, mock_overlay, client):
        mock_overlay.return_value = (
            {"GENE1": "#b2182b"},
            {"GENE1": 2.5}
        )
        
        dataset_id = uuid4()
        response = client.post(
            f"/api/v1/datasets/{dataset_id}/expression-overlay",
            json={"gene_symbols": ["GENE1"], "colormap": "diverging"}
        )
        assert response.status_code == 200
        data = response.json()
        assert "node_colors" in data
