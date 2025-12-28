"""API tests for Bayesian endpoints."""
from __future__ import annotations

import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock

from amprenta_rag.api.main import app

client = TestClient(app)


def test_dose_response_endpoint_success():
    """Test dose-response endpoint with valid data."""
    mock_result = {
        "ec50_mean": 1.5,
        "ec50_ci": (0.8, 2.2),
        "hill_slope": 1.2,
    }
    with patch("amprenta_rag.api.routers.bayesian.fit_bayesian_dose_response", return_value=mock_result):
        response = client.post("/api/v1/bayesian-optimization/dose-response", json={
            "concentrations": [0.1, 0.3, 1.0, 3.0, 10.0, 30.0],
            "responses": [5, 15, 45, 70, 90, 95],
        })
        assert response.status_code == 200
        data = response.json()
        assert "ec50_mean" in data


def test_dose_response_endpoint_with_custom_prior():
    """Test dose-response endpoint with custom prior config."""
    mock_result = {"ec50_mean": 2.0, "ec50_ci": (1.0, 3.0), "hill_slope": 1.0}
    with patch("amprenta_rag.api.routers.bayesian.fit_bayesian_dose_response", return_value=mock_result):
        response = client.post("/api/v1/bayesian-optimization/dose-response", json={
            "concentrations": [0.1, 0.3, 1.0, 3.0, 10.0, 30.0],
            "responses": [5, 15, 45, 70, 90, 95],
            "prior_config": {
                "ec50_prior_mean": -6.0,
                "ec50_prior_sd": 1.5,
                "hill_prior_mean": 1.0,
                "hill_prior_sd": 1.0,
            }
        })
        assert response.status_code == 200


def test_dose_response_validation_error():
    """Test dose-response endpoint with mismatched lengths."""
    response = client.post("/api/v1/bayesian-optimization/dose-response", json={
        "concentrations": [0.1, 0.3, 1.0],
        "responses": [5, 15],  # Mismatched length
    })
    assert response.status_code in [400, 500]  # Either validation or internal error


def test_variable_selection_endpoint_success():
    """Test variable selection endpoint with valid data."""
    mock_result = {
        "pips": {"gene_A": 0.8, "gene_B": 0.3},
        "included_features": ["gene_A"],
        "coefficients": {"gene_A": 1.2, "gene_B": 0.1},
        "n_included": 1,
        "pip_threshold": 0.5,
    }
    with patch("amprenta_rag.analysis.bayesian_variable_selection.compute_posterior_inclusion_probabilities", return_value=mock_result):
        response = client.post("/api/v1/bayesian-optimization/variable-selection", json={
            "features": [
                {"name": "gene_A", "values": [1.0, 2.0, 3.0, 4.0, 5.0]},
                {"name": "gene_B", "values": [0.5, 1.5, 2.5, 3.5, 4.5]},
            ],
            "response": [10, 20, 30, 40, 50],
            "pip_threshold": 0.5,
        })
        assert response.status_code == 200
        data = response.json()
        assert "pips" in data
        assert "included_features" in data


def test_variable_selection_validation_error():
    """Test variable selection with invalid features."""
    response = client.post("/api/v1/bayesian-optimization/variable-selection", json={
        "features": [],  # Empty features
        "response": [10, 20, 30],
    })
    assert response.status_code == 400 or response.status_code == 500
