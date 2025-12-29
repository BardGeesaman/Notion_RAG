"""Tests for experiment planner API endpoints."""

from __future__ import annotations

from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestPowerEndpoint:
    """Tests for POST /api/v1/planner/power endpoint."""

    @patch("amprenta_rag.api.routers.planner.calculate_sample_size")
    def test_calculate_power_success(self, mock_calc):
        """Test successful power calculation."""
        mock_calc.return_value = 64
        
        response = client.post(
            "/api/v1/planner/power",
            json={
                "effect_size": 0.5,
                "alpha": 0.05,
                "power": 0.80,
                "test_type": "t-test",
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["n_per_group"] == 64


class TestEffectSizeEndpoint:
    """Tests for POST /api/v1/planner/effect-size endpoint."""

    @patch("amprenta_rag.api.routers.planner.estimate_effect_size_from_data")
    def test_estimate_effect_size_success(self, mock_estimate):
        """Test successful effect size estimation."""
        mock_estimate.return_value = 0.75
        
        response = client.post(
            "/api/v1/planner/effect-size",
            json={
                "group1": [10.0, 12.0, 11.0],
                "group2": [8.0, 9.0, 7.0],
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["effect_size"] == 0.75


class TestPlatesEndpoint:
    """Tests for POST /api/v1/planner/plates endpoint."""

    @patch("amprenta_rag.api.routers.planner.calculate_plate_layout")
    def test_calculate_plates_success(self, mock_calc):
        """Test successful plate layout calculation."""
        mock_calc.return_value = {
            "plates_needed": 2,
            "wells_used": 100,
            "empty_wells": 92,
        }
        
        response = client.post(
            "/api/v1/planner/plates",
            json={"n": 100, "plate_format": 96},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["plates_needed"] == 2
        assert data["wells_used"] == 100


class TestCostEndpoint:
    """Tests for POST /api/v1/planner/cost endpoint."""

    @patch("amprenta_rag.api.routers.planner.estimate_experiment_cost")
    def test_estimate_cost_success(self, mock_estimate):
        """Test successful cost estimation."""
        mock_estimate.return_value = {
            "sample_cost": 1000.0,
            "overhead": 100.0,
            "total": 1100.0,
        }
        
        response = client.post(
            "/api/v1/planner/cost",
            json={
                "n": 100,
                "cost_per_sample": 10.0,
                "overhead_pct": 0.1,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 1100.0

    @patch("amprenta_rag.api.routers.planner.estimate_experiment_cost")
    def test_estimate_cost_custom_overhead(self, mock_estimate):
        """Test cost estimation with custom overhead."""
        mock_estimate.return_value = {
            "sample_cost": 500.0,
            "overhead": 75.0,
            "total": 575.0,
        }
        
        response = client.post(
            "/api/v1/planner/cost",
            json={
                "n": 50,
                "cost_per_sample": 10.0,
                "overhead_pct": 0.15,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["overhead"] == 75.0

