"""
Unit tests for quality API endpoints.

Tests dataset quality checking with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestQualitySummary:
    """Tests for GET /api/v1/quality/summary endpoint."""

    @patch("amprenta_rag.api.routers.quality.get_quality_summary")
    def test_quality_summary_success(self, mock_get_summary):
        """Test successful quality summary retrieval."""
        # Match QualitySummary schema: total, high, medium, low
        mock_get_summary.return_value = {
            "total": 100,
            "high": 75,
            "medium": 20,
            "low": 5,
        }
        
        response = client.get("/api/v1/quality/summary")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 100
        assert data["high"] == 75
        assert data["medium"] == 20
        assert data["low"] == 5


class TestQualityDatasets:
    """Tests for GET /api/v1/quality/datasets endpoint."""

    @patch("amprenta_rag.api.routers.quality.scan_all_datasets")
    def test_quality_datasets_success(self, mock_scan):
        """Test successful dataset quality scan."""
        dataset1_id = uuid4()
        dataset2_id = uuid4()
        
        # Mock quality reports
        mock_report1 = MagicMock()
        mock_report1.dataset_id = dataset1_id
        mock_report1.dataset_name = "Dataset 1"
        mock_report1.score = 90.0
        mock_report1.status = "good"
        mock_report1.issues = []
        mock_report1.metrics = {"completeness": 0.95}
        
        mock_report2 = MagicMock()
        mock_report2.dataset_id = dataset2_id
        mock_report2.dataset_name = "Dataset 2"
        mock_report2.score = 45.0
        mock_report2.status = "poor"
        mock_report2.issues = ["High missing rate"]
        mock_report2.metrics = {"completeness": 0.40}
        
        mock_scan.return_value = [mock_report1, mock_report2]
        
        response = client.get("/api/v1/quality/datasets")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        assert data[0]["score"] == 90.0
        assert data[1]["score"] == 45.0

    @patch("amprenta_rag.api.routers.quality.scan_all_datasets")
    def test_quality_datasets_empty(self, mock_scan):
        """Test quality datasets when no datasets exist."""
        mock_scan.return_value = []
        
        response = client.get("/api/v1/quality/datasets")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []


class TestQualityAlerts:
    """Tests for GET /api/v1/quality/alerts endpoint."""

    @patch("amprenta_rag.api.routers.quality.get_low_quality_datasets")
    def test_quality_alerts_default_threshold(self, mock_get_low):
        """Test quality alerts with default threshold."""
        dataset_id = uuid4()
        
        # Mock low quality report
        mock_report = MagicMock()
        mock_report.dataset_id = dataset_id
        mock_report.dataset_name = "Low Quality Dataset"
        mock_report.score = 30.0
        mock_report.status = "critical"
        mock_report.issues = ["High missing rate", "Outliers detected"]
        mock_report.metrics = {"completeness": 0.25}
        
        mock_get_low.return_value = [mock_report]
        
        response = client.get("/api/v1/quality/alerts")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["score"] == 30.0
        assert data[0]["status"] == "critical"

    @patch("amprenta_rag.api.routers.quality.get_low_quality_datasets")
    def test_quality_alerts_custom_threshold(self, mock_get_low):
        """Test quality alerts with custom threshold."""
        mock_get_low.return_value = []
        
        response = client.get("/api/v1/quality/alerts?threshold=70.0")
        
        assert response.status_code == 200
        # Verify threshold was passed (check mock was called with it)
        mock_get_low.assert_called_once_with(threshold=70.0)

    @patch("amprenta_rag.api.routers.quality.get_low_quality_datasets")
    def test_quality_alerts_empty(self, mock_get_low):
        """Test quality alerts when no low quality datasets exist."""
        mock_get_low.return_value = []
        
        response = client.get("/api/v1/quality/alerts")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []

