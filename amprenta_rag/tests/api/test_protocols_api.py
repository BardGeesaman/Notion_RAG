"""
Unit tests for protocols API endpoints.

Tests protocol history, diffing, and deviation auditing with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestProtocolHistory:
    """Tests for GET /api/v1/protocols/{protocol_id}/history endpoint."""

    @patch("amprenta_rag.api.routers.protocols.get_protocol_history")
    def test_protocol_history_success(self, mock_get_history):
        """Test successful protocol history retrieval."""
        protocol_id = uuid4()
        parent_id = uuid4()
        
        # Mock history items matching ProtocolHistoryItem schema
        mock_item1 = MagicMock()
        mock_item1.asdict.return_value = {
            "protocol_id": protocol_id,
            "version": 1,
            "parent_id": None,
            "name": "Initial Protocol",
        }
        
        mock_item2 = MagicMock()
        mock_item2.asdict.return_value = {
            "protocol_id": protocol_id,
            "version": 2,
            "parent_id": parent_id,
            "name": "Updated Protocol",
        }
        
        mock_get_history.return_value = [mock_item1, mock_item2]
        
        response = client.get(f"/api/v1/protocols/{protocol_id}/history")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        assert data[0]["version"] == 1
        assert data[1]["version"] == 2
        assert data[0]["name"] == "Initial Protocol"

    @patch("amprenta_rag.api.routers.protocols.get_protocol_history")
    def test_protocol_history_empty(self, mock_get_history):
        """Test protocol history when no versions exist."""
        protocol_id = uuid4()
        
        mock_get_history.return_value = []
        
        response = client.get(f"/api/v1/protocols/{protocol_id}/history")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []


class TestProtocolDiff:
    """Tests for GET /api/v1/protocols/{protocol_id}/diff/{other_id} endpoint."""

    @patch("amprenta_rag.api.routers.protocols.diff_protocols")
    @patch("amprenta_rag.api.routers.protocols.db_session")
    def test_protocol_diff_success(self, mock_db_session, mock_diff):
        """Test successful protocol diff."""
        protocol1_id = uuid4()
        protocol2_id = uuid4()
        
        # Mock protocols
        mock_protocol1 = MagicMock()
        mock_protocol1.id = protocol1_id
        
        mock_protocol2 = MagicMock()
        mock_protocol2.id = protocol2_id
        
        # Mock database session
        mock_db = MagicMock()
        
        def mock_query_side_effect(*args):
            query = MagicMock()
            
            def filter_side_effect(*filter_args):
                # Return different protocols based on filter
                first_mock = MagicMock()
                # Crude check: return protocol1 first, then protocol2
                if not hasattr(mock_query_side_effect, 'call_count'):
                    mock_query_side_effect.call_count = 0
                
                if mock_query_side_effect.call_count == 0:
                    first_mock.first.return_value = mock_protocol1
                else:
                    first_mock.first.return_value = mock_protocol2
                
                mock_query_side_effect.call_count += 1
                return first_mock
            
            query.filter.side_effect = filter_side_effect
            return query
        
        mock_db.query.side_effect = mock_query_side_effect
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock diff result matching ProtocolDiff schema
        mock_diff_result = MagicMock()
        mock_diff_result.asdict.return_value = {
            "protocol_id": protocol1_id,
            "other_id": protocol2_id,
            "added_steps": [],
            "removed_steps": [],
            "changed_steps": [{"step_id": "step_3", "old": "Mix for 5min", "new": "Mix for 10min"}],
            "materials_added": [],
            "materials_removed": [],
            "parameters_changed": [],
        }
        mock_diff.return_value = mock_diff_result
        
        response = client.get(f"/api/v1/protocols/{protocol1_id}/diff/{protocol2_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["protocol_id"] == str(protocol1_id)
        assert data["other_id"] == str(protocol2_id)
        assert len(data["changed_steps"]) == 1

    @patch("amprenta_rag.api.routers.protocols.db_session")
    def test_protocol_diff_not_found(self, mock_db_session):
        """Test diff with non-existent protocol."""
        protocol1_id = uuid4()
        protocol2_id = uuid4()
        
        # Mock database session - first protocol not found
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get(f"/api/v1/protocols/{protocol1_id}/diff/{protocol2_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()


class TestExperimentDeviations:
    """Tests for GET /api/v1/experiments/{experiment_id}/deviations endpoint."""

    @patch("amprenta_rag.api.routers.protocols.audit_deviations")
    def test_experiment_deviations_success(self, mock_audit):
        """Test successful deviation audit."""
        experiment_id = uuid4()
        protocol_id = uuid4()
        
        # Mock deviation reports matching DeviationReport schema
        mock_report1 = MagicMock()
        mock_report1.asdict.return_value = {
            "experiment_id": experiment_id,
            "protocol_id": protocol_id,
            "protocol_name": "Standard Protocol",
            "protocol_version": 2,
            "deviations": [
                {"step": "Step 3", "expected": "Mix for 10min", "actual": "Mix for 5min", "severity": "minor"}
            ],
        }
        
        mock_report2 = MagicMock()
        mock_report2.asdict.return_value = {
            "experiment_id": experiment_id,
            "protocol_id": protocol_id,
            "protocol_name": "Standard Protocol",
            "protocol_version": 2,
            "deviations": [
                {"step": "Step 5", "expected": "37°C", "actual": "25°C", "severity": "major"}
            ],
        }
        
        mock_audit.return_value = [mock_report1, mock_report2]
        
        response = client.get(f"/api/v1/experiments/{experiment_id}/deviations")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        assert data[0]["experiment_id"] == str(experiment_id)
        assert data[0]["protocol_name"] == "Standard Protocol"
        assert len(data[0]["deviations"]) == 1

    @patch("amprenta_rag.api.routers.protocols.audit_deviations")
    def test_experiment_deviations_none(self, mock_audit):
        """Test deviation audit when no deviations found."""
        experiment_id = uuid4()
        
        mock_audit.return_value = []
        
        response = client.get(f"/api/v1/experiments/{experiment_id}/deviations")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []

