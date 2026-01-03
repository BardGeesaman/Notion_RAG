"""Tests for compound inventory API endpoints."""

from __future__ import annotations

from datetime import date, timedelta
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestCompoundSampleEndpoints:
    """Tests for compound sample API endpoints."""

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_create_compound_sample_success(self, mock_get_db, mock_get_user, mock_service):
        """Test successful compound sample creation."""
        # Mock user and db
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        mock_get_db.return_value = MagicMock()

        # Mock service response
        mock_sample = MagicMock()
        mock_sample.id = uuid4()
        mock_sample.barcode = "COMP-20250103-ABC123"
        mock_sample.compound_id = uuid4()
        mock_sample.quantity = 500.0
        mock_sample.unit = "µL"
        mock_sample.status = "available"
        mock_sample.sample_type = "compound"
        mock_service.create_compound_sample.return_value = mock_sample

        response = client.post(
            "/api/v1/inventory/samples",
            json={
                "compound_id": str(uuid4()),
                "quantity": 500.0,
                "quantity_unit": "µL",
                "concentration": 10.0,
                "concentration_unit": "mM",
                "solvent": "DMSO",
            }
        )

        assert response.status_code == 201
        data = response.json()
        assert data["barcode"] == "COMP-20250103-ABC123"
        assert data["quantity"] == 500.0

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_create_sample_barcode_collision(self, mock_get_db, mock_get_user, mock_service):
        """Test barcode collision error."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()
        
        mock_service.create_compound_sample.side_effect = ValueError("Barcode TEST-123 already in use")

        response = client.post(
            "/api/v1/inventory/samples",
            json={
                "compound_id": str(uuid4()),
                "quantity": 100.0,
                "barcode": "TEST-123"
            }
        )

        assert response.status_code == 400
        assert "already in use" in response.json()["detail"]

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_get_compound_sample(self, mock_get_db, mock_get_user, mock_service):
        """Test get compound sample by ID."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        sample_id = uuid4()
        mock_sample = MagicMock()
        mock_sample.id = sample_id
        mock_sample.barcode = "COMP-20250103-XYZ789"
        mock_service.get_compound_sample.return_value = mock_sample

        response = client.get(f"/api/v1/inventory/samples/{sample_id}")

        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(sample_id)

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_get_sample_not_found(self, mock_get_db, mock_get_user, mock_service):
        """Test get non-existent sample returns 404."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()
        mock_service.get_compound_sample.return_value = None

        response = client.get(f"/api/v1/inventory/samples/{uuid4()}")

        assert response.status_code == 404
        assert "not found" in response.json()["detail"]

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_list_compound_samples(self, mock_get_db, mock_get_user, mock_service):
        """Test list compound samples with filters."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        mock_samples = [MagicMock(), MagicMock()]
        mock_service.list_compound_samples.return_value = mock_samples

        response = client.get("/api/v1/inventory/samples?status=available&limit=50")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        mock_service.list_compound_samples.assert_called_once()

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_dispense_from_sample(self, mock_get_db, mock_get_user, mock_service):
        """Test dispensing quantity from sample."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        sample_id = uuid4()
        mock_sample = MagicMock()
        mock_sample.id = sample_id
        mock_sample.quantity = 100.0
        mock_service.get_compound_sample.return_value = mock_sample

        mock_updated = MagicMock()
        mock_updated.quantity = 70.0
        mock_service.update_compound_sample.return_value = mock_updated

        response = client.post(
            f"/api/v1/inventory/samples/{sample_id}/dispense",
            json={"quantity": 30.0, "notes": "Dispensed for assay"}
        )

        assert response.status_code == 200
        data = response.json()
        assert data["quantity"] == 70.0

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_dispense_insufficient_quantity(self, mock_get_db, mock_get_user, mock_service):
        """Test dispensing more than available quantity fails."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        sample_id = uuid4()
        mock_sample = MagicMock()
        mock_sample.quantity = 20.0  # Less than requested
        mock_service.get_compound_sample.return_value = mock_sample

        response = client.post(
            f"/api/v1/inventory/samples/{sample_id}/dispense",
            json={"quantity": 50.0}
        )

        assert response.status_code == 400
        assert "Insufficient quantity" in response.json()["detail"]


class TestCompoundPlateEndpoints:
    """Tests for compound plate API endpoints."""

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_create_plate_success(self, mock_get_db, mock_get_user, mock_service):
        """Test successful plate creation."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        mock_get_db.return_value = MagicMock()

        mock_plate = MagicMock()
        mock_plate.id = uuid4()
        mock_plate.barcode = "PLATE-20250103-DEF456"
        mock_plate.plate_format = "384"
        mock_plate.plate_type = "mother"
        mock_service.create_compound_plate.return_value = mock_plate

        response = client.post(
            "/api/v1/inventory/plates",
            json={
                "barcode": "PLATE-20250103-DEF456",
                "plate_format": "384",
                "plate_type": "mother"
            }
        )

        assert response.status_code == 201
        data = response.json()
        assert data["barcode"] == "PLATE-20250103-DEF456"
        assert data["plate_format"] == "384"

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_get_plate_with_contents(self, mock_get_db, mock_get_user, mock_service):
        """Test get plate with well contents."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        plate_id = uuid4()
        mock_plate = MagicMock()
        mock_plate.id = plate_id
        mock_service.get_plate.return_value = mock_plate

        mock_contents = [MagicMock(), MagicMock()]
        mock_service.get_plate_contents.return_value = mock_contents

        response = client.get(f"/api/v1/inventory/plates/{plate_id}")

        assert response.status_code == 200
        data = response.json()
        assert "plate" in data
        assert "wells" in data
        assert len(data["wells"]) == 2


class TestRequestWorkflowEndpoints:
    """Tests for compound request workflow API endpoints."""

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_create_request_success(self, mock_get_db, mock_get_user, mock_service):
        """Test successful request creation."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        mock_get_db.return_value = MagicMock()

        mock_request = MagicMock()
        mock_request.id = uuid4()
        mock_request.status = "requested"
        mock_request.requested_quantity = 50.0
        mock_service.create_request.return_value = mock_request

        response = client.post(
            "/api/v1/inventory/requests",
            json={
                "compound_id": str(uuid4()),
                "requested_quantity": 50.0,
                "priority": "high",
                "purpose": "screening"
            }
        )

        assert response.status_code == 201
        data = response.json()
        assert data["status"] == "requested"
        assert data["requested_quantity"] == 50.0

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_approve_request(self, mock_get_db, mock_get_user, mock_service):
        """Test request approval."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        mock_get_db.return_value = MagicMock()

        request_id = uuid4()
        mock_request = MagicMock()
        mock_request.status = "approved"
        mock_service.approve_request.return_value = mock_request

        response = client.put(f"/api/v1/inventory/requests/{request_id}/approve")

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "approved"

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_fulfill_request(self, mock_get_db, mock_get_user, mock_service):
        """Test request fulfillment."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        mock_get_db.return_value = MagicMock()

        request_id = uuid4()
        mock_request = MagicMock()
        mock_request.status = "fulfilled"
        mock_service.fulfill_request.return_value = mock_request

        response = client.put(
            f"/api/v1/inventory/requests/{request_id}/fulfill",
            json={"sample_id": str(uuid4())}
        )

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "fulfilled"

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_reject_request(self, mock_get_db, mock_get_user, mock_service):
        """Test request rejection."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        mock_get_db.return_value = MagicMock()

        request_id = uuid4()
        mock_request = MagicMock()
        mock_request.status = "rejected"
        mock_request.rejection_reason = "Out of stock"
        mock_service.reject_request.return_value = mock_request

        response = client.put(
            f"/api/v1/inventory/requests/{request_id}/reject",
            json={"rejection_reason": "Out of stock"}
        )

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "rejected"

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_list_pending_requests(self, mock_get_db, mock_get_user, mock_service):
        """Test list pending requests."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        mock_requests = [MagicMock(), MagicMock()]
        mock_service.list_pending_requests.return_value = mock_requests

        response = client.get("/api/v1/inventory/requests?pending_only=true")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2


class TestBarcodeLookupEndpoint:
    """Tests for barcode lookup API endpoint."""

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_lookup_sample_barcode(self, mock_get_db, mock_get_user, mock_service):
        """Test barcode lookup returns sample."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        mock_sample = MagicMock()
        mock_sample.id = uuid4()
        mock_service.lookup_barcode.return_value = {"type": "sample", "entity": mock_sample}

        response = client.get("/api/v1/inventory/barcode/COMP-20250103-ABC123")

        assert response.status_code == 200
        data = response.json()
        assert data["type"] == "sample"
        assert data["sample"] is not None
        assert data["plate"] is None

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_lookup_plate_barcode(self, mock_get_db, mock_get_user, mock_service):
        """Test barcode lookup returns plate."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        mock_plate = MagicMock()
        mock_plate.id = uuid4()
        mock_service.lookup_barcode.return_value = {"type": "plate", "entity": mock_plate}

        response = client.get("/api/v1/inventory/barcode/PLATE-20250103-XYZ789")

        assert response.status_code == 200
        data = response.json()
        assert data["type"] == "plate"
        assert data["sample"] is None
        assert data["plate"] is not None

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_lookup_unknown_barcode(self, mock_get_db, mock_get_user, mock_service):
        """Test barcode lookup returns null for unknown barcode."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        mock_service.lookup_barcode.return_value = {"type": None, "entity": None}

        response = client.get("/api/v1/inventory/barcode/UNKNOWN-BARCODE-123")

        assert response.status_code == 200
        data = response.json()
        assert data["type"] is None
        assert data["sample"] is None
        assert data["plate"] is None


class TestAlertEndpoints:
    """Tests for inventory alert API endpoints."""

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_low_stock_alerts(self, mock_get_db, mock_get_user, mock_service):
        """Test low stock alerts endpoint."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        mock_samples = [MagicMock(), MagicMock()]
        mock_service.get_low_stock_samples.return_value = mock_samples

        response = client.get("/api/v1/inventory/alerts/low-stock?threshold=50.0")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        mock_service.get_low_stock_samples.assert_called_once()

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_expiring_alerts(self, mock_get_db, mock_get_user, mock_service):
        """Test expiring samples alerts endpoint."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        mock_samples = [MagicMock()]
        mock_service.get_expiring_samples.return_value = mock_samples

        response = client.get("/api/v1/inventory/alerts/expiring?days=15")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1

    @patch("amprenta_rag.api.routers.inventory.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_expired_alerts(self, mock_get_db, mock_get_user, mock_service):
        """Test expired samples alerts endpoint."""
        mock_get_user.return_value = MagicMock()
        mock_get_db.return_value = MagicMock()

        mock_samples = []
        mock_service.get_expired_samples.return_value = mock_samples

        response = client.get("/api/v1/inventory/alerts/expired")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 0


class TestInventoryAPIIntegration:
    """Integration tests for inventory API structure."""

    def test_inventory_router_registered(self):
        """Test that inventory router is registered in main app."""
        # Check if inventory endpoints are available
        from amprenta_rag.api.main import app
        routes = [route.path for route in app.routes]
        
        # Should have inventory endpoints
        inventory_routes = [r for r in routes if "/inventory" in r]
        assert len(inventory_routes) > 0

    def test_endpoint_count(self):
        """Test that all expected inventory endpoints exist."""
        from amprenta_rag.api.main import app
        
        inventory_endpoints = []
        for route in app.routes:
            if hasattr(route, 'path') and '/inventory' in route.path:
                if hasattr(route, 'methods'):
                    for method in route.methods:
                        if method != 'HEAD':  # Skip HEAD methods
                            inventory_endpoints.append(f"{method} {route.path}")
        
        # Should have at least 15 endpoints (samples, plates, requests, alerts, barcode)
        assert len(inventory_endpoints) >= 15

    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_authentication_required(self, mock_get_db, mock_get_user):
        """Test that inventory endpoints require authentication."""
        # Mock authentication failure
        from fastapi import HTTPException
        mock_get_user.side_effect = HTTPException(status_code=401, detail="Not authenticated")
        
        response = client.get("/api/v1/inventory/samples")
        assert response.status_code == 401
