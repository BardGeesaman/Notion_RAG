"""Tests for IP API endpoints."""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestCreateDisclosure:
    """Tests for POST /api/v1/ip/disclosures endpoint."""

    @patch("amprenta_rag.api.routers.ip.create_disclosure")
    def test_create_disclosure_success(self, mock_create):
        """Test successful disclosure creation."""
        mock_disc = MagicMock()
        mock_disc.id = uuid4()
        mock_disc.title = "Test Invention"
        mock_disc.description = "Description"
        mock_disc.status = "draft"
        mock_disc.technical_field = "Chemistry"
        mock_disc.priority_date = None
        mock_disc.created_at = datetime(2024, 1, 1)
        
        mock_create.return_value = mock_disc
        
        response = client.post(
            "/api/v1/ip/disclosures",
            json={
                "title": "Test Invention",
                "description": "Description",
                "technical_field": "Chemistry",
                "inventors": [{"user_id": str(uuid4()), "is_primary": True}],
                "user_id": str(uuid4()),
                "company_id": str(uuid4()),
            },
        )
        
        assert response.status_code == 201
        data = response.json()
        assert data["title"] == "Test Invention"


class TestListDisclosures:
    """Tests for GET /api/v1/ip/disclosures endpoint."""

    @patch("amprenta_rag.api.routers.ip.get_disclosure_portfolio")
    def test_list_with_filter(self, mock_get):
        """Test listing disclosures with status filter."""
        mock_disc = MagicMock()
        mock_disc.id = uuid4()
        mock_disc.title = "Filed Invention"
        mock_disc.description = "Description"
        mock_disc.status = "filed"
        mock_disc.technical_field = None
        mock_disc.priority_date = None
        mock_disc.created_at = datetime(2024, 1, 1)
        
        mock_get.return_value = [mock_disc]
        
        company_id = uuid4()
        response = client.get(f"/api/v1/ip/disclosures?company_id={company_id}&status=filed")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1


class TestGetDisclosure:
    """Tests for GET /api/v1/ip/disclosures/{id} endpoint."""

    def test_get_disclosure_details(self):
        """Test getting disclosure details."""
        disclosure_id = uuid4()
        
        mock_disc = MagicMock()
        mock_disc.id = disclosure_id
        mock_disc.title = "Disclosure"
        mock_disc.description = "Description"
        mock_disc.status = "draft"
        mock_disc.technical_field = None
        mock_disc.priority_date = None
        mock_disc.created_at = datetime(2024, 1, 1)
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_disc
        
        def mock_get_db():
            yield mock_db
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get(f"/api/v1/ip/disclosures/{disclosure_id}")
            
            assert response.status_code == 200
        finally:
            app.dependency_overrides.clear()


class TestUpdateStatus:
    """Tests for PUT /api/v1/ip/disclosures/{id}/status endpoint."""

    @patch("amprenta_rag.api.routers.ip.update_disclosure_status")
    def test_update_status(self, mock_update):
        """Test disclosure status update."""
        disclosure_id = uuid4()
        
        mock_disc = MagicMock()
        mock_disc.id = disclosure_id
        mock_disc.title = "Disclosure"
        mock_disc.description = "Description"
        mock_disc.status = "submitted"
        mock_disc.technical_field = None
        mock_disc.priority_date = None
        mock_disc.created_at = datetime(2024, 1, 1)
        
        mock_update.return_value = mock_disc
        
        response = client.put(
            f"/api/v1/ip/disclosures/{disclosure_id}/status",
            json={"new_status": "submitted", "user_id": str(uuid4())},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "submitted"


class TestFilePatent:
    """Tests for POST /api/v1/ip/disclosures/{id}/file-patent endpoint."""

    @patch("amprenta_rag.api.routers.ip.file_patent")
    def test_file_patent_success(self, mock_file):
        """Test patent filing."""
        disclosure_id = uuid4()
        
        mock_patent = MagicMock()
        mock_patent.id = uuid4()
        mock_patent.application_number = "US12345678"
        mock_patent.jurisdiction = "US"
        mock_patent.status = "pending"
        mock_patent.filing_date = datetime(2024, 1, 1)
        mock_patent.grant_date = None
        
        mock_file.return_value = mock_patent
        
        response = client.post(
            f"/api/v1/ip/disclosures/{disclosure_id}/file-patent",
            json={
                "application_number": "US12345678",
                "jurisdiction": "US",
                "filing_date": "2024-01-01T00:00:00",
            },
        )
        
        assert response.status_code == 201
        data = response.json()
        assert data["application_number"] == "US12345678"


class TestListPatents:
    """Tests for GET /api/v1/ip/patents endpoint."""

    def test_list_patents(self):
        """Test listing patents."""
        mock_patent = MagicMock()
        mock_patent.id = uuid4()
        mock_patent.application_number = "US11111111"
        mock_patent.jurisdiction = "US"
        mock_patent.status = "pending"
        mock_patent.filing_date = datetime(2024, 1, 1)
        mock_patent.grant_date = None
        
        mock_db = MagicMock()
        mock_db.query.return_value.all.return_value = [mock_patent]
        
        def mock_get_db():
            yield mock_db
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get("/api/v1/ip/patents")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
        finally:
            app.dependency_overrides.clear()


class TestGetPatentClaims:
    """Tests for GET /api/v1/ip/patents/{id}/claims endpoint."""

    def test_get_claims(self):
        """Test getting patent claims."""
        patent_id = uuid4()
        
        mock_claim = MagicMock()
        mock_claim.id = uuid4()
        mock_claim.claim_number = 1
        mock_claim.claim_text = "A compound comprising..."
        mock_claim.claim_type = "independent"
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.order_by.return_value.all.return_value = [mock_claim]
        
        def mock_get_db():
            yield mock_db
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get(f"/api/v1/ip/patents/{patent_id}/claims")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
        finally:
            app.dependency_overrides.clear()


class TestIPLinks:
    """Tests for IP link endpoints."""

    @patch("amprenta_rag.api.routers.ip.link_entity_to_disclosure")
    def test_create_link(self, mock_link):
        """Test creating IP link."""
        mock_link_obj = MagicMock()
        mock_link_obj.id = uuid4()
        mock_link_obj.entity_type = "compound"
        mock_link_obj.entity_id = uuid4()
        mock_link_obj.disclosure_id = uuid4()
        mock_link_obj.link_type = "embodiment"
        mock_link_obj.notes = None
        
        mock_link.return_value = mock_link_obj
        
        response = client.post(
            "/api/v1/ip/links",
            json={
                "entity_type": "compound",
                "entity_id": str(uuid4()),
                "disclosure_id": str(uuid4()),
                "link_type": "embodiment",
            },
        )
        
        assert response.status_code == 201

    @patch("amprenta_rag.api.routers.ip.get_entity_ip_links")
    def test_get_entity_links(self, mock_get):
        """Test getting entity IP links."""
        mock_link = MagicMock()
        mock_link.id = uuid4()
        mock_link.entity_type = "compound"
        mock_link.entity_id = uuid4()
        mock_link.disclosure_id = uuid4()
        mock_link.link_type = "embodiment"
        mock_link.notes = "Test notes"
        
        mock_get.return_value = [mock_link]
        
        entity_id = uuid4()
        response = client.get(f"/api/v1/ip/links/compound/{entity_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1

