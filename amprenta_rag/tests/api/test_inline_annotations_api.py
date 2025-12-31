"""API tests for inline annotations endpoints."""

from __future__ import annotations

import uuid
from unittest.mock import MagicMock, Mock, patch
from uuid import UUID

import pytest


# Test user fixture
TEST_USER_ID = UUID("00000000-0000-0000-0000-000000000001")


def _auth_headers():
    """Helper to create auth headers for API requests."""
    return {"X-User-Id": str(TEST_USER_ID)}


def mock_current_user():
    """Mock user for dependency override."""
    class FakeUser:
        id = TEST_USER_ID
        username = "testuser"
        email = "test@example.com"
    return FakeUser()


def test_create_annotation_success():
    """Test creating an annotation via API."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    
    # Mock database session
    mock_db = MagicMock()
    
    # Override FastAPI dependencies
    app.dependency_overrides[get_db] = lambda: mock_db
    app.dependency_overrides[get_current_user] = mock_current_user
    
    try:
        client = TestClient(app)
        
        entity_id = str(uuid.uuid4())
        payload = {
            "entity_type": "notebook",
            "entity_id": entity_id,
            "position_type": "cell",
            "position_data": {"cell_index": 5},
            "content": "This cell needs review",
            "parent_id": None,
        }
        
        with patch("amprenta_rag.api.routers.inline_annotations.create_annotation") as mock_create:
            mock_annotation = Mock()
            mock_annotation.id = uuid.uuid4()
            mock_annotation.entity_type = "notebook"
            mock_annotation.entity_id = UUID(entity_id)
            mock_annotation.position_type = "cell"
            mock_annotation.position_data = {"cell_index": 5}
            mock_annotation.content = "This cell needs review"
            mock_annotation.status = "open"
            mock_annotation.parent_id = None
            mock_annotation.created_by_id = TEST_USER_ID
            mock_annotation.resolved_by_id = None
            mock_annotation.created_at = "2025-12-31T12:00:00Z"
            mock_annotation.updated_at = None
            mock_annotation.resolved_at = None
            mock_create.return_value = mock_annotation
            
            response = client.post("/api/v1/annotations", json=payload, headers=_auth_headers())
            
            assert response.status_code == 201
            data = response.json()
            assert data["entity_type"] == "notebook"
            assert data["position_type"] == "cell"
            assert data["content"] == "This cell needs review"
            assert data["status"] == "open"
            
    finally:
        app.dependency_overrides.clear()


def test_create_annotation_invalid_position():
    """Test creating annotation with invalid position data."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    app.dependency_overrides[get_current_user] = mock_current_user
    
    try:
        client = TestClient(app)
        
        payload = {
            "entity_type": "notebook",
            "entity_id": str(uuid.uuid4()),
            "position_type": "cell",
            "position_data": {"invalid_field": 5},  # Wrong field for cell type
            "content": "This should fail",
            "parent_id": None,
        }
        
        with patch("amprenta_rag.api.routers.inline_annotations.create_annotation") as mock_create:
            mock_create.side_effect = ValueError("Invalid position_data for position_type 'cell'")
            
            response = client.post("/api/v1/annotations", json=payload, headers=_auth_headers())
            
            assert response.status_code == 400
            assert "Invalid position_data" in response.json()["detail"]
            
    finally:
        app.dependency_overrides.clear()


def test_list_annotations():
    """Test listing annotations for an entity."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    
    try:
        client = TestClient(app)
        
        entity_id = str(uuid.uuid4())
        
        with patch("amprenta_rag.api.routers.inline_annotations.get_annotations") as mock_get:
            with patch("amprenta_rag.api.routers.inline_annotations.get_annotation_count") as mock_count:
                mock_annotation = Mock()
                mock_annotation.id = uuid.uuid4()
                mock_annotation.entity_type = "dataset"
                mock_annotation.entity_id = UUID(entity_id)
                mock_annotation.position_type = "column"
                mock_annotation.position_data = {"column": "gene_name"}
                mock_annotation.content = "This column has issues"
                mock_annotation.status = "open"
                mock_annotation.parent_id = None
                mock_annotation.created_by_id = None
                mock_annotation.resolved_by_id = None
                mock_annotation.created_at = "2025-12-31T12:00:00Z"
                mock_annotation.updated_at = None
                mock_annotation.resolved_at = None
                
                mock_get.return_value = [mock_annotation]
                mock_count.return_value = {"open": 1, "resolved": 0, "total": 1}
                
                response = client.get(
                    f"/api/v1/annotations?entity_type=dataset&entity_id={entity_id}",
                    headers=_auth_headers()
                )
                
                assert response.status_code == 200
                data = response.json()
                assert data["total"] == 1
                assert data["open_count"] == 1
                assert data["resolved_count"] == 0
                assert len(data["annotations"]) == 1
                assert data["annotations"][0]["content"] == "This column has issues"
                
    finally:
        app.dependency_overrides.clear()


def test_list_annotations_filtered():
    """Test listing annotations with status filter."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    
    try:
        client = TestClient(app)
        
        entity_id = str(uuid.uuid4())
        
        with patch("amprenta_rag.api.routers.inline_annotations.get_annotations") as mock_get:
            with patch("amprenta_rag.api.routers.inline_annotations.get_annotation_count") as mock_count:
                mock_get.return_value = []  # No open annotations
                mock_count.return_value = {"open": 0, "resolved": 2, "total": 2}
                
                response = client.get(
                    f"/api/v1/annotations?entity_type=dataset&entity_id={entity_id}&status=open",
                    headers=_auth_headers()
                )
                
                assert response.status_code == 200
                data = response.json()
                assert data["total"] == 2  # Total for entity
                assert data["open_count"] == 0
                assert len(data["annotations"]) == 0  # Filtered results
                
                # Verify filter was passed to service
                mock_get.assert_called_once_with(
                    entity_type="dataset",
                    entity_id=UUID(entity_id),
                    db=mock_db,
                    status="open",
                    position_type=None,
                )
                
    finally:
        app.dependency_overrides.clear()


def test_get_annotation():
    """Test getting a single annotation by ID."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    
    try:
        client = TestClient(app)
        
        annotation_id = str(uuid.uuid4())
        
        with patch("amprenta_rag.api.routers.inline_annotations.get_annotation") as mock_get:
            mock_annotation = Mock()
            mock_annotation.id = UUID(annotation_id)
            mock_annotation.entity_type = "experiment"
            mock_annotation.entity_id = uuid.uuid4()
            mock_annotation.position_type = "field"
            mock_annotation.position_data = {"field": "protocol_version"}
            mock_annotation.content = "Protocol needs updating"
            mock_annotation.status = "open"
            mock_annotation.parent_id = None
            mock_annotation.created_by_id = None
            mock_annotation.resolved_by_id = None
            mock_annotation.created_at = "2025-12-31T12:00:00Z"
            mock_annotation.updated_at = None
            mock_annotation.resolved_at = None
            mock_get.return_value = mock_annotation
            
            response = client.get(f"/api/v1/annotations/{annotation_id}", headers=_auth_headers())
            
            assert response.status_code == 200
            data = response.json()
            assert data["id"] == annotation_id
            assert data["content"] == "Protocol needs updating"
            assert data["position_type"] == "field"
            
    finally:
        app.dependency_overrides.clear()


def test_get_annotation_not_found():
    """Test getting a non-existent annotation."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    
    try:
        client = TestClient(app)
        
        annotation_id = str(uuid.uuid4())
        
        with patch("amprenta_rag.api.routers.inline_annotations.get_annotation") as mock_get:
            mock_get.return_value = None
            
            response = client.get(f"/api/v1/annotations/{annotation_id}", headers=_auth_headers())
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
            
    finally:
        app.dependency_overrides.clear()


def test_resolve_annotation():
    """Test resolving an annotation."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    app.dependency_overrides[get_current_user] = mock_current_user
    
    try:
        client = TestClient(app)
        
        annotation_id = str(uuid.uuid4())
        
        with patch("amprenta_rag.api.routers.inline_annotations.resolve_annotation") as mock_resolve:
            mock_annotation = Mock()
            mock_annotation.id = UUID(annotation_id)
            mock_annotation.entity_type = "dataset"
            mock_annotation.entity_id = uuid.uuid4()
            mock_annotation.position_type = "row"
            mock_annotation.position_data = {"row_index": 42}
            mock_annotation.content = "Row has outliers"
            mock_annotation.status = "resolved"
            mock_annotation.parent_id = None
            mock_annotation.created_by_id = None
            mock_annotation.resolved_by_id = TEST_USER_ID
            mock_annotation.created_at = "2025-12-31T12:00:00Z"
            mock_annotation.updated_at = "2025-12-31T12:05:00Z"
            mock_annotation.resolved_at = "2025-12-31T12:05:00Z"
            mock_resolve.return_value = mock_annotation
            
            response = client.patch(f"/api/v1/annotations/{annotation_id}/resolve", headers=_auth_headers())
            
            assert response.status_code == 200
            data = response.json()
            assert data["id"] == annotation_id
            assert data["status"] == "resolved"
            assert data["resolved_by_id"] == str(TEST_USER_ID)
            
    finally:
        app.dependency_overrides.clear()


def test_reopen_annotation():
    """Test reopening a resolved annotation."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    app.dependency_overrides[get_current_user] = mock_current_user
    
    try:
        client = TestClient(app)
        
        annotation_id = str(uuid.uuid4())
        
        with patch("amprenta_rag.api.routers.inline_annotations.reopen_annotation") as mock_reopen:
            mock_annotation = Mock()
            mock_annotation.id = UUID(annotation_id)
            mock_annotation.entity_type = "notebook"
            mock_annotation.entity_id = uuid.uuid4()
            mock_annotation.position_type = "cell"
            mock_annotation.position_data = {"cell_index": 8}
            mock_annotation.content = "Cell needs more work"
            mock_annotation.status = "open"
            mock_annotation.parent_id = None
            mock_annotation.created_by_id = None
            mock_annotation.resolved_by_id = None
            mock_annotation.created_at = "2025-12-31T12:00:00Z"
            mock_annotation.updated_at = "2025-12-31T12:10:00Z"
            mock_annotation.resolved_at = None
            mock_reopen.return_value = mock_annotation
            
            response = client.patch(f"/api/v1/annotations/{annotation_id}/reopen", headers=_auth_headers())
            
            assert response.status_code == 200
            data = response.json()
            assert data["id"] == annotation_id
            assert data["status"] == "open"
            assert data["resolved_by_id"] is None
            assert data["resolved_at"] is None
            
    finally:
        app.dependency_overrides.clear()


def test_reply_to_annotation():
    """Test adding a reply to an annotation."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    app.dependency_overrides[get_current_user] = mock_current_user
    
    try:
        client = TestClient(app)
        
        parent_id = str(uuid.uuid4())
        reply_id = str(uuid.uuid4())
        payload = {"content": "I've fixed this issue"}
        
        with patch("amprenta_rag.api.routers.inline_annotations.reply_to_annotation") as mock_reply:
            mock_reply_annotation = Mock()
            mock_reply_annotation.id = UUID(reply_id)
            mock_reply_annotation.entity_type = "dataset"
            mock_reply_annotation.entity_id = uuid.uuid4()
            mock_reply_annotation.position_type = "column"
            mock_reply_annotation.position_data = {"column": "gene_name"}
            mock_reply_annotation.content = "I've fixed this issue"
            mock_reply_annotation.status = "open"
            mock_reply_annotation.parent_id = UUID(parent_id)
            mock_reply_annotation.created_by_id = TEST_USER_ID
            mock_reply_annotation.resolved_by_id = None
            mock_reply_annotation.created_at = "2025-12-31T12:15:00Z"
            mock_reply_annotation.updated_at = None
            mock_reply_annotation.resolved_at = None
            mock_reply.return_value = mock_reply_annotation
            
            response = client.post(f"/api/v1/annotations/{parent_id}/reply", json=payload, headers=_auth_headers())
            
            assert response.status_code == 201
            data = response.json()
            assert data["id"] == reply_id
            assert data["content"] == "I've fixed this issue"
            assert data["parent_id"] == parent_id
            assert data["created_by_id"] == str(TEST_USER_ID)
            
    finally:
        app.dependency_overrides.clear()


def test_delete_annotation_success():
    """Test deleting an annotation successfully."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    app.dependency_overrides[get_current_user] = mock_current_user
    
    try:
        client = TestClient(app)
        
        annotation_id = str(uuid.uuid4())
        
        with patch("amprenta_rag.api.routers.inline_annotations.delete_annotation") as mock_delete:
            mock_delete.return_value = True  # Successfully deleted
            
            response = client.delete(f"/api/v1/annotations/{annotation_id}", headers=_auth_headers())
            
            assert response.status_code == 204
            
    finally:
        app.dependency_overrides.clear()


def test_delete_annotation_unauthorized():
    """Test deleting an annotation without authorization."""
    pytest.importorskip("fastapi")
    
    from fastapi.testclient import TestClient
    from amprenta_rag.api.main import app
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    
    mock_db = MagicMock()
    app.dependency_overrides[get_db] = lambda: mock_db
    app.dependency_overrides[get_current_user] = mock_current_user
    
    try:
        client = TestClient(app)
        
        annotation_id = str(uuid.uuid4())
        
        with patch("amprenta_rag.api.routers.inline_annotations.delete_annotation") as mock_delete:
            mock_delete.return_value = False  # Not authorized or not found
            
            response = client.delete(f"/api/v1/annotations/{annotation_id}", headers=_auth_headers())
            
            assert response.status_code == 403
            assert "not authorized" in response.json()["detail"].lower()
            
    finally:
        app.dependency_overrides.clear()
