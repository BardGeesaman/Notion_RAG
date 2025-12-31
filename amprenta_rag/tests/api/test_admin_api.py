"""Tests for admin API endpoints."""

import pytest
from unittest.mock import MagicMock, patch
from fastapi.testclient import TestClient
from uuid import uuid4

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user


class TestAdminAPI:
    """Test admin API endpoints."""
    
    def setup_method(self):
        """Set up test environment."""
        self.client = TestClient(app)
    
    def teardown_method(self):
        """Clean up after each test."""
        # Clear any dependency overrides
        app.dependency_overrides.clear()
    
    def _create_mock_user(self, role: str = "admin") -> MagicMock:
        """Create a mock user with specified role."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_user.role = role
        mock_user.username = f"test_{role}"
        mock_user.email = f"test_{role}@example.com"
        return mock_user
    
    def _override_auth(self, user: MagicMock) -> None:
        """Override the authentication dependency."""
        app.dependency_overrides[get_current_user] = lambda: user
    
    def test_list_caches_admin_only(self):
        """Test that cache listing requires admin role."""
        # Test with non-admin user
        non_admin_user = self._create_mock_user(role="researcher")
        self._override_auth(non_admin_user)
        
        response = self.client.get("/api/v1/admin/caches")
        assert response.status_code == 403
        assert "Admin access required" in response.json()["detail"]
        
        # Test with no user (None)
        app.dependency_overrides[get_current_user] = lambda: None
        response = self.client.get("/api/v1/admin/caches")
        assert response.status_code == 403
        assert "Admin access required" in response.json()["detail"]
    
    @patch("amprenta_rag.api.routers.admin.get_cache_stats")
    def test_list_caches_returns_stats(self, mock_get_stats):
        """Test that cache listing returns proper statistics."""
        admin_user = self._create_mock_user(role="admin")
        self._override_auth(admin_user)
        
        # Mock cache stats
        mock_stats = {
            "semantic_cache": {
                "type": "SemanticCache",
                "entries": 5,
                "ttl_seconds": 3600,
                "similarity_threshold": 0.92
            },
            "dataset_feature_cache": {
                "type": "DatasetFeatureCache",
                "entries": 10,
                "ttl_seconds": 3600,
                "max_size": 1000,
                "hits": 25,
                "misses": 5,
                "evictions": 0
            }
        }
        mock_get_stats.return_value = mock_stats
        
        response = self.client.get("/api/v1/admin/caches")
        assert response.status_code == 200
        
        data = response.json()
        assert "semantic_cache" in data
        assert "dataset_feature_cache" in data
        assert data["semantic_cache"]["entries"] == 5
        assert data["dataset_feature_cache"]["hits"] == 25
        
        mock_get_stats.assert_called_once()
    
    @patch("amprenta_rag.api.routers.admin.clear_cache")
    def test_clear_cache_success(self, mock_clear_cache):
        """Test successful cache clearing."""
        admin_user = self._create_mock_user(role="admin")
        self._override_auth(admin_user)
        mock_clear_cache.return_value = True
        
        response = self.client.post("/api/v1/admin/caches/semantic_cache/clear")
        assert response.status_code == 200
        
        data = response.json()
        assert data["status"] == "cleared"
        assert data["cache"] == "semantic_cache"
        assert "Successfully cleared" in data["message"]
        
        mock_clear_cache.assert_called_once_with("semantic_cache")
    
    @patch("amprenta_rag.api.routers.admin.clear_cache")
    def test_clear_cache_not_found(self, mock_clear_cache):
        """Test clearing non-existent cache returns 404."""
        admin_user = self._create_mock_user(role="admin")
        self._override_auth(admin_user)
        mock_clear_cache.return_value = False
        
        response = self.client.post("/api/v1/admin/caches/invalid_cache/clear")
        assert response.status_code == 404
        assert "not found or failed to clear" in response.json()["detail"]
        
        mock_clear_cache.assert_called_once_with("invalid_cache")
    
    @patch("amprenta_rag.api.routers.admin.clear_all_caches")
    def test_clear_all_caches(self, mock_clear_all):
        """Test clearing all caches."""
        admin_user = self._create_mock_user(role="admin")
        self._override_auth(admin_user)
        
        # Mock results with some successes and failures
        mock_results = {
            "semantic_cache": True,
            "dataset_feature_cache": True,
            "enhanced_feature_cache": False
        }
        mock_clear_all.return_value = mock_results
        
        response = self.client.post("/api/v1/admin/caches/clear-all")
        assert response.status_code == 200
        
        data = response.json()
        assert data["status"] == "completed"
        assert data["results"] == mock_results
        assert data["summary"]["total_caches"] == 3
        assert data["summary"]["successful"] == 2
        assert data["summary"]["failed"] == 1
        
        mock_clear_all.assert_called_once()
    
    def test_non_admin_forbidden(self):
        """Test that all admin endpoints reject non-admin users."""
        researcher_user = self._create_mock_user(role="researcher")
        self._override_auth(researcher_user)
        
        endpoints = [
            ("GET", "/api/v1/admin/caches"),
            ("GET", "/api/v1/admin/caches/summary"),
            ("POST", "/api/v1/admin/caches/semantic_cache/clear"),
            ("POST", "/api/v1/admin/caches/clear-all"),
        ]
        
        for method, endpoint in endpoints:
            if method == "GET":
                response = self.client.get(endpoint)
            else:  # POST
                response = self.client.post(endpoint)
            
            assert response.status_code == 403, f"Endpoint {method} {endpoint} should require admin"
            assert "Admin access required" in response.json()["detail"]
    
    @patch("amprenta_rag.api.routers.admin.get_cache_summary")
    def test_cache_summary_endpoint(self, mock_get_summary):
        """Test the cache summary endpoint."""
        admin_user = self._create_mock_user(role="admin")
        self._override_auth(admin_user)
        
        mock_summary = {
            "total_caches": 3,
            "total_entries": 15,
            "caches_with_errors": [],
            "cache_names": ["semantic_cache", "dataset_feature_cache", "enhanced_feature_cache"]
        }
        mock_get_summary.return_value = mock_summary
        
        response = self.client.get("/api/v1/admin/caches/summary")
        assert response.status_code == 200
        
        data = response.json()
        assert data["total_caches"] == 3
        assert data["total_entries"] == 15
        assert len(data["cache_names"]) == 3
        
        mock_get_summary.assert_called_once()
    
    @patch("amprenta_rag.api.routers.admin.clear_cache")
    def test_clear_cache_admin_role_required(self, mock_clear_cache):
        """Test that cache clearing specifically requires admin role."""
        # Test with different non-admin roles
        roles_to_test = ["researcher", "viewer", "member", "guest"]
        
        for role in roles_to_test:
            user = self._create_mock_user(role=role)
            self._override_auth(user)
            
            response = self.client.post("/api/v1/admin/caches/semantic_cache/clear")
            assert response.status_code == 403, f"Role '{role}' should be forbidden"
            assert "Admin access required" in response.json()["detail"]
        
        # Ensure mock was never called due to authorization failures
        mock_clear_cache.assert_not_called()
    
    @patch("amprenta_rag.api.routers.admin.get_extended_system_info")
    def test_system_health_endpoint(self, mock_get_system_info):
        """Test the system health endpoint."""
        admin_user = self._create_mock_user(role="admin")
        self._override_auth(admin_user)
        
        mock_system_info = {
            "cpu": {
                "cpu_percent": 25.5,
                "cpu_count": 8,
                "cpu_count_logical": 16
            },
            "memory": {
                "total_gb": 16.0,
                "available_gb": 8.5,
                "used_gb": 7.5,
                "percent": 46.9
            },
            "disk": {
                "total_gb": 500.0,
                "free_gb": 250.0,
                "used_gb": 250.0,
                "percent": 50.0
            },
            "timestamp": 1640995200.0
        }
        mock_get_system_info.return_value = mock_system_info
        
        response = self.client.get("/api/v1/admin/health/system")
        assert response.status_code == 200
        
        data = response.json()
        assert "cpu" in data
        assert "memory" in data
        assert "disk" in data
        assert data["cpu"]["cpu_percent"] == 25.5
        assert data["memory"]["total_gb"] == 16.0
        assert data["disk"]["percent"] == 50.0
        
        mock_get_system_info.assert_called_once()
    
    @patch("amprenta_rag.api.routers.admin.get_celery_queue_stats")
    def test_queue_health_endpoint(self, mock_get_queue_stats):
        """Test the queue health endpoint."""
        admin_user = self._create_mock_user(role="admin")
        self._override_auth(admin_user)
        
        mock_queue_stats = {
            "active": 5,
            "reserved": 2,
            "scheduled": 10,
            "workers": ["worker1@host", "worker2@host"],
            "worker_count": 2,
            "status": "available"
        }
        mock_get_queue_stats.return_value = mock_queue_stats
        
        response = self.client.get("/api/v1/admin/health/queues")
        assert response.status_code == 200
        
        data = response.json()
        assert data["active"] == 5
        assert data["reserved"] == 2
        assert data["scheduled"] == 10
        assert data["worker_count"] == 2
        assert data["status"] == "available"
        
        mock_get_queue_stats.assert_called_once()
    
    @patch("amprenta_rag.api.routers.admin.get_connection_status")
    def test_connection_health_endpoint(self, mock_get_connection_status):
        """Test the connection health endpoint."""
        admin_user = self._create_mock_user(role="admin")
        self._override_auth(admin_user)
        
        mock_connection_status = {
            "postgresql": {
                "status": "connected",
                "type": "PostgreSQL"
            },
            "redis": {
                "status": "connected",
                "type": "Redis",
                "broker_url": "localhost:6379"
            }
        }
        mock_get_connection_status.return_value = mock_connection_status
        
        response = self.client.get("/api/v1/admin/health/connections")
        assert response.status_code == 200
        
        data = response.json()
        assert "postgresql" in data
        assert "redis" in data
        assert data["postgresql"]["status"] == "connected"
        assert data["redis"]["status"] == "connected"
        
        mock_get_connection_status.assert_called_once()
    
    def test_health_endpoints_admin_only(self):
        """Test that all health endpoints require admin access."""
        researcher_user = self._create_mock_user(role="researcher")
        self._override_auth(researcher_user)
        
        health_endpoints = [
            "/api/v1/admin/health/system",
            "/api/v1/admin/health/queues",
            "/api/v1/admin/health/connections",
        ]
        
        for endpoint in health_endpoints:
            response = self.client.get(endpoint)
            assert response.status_code == 403, f"Health endpoint {endpoint} should require admin"
            assert "Admin access required" in response.json()["detail"]
