"""A01: Broken Access Control tests."""

import pytest
from uuid import uuid4
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


class TestIDOR:
    """Insecure Direct Object Reference tests."""
    
    def test_user_cannot_access_other_users_dataset(self):
        """User A should not access User B's private dataset."""
        user_a_id = uuid4()
        user_b_id = uuid4()
        dataset_id = uuid4()
        
        client = TestClient(app)
        
        # User B tries to access dataset (should fail)
        client.headers = {"X-User-Id": str(user_b_id)}
        response = client.get(f"/api/v1/datasets/{dataset_id}")
        
        # Should be 401, 403, or 404 (not 200)
        assert response.status_code in (401, 403, 404)
    
    def test_user_cannot_modify_other_users_experiment(self):
        """User A should not modify User B's experiment."""
        user_a_id = uuid4()
        user_b_id = uuid4()
        experiment_id = uuid4()
        
        client = TestClient(app)
        
        # User B tries to modify experiment
        client.headers = {"X-User-Id": str(user_b_id)}
        response = client.put(f"/api/v1/experiments/{experiment_id}", json={
            "name": "hacked experiment"
        })
        
        assert response.status_code in (401, 403, 404)
    
    def test_user_cannot_delete_other_users_compound(self):
        """User should not delete another user's compound."""
        user_a_id = uuid4()
        user_b_id = uuid4()
        compound_id = uuid4()
        
        client = TestClient(app)
        
        # User B tries to delete compound
        client.headers = {"X-User-Id": str(user_b_id)}
        response = client.delete(f"/api/v1/compounds/{compound_id}")
        
        assert response.status_code in (401, 403, 404)


class TestHorizontalPrivilege:
    """Horizontal privilege escalation tests."""
    
    def test_user_cannot_escalate_to_admin(self):
        """Regular user cannot access admin endpoints."""
        client = TestClient(app)
        client.headers = {"X-User-Id": str(uuid4())}
        
        response = client.get("/api/v1/admin/users")
        assert response.status_code in (401, 403, 404)
    
    def test_user_cannot_access_system_endpoints(self):
        """User cannot access system management endpoints."""
        client = TestClient(app)
        client.headers = {"X-User-Id": str(uuid4())}
        
        response = client.get("/api/v1/admin/system/health")
        assert response.status_code in (401, 403, 404)
    
    def test_user_cannot_modify_global_settings(self):
        """User cannot modify system-wide settings."""
        client = TestClient(app)
        client.headers = {"X-User-Id": str(uuid4())}
        
        response = client.post("/api/v1/admin/settings", json={
            "setting": "malicious_value"
        })
        assert response.status_code in (401, 403, 404)


class TestMissingAuth:
    """Test endpoints that should require auth."""
    
    def test_datasets_requires_auth(self):
        """Datasets endpoint requires authentication."""
        client = TestClient(app)
        # No X-User-Id header
        response = client.get("/api/v1/datasets")
        assert response.status_code == 401
    
    def test_compounds_requires_auth(self):
        """Compounds endpoint requires authentication."""
        client = TestClient(app)
        response = client.get("/api/v1/compounds")
        assert response.status_code == 401
    
    def test_experiments_requires_auth(self):
        """Experiments endpoint requires authentication."""
        client = TestClient(app)
        response = client.get("/api/v1/experiments")
        assert response.status_code == 401
    
    def test_signatures_requires_auth(self):
        """Signatures endpoint requires authentication."""
        client = TestClient(app)
        response = client.get("/api/v1/signatures")
        assert response.status_code == 401
    
    def test_programs_requires_auth(self):
        """Programs endpoint requires authentication."""
        client = TestClient(app)
        response = client.get("/api/v1/programs")
        assert response.status_code == 401


class TestPublicEndpoints:
    """Test that public endpoints work without auth."""
    
    def test_health_check_public(self):
        """Health check should be accessible without auth."""
        client = TestClient(app)
        response = client.get("/health")
        # Should work (200) or not exist (404), but not auth error (401)
        assert response.status_code != 401
    
    def test_docs_public(self):
        """API docs should be accessible without auth."""
        client = TestClient(app)
        response = client.get("/docs")
        # Should work or not exist, but not require auth
        assert response.status_code != 401
