"""Integration tests for GEO sync adapter registration and API endpoint."""

from __future__ import annotations

import os
import pytest
from unittest.mock import MagicMock, patch
from uuid import uuid4
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.sync.adapters.geo import GEOSyncAdapter


class TestGEOAdapterRegistration:
    """Test that GEO adapter is properly registered with the sync manager."""
    
    def test_geo_adapter_registered_in_make_manager(self):
        """Verify GEO adapter is registered when _make_manager is called."""
        from amprenta_rag.api.routers.sync import _make_manager
        
        with patch("amprenta_rag.api.routers.sync.GEORepository") as mock_repo_class:
            mock_repo = MagicMock()
            mock_repo_class.return_value = mock_repo
            
            mgr = _make_manager()
            
            # Verify GEO adapter is registered
            assert "geo" in mgr.adapters
            assert isinstance(mgr.adapters["geo"], GEOSyncAdapter)
            
            # Verify all three sources are registered
            assert "chembl" in mgr.adapters
            assert "pubchem" in mgr.adapters
            assert "geo" in mgr.adapters
    
    def test_geo_adapter_uses_repository(self):
        """Verify GEO adapter composes with GEORepository."""
        from amprenta_rag.api.routers.sync import _make_manager
        
        with patch("amprenta_rag.api.routers.sync.GEORepository") as mock_repo_class:
            mock_repo = MagicMock()
            mock_repo_class.return_value = mock_repo
            
            mgr = _make_manager()
            
            # Verify repository was created with correct params
            mock_repo_class.assert_called_once()
            
            # Verify adapter has the repository
            geo_adapter = mgr.adapters["geo"]
            assert geo_adapter._repo == mock_repo


class TestGEOSyncAPIEndpoint:
    """Test the sync API endpoint with GEO source."""
    
    def setup_method(self):
        """Set up test environment."""
        self.client = TestClient(app)
        
    def teardown_method(self):
        """Clean up after each test."""
        app.dependency_overrides.clear()
        # Reset env vars
        if "USE_CELERY" in os.environ:
            del os.environ["USE_CELERY"]
    
    def _create_mock_user(self, role: str = "researcher") -> MagicMock:
        """Create a mock user."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_user.role = role
        mock_user.username = f"test_{role}"
        mock_user.email = f"test_{role}@example.com"
        return mock_user
    
    def _override_auth(self, user: MagicMock) -> None:
        """Override the authentication dependency."""
        app.dependency_overrides[get_current_user] = lambda: user
    
    @patch("amprenta_rag.api.routers.sync.db_session")
    @patch("amprenta_rag.api.routers.sync.GEORepository")
    def test_sync_geo_creates_job(self, mock_repo_class, mock_db_session):
        """Test that POST /sync/run with source=geo creates a sync job."""
        user = self._create_mock_user()
        self._override_auth(user)
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__ = MagicMock(return_value=mock_db)
        mock_db_session.return_value.__exit__ = MagicMock(return_value=False)
        
        # Mock the job creation
        mock_job = MagicMock()
        mock_job.id = uuid4()
        mock_job.status = "pending"
        mock_db.add = MagicMock()
        mock_db.commit = MagicMock()
        mock_db.refresh = MagicMock(side_effect=lambda j: setattr(j, 'id', mock_job.id))
        
        # Set env var and mock Celery task
        os.environ["USE_CELERY"] = "true"
        with patch("amprenta_rag.jobs.tasks.sync.run_sync_job") as mock_task:
            mock_task.delay = MagicMock()
            
            response = self.client.post(
                "/api/sync/run",
                json={"source": "geo", "sync_type": "incremental"}
            )
        
        # Should create job successfully
        assert response.status_code == 200
        data = response.json()
        assert "job_id" in data
        assert data["status"] == "pending"
    
    @patch("amprenta_rag.api.routers.sync.db_session")
    @patch("amprenta_rag.api.routers.sync.GEORepository")
    def test_sync_geo_returns_job_id(self, mock_repo_class, mock_db_session):
        """Test that GEO sync response contains a valid job ID."""
        user = self._create_mock_user()
        self._override_auth(user)
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__ = MagicMock(return_value=mock_db)
        mock_db_session.return_value.__exit__ = MagicMock(return_value=False)
        
        # Mock the job with a specific ID
        expected_job_id = uuid4()
        mock_db.refresh = MagicMock(side_effect=lambda j: setattr(j, 'id', expected_job_id))
        
        os.environ["USE_CELERY"] = "false"
        with patch("amprenta_rag.api.routers.sync.asyncio.create_task"):
            response = self.client.post(
                "/api/sync/run",
                json={"source": "geo"}
            )
        
        assert response.status_code == 200
        data = response.json()
        assert data["job_id"] == str(expected_job_id)
    
    def test_sync_invalid_source_returns_422(self):
        """Test that empty source returns 422 validation error."""
        user = self._create_mock_user()
        self._override_auth(user)
        
        # Empty source should fail Pydantic validation (min_length=1)
        response = self.client.post(
            "/api/sync/run",
            json={"source": "", "sync_type": "incremental"}
        )
        
        # Pydantic validation should reject empty string
        assert response.status_code == 422
    
    @patch("amprenta_rag.api.routers.sync.db_session")
    @patch("amprenta_rag.api.routers.sync.GEORepository")
    def test_sync_geo_triggers_celery_task(self, mock_repo_class, mock_db_session):
        """Test that GEO sync triggers a Celery task when USE_CELERY=true."""
        user = self._create_mock_user()
        self._override_auth(user)
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__ = MagicMock(return_value=mock_db)
        mock_db_session.return_value.__exit__ = MagicMock(return_value=False)
        
        job_id = uuid4()
        mock_db.refresh = MagicMock(side_effect=lambda j: setattr(j, 'id', job_id))
        
        os.environ["USE_CELERY"] = "true"
        with patch("amprenta_rag.jobs.tasks.sync.run_sync_job") as mock_task:
            mock_task.delay = MagicMock()
            
            response = self.client.post(
                "/api/sync/run",
                json={"source": "geo", "sync_type": "incremental"}
            )
        
        assert response.status_code == 200
        # Verify Celery task was called with job ID
        mock_task.delay.assert_called_once_with(str(job_id))
    
    @patch("amprenta_rag.api.routers.sync.db_session")
    @patch("amprenta_rag.api.routers.sync.GEORepository")
    def test_sync_geo_with_full_sync_type(self, mock_repo_class, mock_db_session):
        """Test GEO sync with sync_type=full."""
        user = self._create_mock_user()
        self._override_auth(user)
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__ = MagicMock(return_value=mock_db)
        mock_db_session.return_value.__exit__ = MagicMock(return_value=False)
        
        job_id = uuid4()
        mock_db.refresh = MagicMock(side_effect=lambda j: setattr(j, 'id', job_id))
        
        os.environ["USE_CELERY"] = "false"
        with patch("amprenta_rag.api.routers.sync.asyncio.create_task"):
            response = self.client.post(
                "/api/sync/run",
                json={"source": "geo", "sync_type": "full"}
            )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "pending"


class TestGEOJobListing:
    """Test job listing with GEO source filter."""
    
    def setup_method(self):
        """Set up test environment."""
        self.client = TestClient(app)
        
    def teardown_method(self):
        """Clean up after each test."""
        app.dependency_overrides.clear()
    
    def _create_mock_user(self) -> MagicMock:
        """Create a mock user."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_user.role = "researcher"
        return mock_user
    
    def _override_auth(self, user: MagicMock) -> None:
        """Override the authentication dependency."""
        app.dependency_overrides[get_current_user] = lambda: user
    
    @patch("amprenta_rag.api.routers.sync.db_session")
    def test_list_jobs_filter_by_geo_source(self, mock_db_session):
        """Test listing jobs filtered by GEO source."""
        user = self._create_mock_user()
        self._override_auth(user)
        
        # Mock database session and query
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__ = MagicMock(return_value=mock_db)
        mock_db_session.return_value.__exit__ = MagicMock(return_value=False)
        
        # Mock query chain
        mock_query = MagicMock()
        mock_db.query.return_value = mock_query
        mock_query.order_by.return_value = mock_query
        mock_query.filter.return_value = mock_query
        mock_query.count.return_value = 0
        mock_query.offset.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = []
        
        response = self.client.get("/api/sync/jobs?source=geo")
        
        assert response.status_code == 200
        data = response.json()
        assert "jobs" in data
        assert data["total"] == 0
