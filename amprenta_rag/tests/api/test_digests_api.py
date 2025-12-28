"""
Unit tests for digests API endpoints.

Tests digest schedule CRUD operations with mocked dependencies.
"""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.database.models import User


# Test constants
TEST_USER_ID = uuid4()


def mock_current_user():
    """Mock current user for authentication."""
    user = User(
        id=TEST_USER_ID,
        username=f"test_user_{uuid4().hex[:8]}",
        email=f"test_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )
    # Mock the company_id attribute that may be accessed
    user.company_id = uuid4()
    return user


@pytest.fixture
def auth_client():
    """Test client with mocked authentication."""
    from amprenta_rag.api.dependencies import get_current_user_with_company
    
    app.dependency_overrides[get_current_user_with_company] = mock_current_user
    
    try:
        with TestClient(app) as test_client:
            yield test_client
    finally:
        app.dependency_overrides.clear()


class TestDigestSchedulesCreate:
    """Tests for POST /api/digests endpoint."""

    @patch("amprenta_rag.api.routers.digests.DigestScheduler")
    def test_create_digest_schedule_success(self, mock_scheduler_class, auth_client):
        """Test successful digest schedule creation."""
        program_id = uuid4()
        schedule_id = uuid4()
        
        # Mock program exists
        mock_program = MagicMock()
        mock_program.id = program_id
        mock_program.name = "Test Program"
        
        # Mock created schedule
        mock_schedule = MagicMock()
        mock_schedule.id = schedule_id
        mock_schedule.program_id = program_id
        mock_schedule.notebook_path = "notebooks/weekly_report.ipynb"
        mock_schedule.schedule_cron = "0 9 * * mon"
        mock_schedule.recipients = ["scientist@example.com"]
        mock_schedule.enabled = True
        mock_schedule.last_run_at = None
        mock_schedule.last_status = None
        
        # Mock database session
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_program
        
        def mock_add(obj):
            obj.id = schedule_id
        
        mock_session.add = mock_add
        
        def mock_refresh(obj):
            for k, v in vars(mock_schedule).items():
                if not k.startswith("_"):
                    setattr(obj, k, v)
        
        mock_session.refresh = mock_refresh
        
        def mock_get_db():
            yield mock_session
        
        # Use dependency override
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.post(
                "/api/digests",
                json={
                    "program_id": str(program_id),
                    "notebook_path": "notebooks/weekly_report.ipynb",
                    "schedule_cron": "0 9 * * mon",
                    "recipients": ["scientist@example.com"],
                    "enabled": True,
                },
            )
            
            assert response.status_code == 201
            data = response.json()
            assert data["program_id"] == str(program_id)
            assert data["notebook_path"] == "notebooks/weekly_report.ipynb"
            assert data["schedule_cron"] == "0 9 * * mon"
            assert data["recipients"] == ["scientist@example.com"]
            assert data["enabled"] is True
        finally:
            app.dependency_overrides.clear()

    def test_create_digest_schedule_program_not_found(self, auth_client):
        """Test digest schedule creation with non-existent program."""
        program_id = uuid4()
        
        # Mock database session - program not found
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.post(
                "/api/digests",
                json={
                    "program_id": str(program_id),
                    "notebook_path": "notebooks/weekly_report.ipynb",
                    "schedule_cron": "0 9 * * mon",
                    "recipients": [],
                    "enabled": True,
                },
            )
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()


class TestDigestSchedulesList:
    """Tests for GET /api/digests endpoint."""

    def test_list_digest_schedules_success(self, auth_client):
        """Test successful digest schedules list retrieval."""
        schedule1_id = uuid4()
        schedule2_id = uuid4()
        program_id = uuid4()
        
        # Mock schedules
        mock_schedule1 = MagicMock()
        mock_schedule1.id = schedule1_id
        mock_schedule1.program_id = program_id
        mock_schedule1.notebook_path = "notebooks/report1.ipynb"
        mock_schedule1.schedule_cron = "0 9 * * mon"
        mock_schedule1.recipients = ["user1@example.com"]
        mock_schedule1.enabled = True
        mock_schedule1.last_run_at = None
        mock_schedule1.last_status = None
        
        mock_schedule2 = MagicMock()
        mock_schedule2.id = schedule2_id
        mock_schedule2.program_id = program_id
        mock_schedule2.notebook_path = "notebooks/report2.ipynb"
        mock_schedule2.schedule_cron = "0 10 * * tue"
        mock_schedule2.recipients = ["user2@example.com"]
        mock_schedule2.enabled = False
        mock_schedule2.last_run_at = datetime(2024, 1, 1, 10, 0, 0)
        mock_schedule2.last_status = "success"
        
        # Mock database session
        mock_session = MagicMock()
        mock_session.query.return_value.order_by.return_value.all.return_value = [
            mock_schedule1,
            mock_schedule2,
        ]
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.get("/api/digests")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 2
            assert data[0]["id"] == str(schedule1_id)
            assert data[0]["enabled"] is True
            assert data[1]["id"] == str(schedule2_id)
            assert data[1]["enabled"] is False
        finally:
            app.dependency_overrides.clear()

    def test_list_digest_schedules_empty(self, auth_client):
        """Test listing digest schedules when none exist."""
        # Mock database session - empty result
        mock_session = MagicMock()
        mock_session.query.return_value.order_by.return_value.all.return_value = []
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.get("/api/digests")
            
            assert response.status_code == 200
            data = response.json()
            assert data == []
        finally:
            app.dependency_overrides.clear()


class TestDigestSchedulesUpdate:
    """Tests for PUT /api/digests/{schedule_id} endpoint."""

    @patch("amprenta_rag.api.routers.digests.DigestScheduler")
    def test_update_digest_schedule_success(self, mock_scheduler_class, auth_client):
        """Test successful digest schedule update."""
        schedule_id = uuid4()
        program_id = uuid4()
        
        # Mock existing schedule
        mock_schedule = MagicMock()
        mock_schedule.id = schedule_id
        mock_schedule.program_id = program_id
        mock_schedule.notebook_path = "notebooks/old_report.ipynb"
        mock_schedule.schedule_cron = "0 9 * * mon"
        mock_schedule.recipients = ["old@example.com"]
        mock_schedule.enabled = True
        mock_schedule.last_run_at = None
        mock_schedule.last_status = None
        
        # Mock database session
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_schedule
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.put(
                f"/api/digests/{schedule_id}",
                json={
                    "notebook_path": "notebooks/new_report.ipynb",
                    "enabled": False,
                },
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["id"] == str(schedule_id)
            # Note: Actual values depend on mock implementation
        finally:
            app.dependency_overrides.clear()

    def test_update_digest_schedule_not_found(self, auth_client):
        """Test updating non-existent digest schedule."""
        schedule_id = uuid4()
        
        # Mock database session - schedule not found
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.put(
                f"/api/digests/{schedule_id}",
                json={"enabled": False},
            )
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()


class TestDigestSchedulesDelete:
    """Tests for DELETE /api/digests/{schedule_id} endpoint."""

    @patch("amprenta_rag.api.routers.digests.DigestScheduler")
    def test_delete_digest_schedule_success(self, mock_scheduler_class, auth_client):
        """Test successful digest schedule deletion."""
        schedule_id = uuid4()
        
        # Mock existing schedule
        mock_schedule = MagicMock()
        mock_schedule.id = schedule_id
        
        # Mock database session
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_schedule
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.delete(f"/api/digests/{schedule_id}")
            
            assert response.status_code == 200
            data = response.json()
            assert data["deleted"] is True
        finally:
            app.dependency_overrides.clear()

    def test_delete_digest_schedule_not_found(self, auth_client):
        """Test deleting non-existent digest schedule."""
        schedule_id = uuid4()
        
        # Mock database session - schedule not found
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.delete(f"/api/digests/{schedule_id}")
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()


class TestDigestSchedulesRun:
    """Tests for POST /api/digests/{schedule_id}/run endpoint."""

    @patch("amprenta_rag.api.routers.digests.DigestScheduler")
    def test_run_digest_schedule_success(self, mock_scheduler_class, auth_client):
        """Test successful digest schedule execution."""
        schedule_id = uuid4()
        
        # Mock existing schedule
        mock_schedule = MagicMock()
        mock_schedule.id = schedule_id
        
        # Mock scheduler run result
        mock_result = MagicMock()
        mock_result.status = "success"
        mock_result.output_html = "/path/to/output.html"
        
        mock_scheduler = MagicMock()
        mock_scheduler.run_digest.return_value = mock_result
        mock_scheduler_class.return_value = mock_scheduler
        
        # Mock database session
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_schedule
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.post(f"/api/digests/{schedule_id}/run")
            
            assert response.status_code == 200
            data = response.json()
            assert data["job"] == "ran"
            assert data["status"] == "success"
            assert "output_html" in data
        finally:
            app.dependency_overrides.clear()

    def test_run_digest_schedule_not_found(self, auth_client):
        """Test running non-existent digest schedule."""
        schedule_id = uuid4()
        
        # Mock database session - schedule not found
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.post(f"/api/digests/{schedule_id}/run")
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()


class TestDigestSchedulesOutput:
    """Tests for GET /api/digests/{schedule_id}/output endpoint."""

    @patch("amprenta_rag.api.routers.digests.FileResponse")
    @patch("amprenta_rag.api.routers.digests.Path")
    def test_get_digest_output_success(self, mock_path_class, mock_file_response, auth_client):
        """Test successful digest output retrieval."""
        schedule_id = uuid4()
        
        # Mock existing schedule
        mock_schedule = MagicMock()
        mock_schedule.id = schedule_id
        
        # Mock Path to return a path that exists
        mock_html_path = MagicMock()
        mock_html_path.exists.return_value = True
        mock_html_path.__str__.return_value = "/fake/path/to/digest.html"
        
        # Create a proper mock chain for Path operations
        mock_root = MagicMock()
        mock_data_dir = MagicMock()
        mock_digests_dir = MagicMock()
        mock_schedule_dir = MagicMock()
        
        mock_root.__truediv__.return_value = mock_data_dir
        mock_data_dir.__truediv__.return_value = mock_digests_dir
        mock_digests_dir.__truediv__.return_value = mock_schedule_dir
        mock_schedule_dir.__truediv__.return_value = mock_html_path
        
        # Path(__file__).resolve().parents[3]
        mock_path_inst = MagicMock()
        mock_path_inst.resolve.return_value.parents.__getitem__.return_value = mock_root
        mock_path_class.return_value = mock_path_inst
        
        # Mock FileResponse to return a mock response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_file_response.return_value = mock_response
        
        # Mock database session
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_schedule
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.get(f"/api/digests/{schedule_id}/output")
            
            # Should return FileResponse (status 200)
            assert response.status_code == 200
        finally:
            app.dependency_overrides.clear()

    def test_get_digest_output_not_found_schedule(self, auth_client):
        """Test getting output for non-existent digest schedule."""
        schedule_id = uuid4()
        
        # Mock database session - schedule not found
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = auth_client.get(f"/api/digests/{schedule_id}/output")
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()

