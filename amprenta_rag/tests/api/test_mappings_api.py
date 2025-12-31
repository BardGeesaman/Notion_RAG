"""Tests for ID mapping API endpoints."""

import pytest
from uuid import uuid4
from unittest.mock import MagicMock, patch

from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user, get_database_session


class TestMappingsAPI:
    """Test cases for ID mapping API endpoints."""

    @pytest.fixture
    def mock_user(self):
        """Mock authenticated user."""
        user = MagicMock()
        user.id = uuid4()
        user.role = "user"
        return user

    @pytest.fixture
    def mock_admin_user(self):
        """Mock authenticated admin user."""
        user = MagicMock()
        user.id = uuid4()
        user.role = "admin"
        return user

    @pytest.fixture
    def client(self, mock_user):
        """Test client with mocked authentication."""
        def override_get_current_user():
            return mock_user

        def override_get_database_session():
            return MagicMock()

        app.dependency_overrides[get_current_user] = override_get_current_user
        app.dependency_overrides[get_database_session] = override_get_database_session
        
        try:
            yield TestClient(app)
        finally:
            app.dependency_overrides.clear()

    @pytest.fixture
    def admin_client(self, mock_admin_user):
        """Test client with mocked admin authentication."""
        def override_get_current_user():
            return mock_admin_user

        def override_get_database_session():
            return MagicMock()

        app.dependency_overrides[get_current_user] = override_get_current_user
        app.dependency_overrides[get_database_session] = override_get_database_session
        
        try:
            yield TestClient(app)
        finally:
            app.dependency_overrides.clear()

    def test_refresh_endpoint_admin_only(self, client):
        """Test that refresh endpoint rejects non-admin users."""
        response = client.post(
            "/api/v1/mappings/refresh",
            json={"source": "uniprot"}
        )
        
        assert response.status_code == 403
        assert "Admin access required" in response.json()["detail"]

    @patch('amprenta_rag.api.routers.mappings.refresh_uniprot_mappings_task')
    def test_refresh_endpoint_queues_task(self, mock_task, admin_client):
        """Test that admin can trigger refresh job."""
        # Mock Celery task
        mock_task_instance = MagicMock()
        mock_task_instance.id = "test-task-123"
        mock_task.delay.return_value = mock_task_instance
        
        response = admin_client.post(
            "/api/v1/mappings/refresh",
            json={"source": "uniprot"}
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "queued"
        assert data["job_id"] == "test-task-123"
        assert "UniProt refresh job queued" in data["message"]
        mock_task.delay.assert_called_once()

    @patch('amprenta_rag.api.routers.mappings.get_mapping_stats')
    def test_status_endpoint_returns_stats(self, mock_stats, client):
        """Test status endpoint returns mapping statistics."""
        mock_stats.return_value = {
            "total_mappings": 50000,
            "by_source_type": {"gene": 25000, "protein": 25000},
            "expired_mappings": 100
        }
        
        response = client.get("/api/v1/mappings/status")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total_mappings"] == 50000
        assert data["mappings_by_type"]["gene"] == 25000
        assert data["expired_count"] == 100

    @patch('amprenta_rag.api.routers.mappings.get_mapping_stats')
    def test_stats_endpoint_returns_metrics(self, mock_stats, client):
        """Test detailed stats endpoint returns comprehensive metrics."""
        mock_stats.return_value = {
            "total_mappings": 75000,
            "by_source_type": {"gene": 30000, "protein": 25000, "metabolite": 20000},
            "by_target_type": {"uniprot": 55000, "kegg_gene": 15000, "kegg_compound": 5000},
            "expired_mappings": 150,
            "permanent_mappings": 74850
        }
        
        response = client.get("/api/v1/mappings/stats")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 75000
        assert data["by_source_type"]["gene"] == 30000
        assert data["by_target_type"]["uniprot"] == 55000
        assert data["expired"] == 150
        assert data["permanent"] == 74850

    @patch('amprenta_rag.api.routers.mappings.get_mapping')
    def test_lookup_single_found(self, mock_get_mapping, client):
        """Test single ID lookup returns mappings when found."""
        # Mock successful mappings
        def mock_mapping_side_effect(source_type, source_id, target_type, **kwargs):
            if target_type == "uniprot":
                return "P53_HUMAN"
            elif target_type == "kegg_gene":
                return "hsa:7157"
            return None
        
        mock_get_mapping.side_effect = mock_mapping_side_effect
        
        response = client.get("/api/v1/mappings/gene/TP53")
        
        assert response.status_code == 200
        data = response.json()
        assert data["source_type"] == "gene"
        assert data["source_id"] == "TP53"
        assert data["mappings"]["uniprot"] == "P53_HUMAN"
        assert data["mappings"]["kegg_gene"] == "hsa:7157"
        assert "kegg_compound" not in data["mappings"]  # Not found

    @patch('amprenta_rag.api.routers.mappings.get_mapping')
    def test_lookup_single_not_found(self, mock_get_mapping, client):
        """Test single ID lookup when no mappings found."""
        mock_get_mapping.return_value = None
        
        response = client.get("/api/v1/mappings/gene/NONEXISTENT")
        
        assert response.status_code == 200
        data = response.json()
        assert data["source_type"] == "gene"
        assert data["source_id"] == "NONEXISTENT"
        assert data["mappings"] == {}

    @patch('amprenta_rag.api.routers.mappings.get_mapping')
    def test_lookup_with_fallback_disabled(self, mock_get_mapping, client):
        """Test lookup with fallback=false parameter."""
        mock_get_mapping.return_value = "P53_HUMAN"
        
        response = client.get("/api/v1/mappings/gene/TP53?fallback=false")
        
        assert response.status_code == 200
        # Verify fallback=False was passed to the service
        mock_get_mapping.assert_called()
        # Check that fallback=False was used in at least one call
        calls = mock_get_mapping.call_args_list
        assert any(call.kwargs.get("fallback") is False for call in calls)

    @patch('amprenta_rag.api.routers.mappings.get_mappings_batch')
    def test_batch_lookup_success(self, mock_batch, client):
        """Test batch lookup with mixed results."""
        mock_batch.return_value = {
            "TP53": "hsa:7157",
            "BRCA1": "hsa:672",
            "NONEXISTENT": None
        }
        
        response = client.post(
            "/api/v1/mappings/batch",
            json={
                "ids": ["TP53", "BRCA1", "NONEXISTENT"],
                "source_type": "gene",
                "target_type": "kegg_gene",
                "organism": "human"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["results"]["TP53"] == "hsa:7157"
        assert data["results"]["BRCA1"] == "hsa:672"
        assert data["results"]["NONEXISTENT"] is None
        assert data["found"] == 2
        assert data["not_found"] == 1

    def test_batch_lookup_max_limit(self, client):
        """Test batch lookup enforces 1000 ID limit."""
        large_id_list = [f"ID_{i}" for i in range(1001)]
        
        response = client.post(
            "/api/v1/mappings/batch",
            json={
                "ids": large_id_list,
                "source_type": "gene",
                "target_type": "uniprot"
            }
        )
        
        assert response.status_code == 400
        assert "Maximum 1000 IDs per request" in response.json()["detail"]

    @patch('amprenta_rag.api.routers.mappings.get_mappings_batch')
    def test_batch_lookup_partial_results(self, mock_batch, client):
        """Test batch lookup handles partial results correctly."""
        mock_batch.return_value = {
            "GENE1": "UNIPROT1",
            "GENE2": None,
            "GENE3": "UNIPROT3",
            "GENE4": None,
            "GENE5": "UNIPROT5"
        }
        
        response = client.post(
            "/api/v1/mappings/batch",
            json={
                "ids": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
                "source_type": "gene",
                "target_type": "uniprot"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["found"] == 3  # GENE1, GENE3, GENE5
        assert data["not_found"] == 2  # GENE2, GENE4
        assert len(data["results"]) == 5
