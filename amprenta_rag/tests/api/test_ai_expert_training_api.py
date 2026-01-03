"""Tests for AI Expert Agents Training API endpoints."""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestTrainingExampleEndpoints:
    """Tests for training example endpoints."""

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_create_training_example(self, mock_get_db, mock_get_user, mock_service):
        """Test POST /experts/{id}/training-examples creates example."""
        mock_get_db.return_value = MagicMock()
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        expert_id = uuid4()
        mock_example = MagicMock()
        mock_example.id = uuid4()
        mock_example.expert_id = expert_id
        mock_example.question = "Test question"
        mock_example.ideal_answer = "Test answer"
        mock_example.is_approved = False
        mock_example.prompt_version = "1.0"
        mock_example.source = "manual"
        
        mock_service.create_training_example.return_value = mock_example

        response = client.post(
            f"/api/v1/experts/experts/{expert_id}/training-examples",
            json={
                "expert_id": str(expert_id),
                "input_text": "Test question",
                "ideal_output": "Test answer"
            }
        )

        assert response.status_code == 201
        data = response.json()
        assert data["question"] == "Test question"
        assert data["ideal_answer"] == "Test answer"
        assert data["is_approved"] is False

    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_list_training_examples(self, mock_get_db, mock_get_user):
        """Test GET /experts/{id}/training-examples lists examples."""
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_user.return_value = MagicMock()
        
        expert_id = uuid4()
        mock_example = MagicMock()
        mock_example.id = uuid4()
        mock_example.expert_id = expert_id
        mock_example.question = "Example question"
        mock_example.ideal_answer = "Example answer"
        mock_example.is_approved = True
        mock_example.prompt_version = "1.0"
        mock_example.source = "manual"
        
        mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = [mock_example]

        response = client.get(f"/api/v1/experts/experts/{expert_id}/training-examples")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["question"] == "Example question"

    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_list_training_examples_approved_only(self, mock_get_db, mock_get_user):
        """Test listing only approved training examples."""
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_user.return_value = MagicMock()
        
        expert_id = uuid4()
        mock_example = MagicMock()
        mock_example.is_approved = True
        
        mock_db.query.return_value.filter.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = [mock_example]

        response = client.get(f"/api/v1/experts/experts/{expert_id}/training-examples?approved_only=true")

        assert response.status_code == 200


class TestKnowledgeManagementEndpoints:
    """Tests for knowledge management endpoints."""

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_upload_knowledge_document(self, mock_get_db, mock_get_user, mock_service):
        """Test POST /experts/{id}/knowledge uploads document."""
        mock_get_db.return_value = MagicMock()
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        expert_id = uuid4()
        mock_doc = MagicMock()
        mock_doc.id = uuid4()
        mock_doc.expert_id = expert_id
        mock_doc.namespace = f"expert_{expert_id}"
        mock_doc.title = "Test Knowledge"
        mock_doc.chunk_index = 0
        mock_doc.source_type = "manual"
        mock_doc.source_url = None
        
        mock_service.add_knowledge_doc.return_value = mock_doc

        response = client.post(
            f"/api/v1/experts/experts/{expert_id}/knowledge",
            json={
                "expert_id": str(expert_id),
                "title": "Test Knowledge",
                "content": "This is test knowledge content for the expert",
                "source_type": "manual"
            }
        )

        assert response.status_code == 201
        data = response.json()
        assert data["title"] == "Test Knowledge"
        assert data["expert_id"] == str(expert_id)
        assert data["source_type"] == "manual"


class TestExpertManagementEndpoints:
    """Tests for expert management endpoints."""

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_update_expert(self, mock_get_db, mock_get_user, mock_service):
        """Test PATCH /experts/{id} updates expert."""
        mock_get_db.return_value = MagicMock()
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        expert_id = uuid4()
        mock_expert = MagicMock()
        mock_expert.id = expert_id
        mock_expert.name = "Dr. Updated"
        mock_expert.role = "Updated Expert"
        mock_expert.specializations = ["updated"]
        mock_expert.is_active = True
        mock_expert.prompt_version = "1.1"
        
        mock_service.update_expert.return_value = mock_expert

        response = client.patch(
            f"/api/v1/experts/experts/{expert_id}",
            json={
                "system_prompt": "Updated system prompt",
                "bump_version": True
            }
        )

        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(expert_id)
        assert data["prompt_version"] == "1.1"

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_export_training_data(self, mock_get_db, mock_get_user, mock_service):
        """Test POST /experts/{id}/export-training exports data."""
        mock_get_db.return_value = MagicMock()
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        expert_id = uuid4()
        mock_training_data = [
            {
                "messages": [
                    {"role": "system", "content": "You are an expert"},
                    {"role": "user", "content": "Question"},
                    {"role": "assistant", "content": "Answer"}
                ]
            }
        ]
        
        mock_service.export_training_data.return_value = mock_training_data

        response = client.post(f"/api/v1/experts/experts/{expert_id}/export-training")

        assert response.status_code == 200
        data = response.json()
        assert data["format"] == "jsonl"
        assert data["expert_id"] == str(expert_id)
        assert data["example_count"] == 1
        assert len(data["data"]) == 1


class TestFeedbackSummaryEndpoints:
    """Tests for feedback summary endpoints."""

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_get_feedback_summary(self, mock_get_db, mock_service):
        """Test GET /experts/{id}/feedback-summary returns metrics."""
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        
        expert_id = uuid4()
        mock_expert = MagicMock()
        mock_expert.name = "Dr. Feedback"
        mock_service.get_expert.return_value = mock_expert
        
        # Mock feedback data
        mock_feedback1 = MagicMock()
        mock_feedback1.rating = 4
        mock_feedback1.correction = "Good response"
        
        mock_feedback2 = MagicMock()
        mock_feedback2.rating = 5
        mock_feedback2.correction = None
        
        mock_db.query.return_value.join.return_value.filter.return_value.all.return_value = [mock_feedback1, mock_feedback2]

        response = client.get(f"/api/v1/experts/experts/{expert_id}/feedback-summary")

        assert response.status_code == 200
        data = response.json()
        assert data["expert_name"] == "Dr. Feedback"
        assert data["avg_rating"] == 4.5
        assert data["feedback_count"] == 2
        assert len(data["recent_corrections"]) == 1

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_get_feedback_summary_no_feedback(self, mock_get_db, mock_service):
        """Test feedback summary with no feedback data."""
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        
        expert_id = uuid4()
        mock_expert = MagicMock()
        mock_expert.name = "Dr. NoFeedback"
        mock_service.get_expert.return_value = mock_expert
        
        mock_db.query.return_value.join.return_value.filter.return_value.all.return_value = []

        response = client.get(f"/api/v1/experts/experts/{expert_id}/feedback-summary")

        assert response.status_code == 200
        data = response.json()
        assert data["avg_rating"] is None
        assert data["feedback_count"] == 0


class TestTrainingAPIIntegration:
    """Integration tests for training API."""

    def test_training_endpoints_registered(self):
        """Test that training endpoints are registered."""
        from amprenta_rag.api.main import app
        routes = [route.path for route in app.routes]
        
        # Should have training endpoints
        training_routes = [r for r in routes if "/experts" in r and ("training" in r or "knowledge" in r or "feedback" in r)]
        assert len(training_routes) > 0

    def test_all_expert_endpoints_count(self):
        """Test total expert endpoint count."""
        from amprenta_rag.api.main import app
        
        expert_endpoints = []
        for route in app.routes:
            if hasattr(route, 'path') and '/experts' in route.path:
                if hasattr(route, 'methods'):
                    for method in route.methods:
                        if method != 'HEAD':  # Skip HEAD methods
                            expert_endpoints.append(f"{method} {route.path}")
        
        # Should have 8+ chat endpoints + 6 training endpoints = 14+ total
        assert len(expert_endpoints) >= 14

    def test_service_functions_available(self):
        """Test that required service functions are available."""
        from amprenta_rag.services import expert_agents as service
        
        training_functions = [
            'create_training_example',
            'approve_training_example',
            'export_training_data',
            'add_knowledge_doc',
            'update_expert',
        ]
        
        for func_name in training_functions:
            assert hasattr(service, func_name)
            assert callable(getattr(service, func_name))
