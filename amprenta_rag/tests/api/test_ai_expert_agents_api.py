"""Tests for AI Expert Agents API endpoints."""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestExpertBrowsingEndpoints:
    """Tests for expert browsing endpoints."""

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_list_experts(self, mock_get_db, mock_service):
        """Test GET /experts returns expert list."""
        mock_get_db.return_value = MagicMock()
        
        mock_expert = MagicMock()
        mock_expert.id = uuid4()
        mock_expert.name = "Dr. Test"
        mock_expert.role = "Test Expert"
        mock_expert.specializations = ["testing", "validation"]
        mock_expert.is_active = True
        mock_expert.prompt_version = "1.0"
        
        mock_service.get_experts.return_value = [mock_expert]

        response = client.get("/api/v1/experts/")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["name"] == "Dr. Test"
        assert data[0]["role"] == "Test Expert"

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_get_expert_details(self, mock_get_db, mock_service):
        """Test GET /experts/{id} returns expert details."""
        mock_get_db.return_value = MagicMock()
        
        expert_id = uuid4()
        mock_expert = MagicMock()
        mock_expert.id = expert_id
        mock_expert.name = "Dr. Detailed"
        mock_expert.role = "Detailed Expert"
        mock_expert.system_prompt = "You are a detailed expert with comprehensive knowledge"
        mock_expert.specializations = ["details", "analysis"]
        mock_expert.is_active = True
        mock_expert.prompt_version = "1.0"
        
        mock_service.get_expert.return_value = mock_expert

        response = client.get(f"/api/v1/experts/{expert_id}")

        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(expert_id)
        assert data["name"] == "Dr. Detailed"

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_get_expert_not_found(self, mock_get_db, mock_service):
        """Test 404 for non-existent expert."""
        mock_get_db.return_value = MagicMock()
        mock_service.get_expert.return_value = None

        response = client.get(f"/api/v1/experts/{uuid4()}")

        assert response.status_code == 404

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_get_expert_statistics(self, mock_get_db, mock_service):
        """Test GET /experts/{id}/stats returns statistics."""
        mock_get_db.return_value = MagicMock()
        
        mock_stats = {
            "expert_name": "Dr. Stats",
            "expert_role": "Statistician",
            "message_count": 100,
            "knowledge_doc_count": 25,
            "training_example_count": 10,
            "average_rating": 4.5,
            "total_feedback": 20,
            "prompt_version": "1.0",
            "specializations": ["statistics"],
        }
        mock_service.get_expert_statistics.return_value = mock_stats

        response = client.get(f"/api/v1/experts/{uuid4()}/stats")

        assert response.status_code == 200
        data = response.json()
        assert data["expert_name"] == "Dr. Stats"
        assert data["message_count"] == 100


class TestConversationEndpoints:
    """Tests for conversation management endpoints."""

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_create_conversation(self, mock_get_db, mock_get_user, mock_service):
        """Test POST /conversations creates conversation."""
        mock_get_db.return_value = MagicMock()
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        mock_conversation = MagicMock()
        mock_conversation.id = uuid4()
        mock_conversation.title = "Test Conversation"
        mock_conversation.is_panel = False
        mock_conversation.participants = [MagicMock()]
        mock_conversation.participants[0].id = uuid4()
        mock_conversation.participants[0].name = "Dr. Test"
        mock_conversation.messages = []
        
        mock_service.create_conversation.return_value = mock_conversation

        response = client.post(
            "/api/v1/experts/conversations",
            json={
                "expert_ids": [str(uuid4())],
                "title": "Test Conversation"
            }
        )

        assert response.status_code == 201
        data = response.json()
        assert data["title"] == "Test Conversation"
        assert data["is_panel"] is False

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_list_user_conversations(self, mock_get_db, mock_get_user, mock_service):
        """Test GET /conversations returns user conversations."""
        mock_get_db.return_value = MagicMock()
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        mock_conversation = MagicMock()
        mock_conversation.id = uuid4()
        mock_conversation.title = "User Conversation"
        mock_conversation.participants = [MagicMock()]
        mock_conversation.participants[0].id = uuid4()
        mock_conversation.participants[0].name = "Dr. User"
        mock_conversation.messages = []
        mock_conversation.is_panel = False
        
        mock_service.get_user_conversations.return_value = [mock_conversation]

        response = client.get("/api/v1/experts/conversations")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["title"] == "User Conversation"

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_send_message(self, mock_get_db, mock_get_user, mock_service):
        """Test POST /conversations/{id}/messages sends message."""
        mock_get_db.return_value = MagicMock()
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        conversation_id = uuid4()
        mock_conversation = MagicMock()
        mock_conversation.user_id = mock_user.id
        mock_conversation.participants = [MagicMock()]
        mock_conversation.participants[0].id = uuid4()
        mock_conversation.participants[0].name = "Dr. Responder"
        
        mock_user_message = MagicMock()
        mock_expert_response = MagicMock()
        mock_expert_response.id = uuid4()
        mock_expert_response.role = "assistant"
        mock_expert_response.content = "Expert response"
        mock_expert_response.expert_id = uuid4()
        mock_expert_response.reasoning = None
        mock_expert_response.token_count = 50
        
        mock_service.get_conversation.return_value = mock_conversation
        mock_service.add_user_message.return_value = mock_user_message
        mock_service.get_expert_response.return_value = mock_expert_response

        response = client.post(
            f"/api/v1/experts/conversations/{conversation_id}/messages",
            json={"content": "Test question"}
        )

        assert response.status_code == 200
        data = response.json()
        assert data["content"] == "Expert response"
        assert data["role"] == "assistant"


class TestPanelEndpoints:
    """Tests for panel discussion endpoints."""

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_panel_discussion(self, mock_get_db, mock_get_user, mock_service):
        """Test POST /panel creates panel discussion."""
        mock_get_db.return_value = MagicMock()
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        mock_conversation = MagicMock()
        mock_conversation.id = uuid4()
        
        mock_panel_result = {
            "consensus": "partial",
            "confidence_score": 0.7,
            "primary_recommendation": "Panel recommends approach A",
            "expert_responses": [
                {
                    "expert_id": str(uuid4()),
                    "expert_name": "Dr. Panel1",
                    "expert_role": "Expert 1",
                    "content": "I recommend A",
                    "message_id": str(uuid4())
                },
                {
                    "expert_id": str(uuid4()),
                    "expert_name": "Dr. Panel2", 
                    "expert_role": "Expert 2",
                    "content": "I also recommend A",
                    "message_id": str(uuid4())
                }
            ],
            "disagreements": []
        }
        
        mock_service.create_conversation.return_value = mock_conversation
        mock_service.get_panel_response.return_value = mock_panel_result

        response = client.post(
            "/api/v1/experts/panel",
            json={
                "expert_ids": [str(uuid4()), str(uuid4())],
                "question": "What is the best approach?"
            }
        )

        assert response.status_code == 200
        data = response.json()
        assert data["consensus"] == "partial"
        assert data["confidence_score"] == 0.7
        assert len(data["responses"]) == 2

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_panel_discussion_insufficient_experts(self, mock_get_db, mock_get_user, mock_service):
        """Test panel requires at least 2 experts."""
        mock_get_db.return_value = MagicMock()
        mock_get_user.return_value = MagicMock()

        response = client.post(
            "/api/v1/experts/panel",
            json={
                "expert_ids": [str(uuid4())],  # Only 1 expert
                "question": "Test question"
            }
        )

        assert response.status_code == 400
        assert "at least 2 experts" in response.json()["detail"]


class TestFeedbackEndpoints:
    """Tests for feedback endpoints."""

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_submit_feedback(self, mock_get_db, mock_get_user, mock_service):
        """Test POST /messages/{id}/feedback submits feedback."""
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_get_user.return_value = mock_user
        
        message_id = uuid4()
        mock_message = MagicMock()
        mock_message.conversation.user_id = mock_user.id
        mock_db.query.return_value.filter.return_value.first.return_value = mock_message
        
        mock_feedback = MagicMock()
        mock_feedback.id = uuid4()
        mock_feedback.message_id = message_id
        mock_feedback.rating = 4
        mock_feedback.correction = "Good response"
        mock_feedback.tags = ["helpful"]
        
        mock_service.record_feedback.return_value = mock_feedback

        response = client.post(
            f"/api/v1/experts/messages/{message_id}/feedback",
            json={
                "rating": 4,
                "correction": "Good response",
                "tags": ["helpful"]
            }
        )

        assert response.status_code == 201
        data = response.json()
        assert data["rating"] == 4
        assert data["correction"] == "Good response"

    @patch("amprenta_rag.api.routers.ai_expert_agents.service")
    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_submit_feedback_invalid_rating(self, mock_get_db, mock_get_user, mock_service):
        """Test invalid rating returns validation error."""
        mock_get_db.return_value = MagicMock()
        mock_get_user.return_value = MagicMock()

        response = client.post(
            f"/api/v1/experts/messages/{uuid4()}/feedback",
            json={"rating": 6}  # Invalid rating
        )

        assert response.status_code == 422  # Validation error


class TestAPIIntegration:
    """Integration tests for expert agents API."""

    def test_expert_endpoints_registered(self):
        """Test that expert endpoints are registered."""
        from amprenta_rag.api.main import app
        routes = [route.path for route in app.routes]
        
        # Should have expert endpoints
        expert_routes = [r for r in routes if "/experts" in r]
        assert len(expert_routes) > 0

    def test_endpoint_count(self):
        """Test that all expected expert endpoints exist."""
        from amprenta_rag.api.main import app
        
        expert_endpoints = []
        for route in app.routes:
            if hasattr(route, 'path') and '/experts' in route.path:
                if hasattr(route, 'methods'):
                    for method in route.methods:
                        if method != 'HEAD':  # Skip HEAD methods
                            expert_endpoints.append(f"{method} {route.path}")
        
        # Should have at least 8 endpoints
        assert len(expert_endpoints) >= 8

    @patch("amprenta_rag.api.dependencies.get_current_user")
    @patch("amprenta_rag.api.dependencies.get_db")
    def test_authentication_required_for_conversations(self, mock_get_db, mock_get_user):
        """Test that conversation endpoints require authentication."""
        from fastapi import HTTPException
        mock_get_user.side_effect = HTTPException(status_code=401, detail="Not authenticated")
        
        response = client.post("/api/v1/experts/conversations", json={"expert_ids": [str(uuid4())]})
        assert response.status_code == 401


class TestServiceIntegration:
    """Test service layer integration."""

    def test_service_imports_successfully(self):
        """Test that expert agents service imports without error."""
        from amprenta_rag.services import expert_agents
        assert expert_agents is not None

    def test_schemas_import_successfully(self):
        """Test that schemas import without error."""
        from amprenta_rag.api import schemas_ai_expert_agents
        assert schemas_ai_expert_agents is not None

    def test_required_service_functions_exist(self):
        """Test that required service functions exist."""
        from amprenta_rag.services import expert_agents as service
        
        required_functions = [
            'get_experts',
            'get_expert',
            'create_conversation',
            'get_conversation',
            'get_user_conversations',
            'add_user_message',
            'get_expert_response',
            'get_panel_response',
            'record_feedback',
            'get_expert_statistics',
        ]
        
        for func_name in required_functions:
            assert hasattr(service, func_name)
            assert callable(getattr(service, func_name))
