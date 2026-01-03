"""Tests for expert agent models."""

import pytest
from datetime import datetime
from uuid import uuid4
from unittest.mock import MagicMock

from amprenta_rag.database.models import (
    ExpertAgent, ExpertConversation, ExpertMessage, ExpertFeedback,
    ExpertTrainingExample, ExpertKnowledgeDoc
)


class TestExpertAgentModel:
    """Test ExpertAgent model."""
    
    def test_expert_agent_creation(self):
        """Test creating an expert agent."""
        expert = ExpertAgent(
            name="Dr. Test",
            role="Test Expert",
            system_prompt="You are a test expert",
            prompt_version="1.0",
            specializations=["testing", "validation"],
        )
        
        assert expert.name == "Dr. Test"
        assert expert.role == "Test Expert"
        assert expert.prompt_version == "1.0"
        assert expert.specializations == ["testing", "validation"]
        # Note: is_active default may not be set until DB insert
    
    def test_expert_agent_prompt_versioning(self):
        """Test prompt versioning functionality."""
        expert = ExpertAgent(
            name="Dr. Version",
            role="Version Expert",
            system_prompt="Version 1 prompt",
            prompt_version="1.0",
            specializations=[],
        )
        
        assert expert.prompt_version == "1.0"
        
        # Update to new version
        expert.prompt_version = "1.1"
        expert.system_prompt = "Version 1.1 prompt"
        
        assert expert.prompt_version == "1.1"
        assert expert.system_prompt == "Version 1.1 prompt"
    
    def test_expert_agent_specializations_jsonb(self):
        """Test specializations stored as JSONB."""
        specializations = ["SAR analysis", "ADMET prediction", "molecular design"]
        expert = ExpertAgent(
            name="Dr. Specialist",
            role="Specialist",
            system_prompt="Prompt",
            specializations=specializations,
        )
        
        assert expert.specializations == specializations
        assert isinstance(expert.specializations, list)


class TestExpertConversationModel:
    """Test ExpertConversation model."""
    
    def test_conversation_creation(self):
        """Test creating a conversation."""
        conversation = ExpertConversation(
            title="Test Conversation",
            context_entity_type="Compound",
            context_entity_id=uuid4(),
            max_messages=50,
            is_panel=False,
        )
        
        assert conversation.title == "Test Conversation"
        assert conversation.context_entity_type == "Compound"
        assert conversation.max_messages == 50
        assert conversation.is_panel is False
        # Note: is_active default may not be set until DB insert
    
    def test_conversation_max_messages_limit(self):
        """Test conversation message limit (P1 fix)."""
        conversation = ExpertConversation(
            title="Limited Conversation",
            max_messages=10,
        )
        
        assert conversation.max_messages == 10
    
    def test_conversation_panel_mode(self):
        """Test panel conversation mode."""
        panel_conversation = ExpertConversation(
            title="Expert Panel",
            is_panel=True,
            max_messages=100,
        )
        
        assert panel_conversation.is_panel is True


class TestExpertMessageModel:
    """Test ExpertMessage model."""
    
    def test_message_creation(self):
        """Test creating a message."""
        message = ExpertMessage(
            conversation_id=uuid4(),
            expert_id=uuid4(),
            role="assistant",
            content="This is a test response",
            prompt_version_used="1.0",
            reasoning="Test reasoning",
            citations=["source1", "source2"],
            token_count=150,
        )
        
        assert message.role == "assistant"
        assert message.content == "This is a test response"
        assert message.prompt_version_used == "1.0"
        assert message.reasoning == "Test reasoning"
        assert message.citations == ["source1", "source2"]
        assert message.token_count == 150
    
    def test_message_token_tracking(self):
        """Test token count tracking (P1 fix)."""
        message = ExpertMessage(
            conversation_id=uuid4(),
            role="assistant",
            content="Short response",
            token_count=25,
        )
        
        assert message.token_count == 25
    
    def test_message_prompt_version_tracking(self):
        """Test prompt version tracking (P1 fix)."""
        message = ExpertMessage(
            conversation_id=uuid4(),
            role="assistant",
            content="Response",
            prompt_version_used="1.2",
        )
        
        assert message.prompt_version_used == "1.2"


class TestExpertFeedbackModel:
    """Test ExpertFeedback model."""
    
    def test_feedback_creation(self):
        """Test creating feedback."""
        feedback = ExpertFeedback(
            message_id=uuid4(),
            user_id=uuid4(),
            rating=4,
            correction="Minor correction needed",
            tags=["helpful", "accurate"],
        )
        
        assert feedback.rating == 4
        assert feedback.correction == "Minor correction needed"
        assert feedback.tags == ["helpful", "accurate"]
    
    def test_feedback_rating_schema(self):
        """Test rating is 1-5 integer (P1 fix)."""
        # Test valid ratings
        for rating in [1, 2, 3, 4, 5]:
            feedback = ExpertFeedback(
                message_id=uuid4(),
                rating=rating,
            )
            assert feedback.rating == rating
            assert isinstance(feedback.rating, int)


class TestExpertTrainingExampleModel:
    """Test ExpertTrainingExample model."""
    
    def test_training_example_creation(self):
        """Test creating training example."""
        example = ExpertTrainingExample(
            expert_id=uuid4(),
            question="What is the best approach for lead optimization?",
            ideal_answer="Focus on ADMET properties while maintaining potency...",
            is_approved=False,
            prompt_version="1.0",
            source="manual",
        )
        
        assert example.question == "What is the best approach for lead optimization?"
        assert example.is_approved is False
        assert example.prompt_version == "1.0"
        assert example.source == "manual"
    
    def test_training_approval_workflow(self):
        """Test training approval workflow (P1 fix)."""
        example = ExpertTrainingExample(
            expert_id=uuid4(),
            question="Test question?",
            ideal_answer="Test answer",
            is_approved=False,
            prompt_version="1.0",
        )
        
        # Initially not approved
        assert example.is_approved is False
        
        # Approve example
        example.is_approved = True
        assert example.is_approved is True


class TestExpertKnowledgeDocModel:
    """Test ExpertKnowledgeDoc model."""
    
    def test_knowledge_doc_creation(self):
        """Test creating knowledge document."""
        doc = ExpertKnowledgeDoc(
            expert_id=uuid4(),
            namespace="chemistry",
            title="SAR Principles",
            content="Structure-activity relationships are...",
            chunk_index=0,
            embedding=[0.1, 0.2, 0.3],
            embedding_model="text-embedding-ada-002",
            source_url="https://example.com/sar",
            source_type="web",
        )
        
        assert doc.namespace == "chemistry"
        assert doc.title == "SAR Principles"
        assert doc.chunk_index == 0
        assert doc.embedding == [0.1, 0.2, 0.3]
        assert doc.embedding_model == "text-embedding-ada-002"
        assert doc.source_type == "web"
    
    def test_knowledge_doc_rag_namespace(self):
        """Test RAG namespace functionality (P1 fix)."""
        doc = ExpertKnowledgeDoc(
            expert_id=uuid4(),
            namespace="admet_prediction",
            title="ADMET Models",
            content="ADMET prediction models...",
            chunk_index=5,
            embedding_model="all-MiniLM-L6-v2",
        )
        
        assert doc.namespace == "admet_prediction"
        assert doc.chunk_index == 5
        assert doc.embedding_model == "all-MiniLM-L6-v2"


class TestModelRelationships:
    """Test model relationships."""
    
    def test_expert_agent_relationships(self):
        """Test ExpertAgent has proper relationships."""
        expert = ExpertAgent(
            name="Dr. Relations",
            role="Relationship Expert",
            system_prompt="Testing relationships",
            specializations=[],
        )
        
        # Check relationship attributes exist
        assert hasattr(expert, 'conversations')
        assert hasattr(expert, 'messages')
        assert hasattr(expert, 'training_examples')
        assert hasattr(expert, 'knowledge_docs')
    
    def test_conversation_relationships(self):
        """Test ExpertConversation has proper relationships."""
        conversation = ExpertConversation(
            title="Test Conversation",
            max_messages=50,
        )
        
        # Check relationship attributes exist
        assert hasattr(conversation, 'participants')
        assert hasattr(conversation, 'messages')
        assert hasattr(conversation, 'user')
    
    def test_message_relationships(self):
        """Test ExpertMessage has proper relationships."""
        message = ExpertMessage(
            conversation_id=uuid4(),
            role="assistant",
            content="Test message",
        )
        
        # Check relationship attributes exist
        assert hasattr(message, 'conversation')
        assert hasattr(message, 'expert')
        assert hasattr(message, 'feedback')
