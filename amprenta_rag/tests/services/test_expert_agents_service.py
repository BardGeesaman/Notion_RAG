"""Unit tests for expert agents service."""

import pytest
from datetime import datetime
from uuid import uuid4
from unittest.mock import MagicMock, patch

from amprenta_rag.services import expert_agents as service


class TestExpertManagement:
    """Test expert management functions."""
    
    def test_get_experts_active_only(self):
        """Test get_experts with active filter."""
        mock_db = MagicMock()
        mock_experts = [MagicMock(), MagicMock()]
        mock_db.query.return_value.filter.return_value.all.return_value = mock_experts
        
        result = service.get_experts(mock_db, active_only=True)
        
        assert result == mock_experts
        mock_db.query.assert_called_once()
    
    def test_get_experts_all(self):
        """Test get_experts without filter."""
        mock_db = MagicMock()
        mock_experts = [MagicMock(), MagicMock(), MagicMock()]
        mock_db.query.return_value.all.return_value = mock_experts
        
        result = service.get_experts(mock_db, active_only=False)
        
        assert result == mock_experts
    
    def test_get_expert_found(self):
        """Test get_expert by ID."""
        mock_db = MagicMock()
        mock_expert = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_expert
        
        result = service.get_expert(mock_db, uuid4())
        
        assert result == mock_expert
    
    def test_get_expert_not_found(self):
        """Test get_expert returns None for non-existent ID."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        result = service.get_expert(mock_db, uuid4())
        
        assert result is None
    
    def test_update_expert_prompt_versioning(self):
        """Test expert update with version bump."""
        mock_db = MagicMock()
        mock_expert = MagicMock()
        mock_expert.prompt_version = "1.0"
        
        with patch.object(service, 'get_expert', return_value=mock_expert):
            result = service.update_expert(
                mock_db,
                uuid4(),
                system_prompt="New prompt",
                bump_version=True
            )
            
            assert result == mock_expert
            assert mock_expert.system_prompt == "New prompt"
            assert mock_expert.prompt_version == "1.1"
            mock_db.commit.assert_called_once()


class TestConversationManagement:
    """Test conversation management functions."""
    
    def test_create_conversation_single_expert(self):
        """Test creating conversation with single expert."""
        mock_db = MagicMock()
        mock_conversation = MagicMock()
        mock_expert = MagicMock()
        
        # Mock the conversation creation
        with patch('amprenta_rag.services.expert_agents.ExpertConversation', return_value=mock_conversation):
            mock_db.query.return_value.filter.return_value.all.return_value = [mock_expert]
            
            result = service.create_conversation(
                mock_db,
                user_id=uuid4(),
                expert_ids=[uuid4()],
                title="Test Conversation"
            )
            
            assert result == mock_conversation
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
    
    def test_create_conversation_panel_mode(self):
        """Test creating panel conversation with multiple experts."""
        mock_db = MagicMock()
        mock_conversation = MagicMock()
        mock_experts = [MagicMock(), MagicMock()]
        
        with patch('amprenta_rag.services.expert_agents.ExpertConversation', return_value=mock_conversation):
            mock_db.query.return_value.filter.return_value.all.return_value = mock_experts
            
            result = service.create_conversation(
                mock_db,
                user_id=uuid4(),
                expert_ids=[uuid4(), uuid4()],
                title="Panel Discussion",
                is_panel=True
            )
            
            assert result == mock_conversation
    
    def test_get_conversation(self):
        """Test get_conversation by ID."""
        mock_db = MagicMock()
        mock_conversation = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_conversation
        
        result = service.get_conversation(mock_db, uuid4())
        
        assert result == mock_conversation
    
    def test_get_user_conversations(self):
        """Test listing user conversations."""
        mock_db = MagicMock()
        mock_conversations = [MagicMock(), MagicMock()]
        mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_conversations
        
        result = service.get_user_conversations(mock_db, uuid4())
        
        assert result == mock_conversations


class TestMessageGeneration:
    """Test message generation functions."""
    
    def test_add_user_message(self):
        """Test adding user message."""
        mock_db = MagicMock()
        mock_message = MagicMock()
        
        with patch('amprenta_rag.services.expert_agents.ExpertMessage', return_value=mock_message):
            result = service.add_user_message(
                mock_db,
                uuid4(),
                "Test message"
            )
            
            assert result == mock_message
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
    
    @patch('amprenta_rag.services.expert_agents.get_llm_response_with_tokens')
    @patch('amprenta_rag.services.expert_agents.retrieve_expert_knowledge')
    @patch('amprenta_rag.services.expert_agents.inject_entity_context')
    def test_get_expert_response(self, mock_inject, mock_retrieve, mock_llm):
        """Test expert response generation."""
        mock_db = MagicMock()
        
        # Mock conversation and expert
        mock_conversation = MagicMock()
        mock_conversation.context_entity_type = None
        mock_expert = MagicMock()
        mock_expert.system_prompt = "You are an expert"
        mock_expert.prompt_version = "1.0"
        
        # Mock dependencies
        mock_retrieve.return_value = []
        mock_inject.return_value = None
        mock_llm.return_value = ("Expert response", 150)
        
        # Mock message history
        mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = []
        
        mock_message = MagicMock()
        with patch.object(service, 'get_conversation', return_value=mock_conversation):
            with patch.object(service, 'get_expert', return_value=mock_expert):
                with patch('amprenta_rag.services.expert_agents.ExpertMessage', return_value=mock_message):
                    result = service.get_expert_response(
                        mock_db,
                        uuid4(),
                        uuid4(),
                        "Test question"
                    )
                    
                    assert result == mock_message
                    mock_llm.assert_called_once()
    
    def test_get_expert_response_extracts_reasoning(self):
        """Test reasoning extraction from response."""
        mock_db = MagicMock()
        mock_conversation = MagicMock()
        mock_conversation.context_entity_type = None
        mock_expert = MagicMock()
        mock_expert.system_prompt = "You are an expert"
        mock_expert.prompt_version = "1.0"
        
        # Mock LLM response with reasoning tags
        llm_response = "Here is my answer. <reasoning>This is my reasoning</reasoning> Final answer."
        
        with patch.object(service, 'get_conversation', return_value=mock_conversation):
            with patch.object(service, 'get_expert', return_value=mock_expert):
                with patch('amprenta_rag.services.expert_agents.get_llm_response_with_tokens', return_value=(llm_response, 100)):
                    with patch('amprenta_rag.services.expert_agents.retrieve_expert_knowledge', return_value=[]):
                        with patch('amprenta_rag.services.expert_agents.ExpertMessage') as mock_message_class:
                            mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = []
                            
                            service.get_expert_response(mock_db, uuid4(), uuid4(), "Question")
                            
                            # Check that reasoning was extracted
                            call_args = mock_message_class.call_args[1]
                            assert call_args['reasoning'] == "This is my reasoning"
                            assert "<reasoning>" not in call_args['content']


class TestRAGIntegration:
    """Test RAG knowledge retrieval."""
    
    def test_retrieve_expert_knowledge(self):
        """Test knowledge retrieval for expert."""
        mock_db = MagicMock()
        mock_expert = MagicMock()
        mock_docs = [MagicMock(), MagicMock()]
        
        with patch.object(service, 'get_expert', return_value=mock_expert):
            with patch('amprenta_rag.services.expert_agents.get_embedding', return_value=[0.1, 0.2]):
                mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_docs
                
                result = service.retrieve_expert_knowledge(mock_db, uuid4(), "test query")
                
                assert result == mock_docs
    
    def test_retrieve_expert_knowledge_fallback(self):
        """Test knowledge retrieval fallback on error."""
        mock_db = MagicMock()
        mock_expert = MagicMock()
        mock_docs = [MagicMock()]
        
        with patch.object(service, 'get_expert', return_value=mock_expert):
            with patch('amprenta_rag.services.expert_agents.get_embedding', side_effect=Exception("Embedding failed")):
                # Mock the fallback query
                mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_docs
                
                result = service.retrieve_expert_knowledge(mock_db, uuid4(), "test query")
                
                assert result == mock_docs
    
    def test_add_knowledge_doc(self):
        """Test adding knowledge document."""
        mock_db = MagicMock()
        mock_doc = MagicMock()
        
        with patch('amprenta_rag.services.expert_agents.get_embedding', return_value=[0.1, 0.2]):
            with patch('amprenta_rag.services.expert_agents.ExpertKnowledgeDoc', return_value=mock_doc):
                result = service.add_knowledge_doc(
                    mock_db,
                    uuid4(),
                    "Test Doc",
                    "Test content",
                    source_type="manual"
                )
                
                assert result == mock_doc
                mock_db.add.assert_called_once()
                mock_db.commit.assert_called_once()


class TestFeedbackAndTraining:
    """Test feedback and training functions."""
    
    def test_record_feedback_valid_rating(self):
        """Test recording feedback with valid rating."""
        mock_db = MagicMock()
        mock_feedback = MagicMock()
        
        with patch('amprenta_rag.services.expert_agents.ExpertFeedback', return_value=mock_feedback):
            result = service.record_feedback(
                mock_db,
                uuid4(),
                uuid4(),
                rating=4,
                correction="Minor improvement",
                tags=["helpful"]
            )
            
            assert result == mock_feedback
            mock_db.add.assert_called_once()
    
    def test_record_feedback_invalid_rating(self):
        """Test invalid rating raises ValueError."""
        mock_db = MagicMock()
        
        with pytest.raises(ValueError, match="Rating must be between 1 and 5"):
            service.record_feedback(mock_db, uuid4(), uuid4(), rating=6)
        
        with pytest.raises(ValueError, match="Rating must be between 1 and 5"):
            service.record_feedback(mock_db, uuid4(), uuid4(), rating=0)
    
    def test_create_training_example(self):
        """Test creating training example."""
        mock_db = MagicMock()
        mock_expert = MagicMock()
        mock_expert.prompt_version = "1.0"
        mock_example = MagicMock()
        
        with patch.object(service, 'get_expert', return_value=mock_expert):
            with patch('amprenta_rag.services.expert_agents.ExpertTrainingExample', return_value=mock_example):
                result = service.create_training_example(
                    mock_db,
                    uuid4(),
                    "Test question?",
                    "Test answer",
                    uuid4()
                )
                
                assert result == mock_example
                mock_db.add.assert_called_once()
    
    def test_approve_training_example(self):
        """Test approving training example."""
        mock_db = MagicMock()
        mock_example = MagicMock()
        mock_example.is_approved = False
        mock_db.query.return_value.filter.return_value.first.return_value = mock_example
        
        result = service.approve_training_example(mock_db, uuid4())
        
        assert result is True
        assert mock_example.is_approved is True
        mock_db.commit.assert_called_once()
    
    def test_export_training_data(self):
        """Test exporting training data in OpenAI format."""
        mock_db = MagicMock()
        
        # Mock training example
        mock_example = MagicMock()
        mock_example.expert.system_prompt = "You are an expert"
        mock_example.question = "Test question?"
        mock_example.ideal_answer = "Test answer"
        
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_example]
        
        result = service.export_training_data(mock_db)
        
        assert len(result) == 1
        assert "messages" in result[0]
        assert len(result[0]["messages"]) == 3  # system, user, assistant
        assert result[0]["messages"][0]["role"] == "system"
        assert result[0]["messages"][1]["role"] == "user"
        assert result[0]["messages"][2]["role"] == "assistant"


class TestEntityContextInjection:
    """Test entity context injection."""
    
    def test_inject_entity_context_compound(self):
        """Test compound context injection."""
        mock_db = MagicMock()
        mock_compound = MagicMock()
        mock_compound.compound_id = "COMP-001"
        mock_compound.smiles = "CCO"
        mock_compound.molecular_weight = 46.07
        
        mock_db.query.return_value.filter.return_value.first.return_value = mock_compound
        
        result = service.inject_entity_context(mock_db, "Compound", uuid4())
        
        assert "COMP-001" in result
        assert "CCO" in result
        assert "46.07" in result
    
    def test_inject_entity_context_dataset(self):
        """Test dataset context injection."""
        mock_db = MagicMock()
        mock_dataset = MagicMock()
        mock_dataset.name = "Test Dataset"
        mock_dataset.omics_type = "transcriptomics"
        mock_dataset.description = "RNA-seq data"
        
        mock_db.query.return_value.filter.return_value.first.return_value = mock_dataset
        
        result = service.inject_entity_context(mock_db, "Dataset", uuid4())
        
        assert "Test Dataset" in result
        assert "transcriptomics" in result
        assert "RNA-seq data" in result
    
    def test_inject_entity_context_not_found(self):
        """Test entity context returns None for unknown entity."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        result = service.inject_entity_context(mock_db, "Compound", uuid4())
        
        assert result is None
    
    def test_inject_entity_context_unknown_type(self):
        """Test unknown entity type returns None."""
        mock_db = MagicMock()
        
        result = service.inject_entity_context(mock_db, "UnknownType", uuid4())
        
        assert result is None


class TestServiceUtilities:
    """Test utility functions."""
    
    def test_get_conversation_messages(self):
        """Test getting conversation messages."""
        mock_db = MagicMock()
        mock_messages = [MagicMock(), MagicMock()]
        mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_messages
        
        result = service.get_conversation_messages(mock_db, uuid4())
        
        assert result == mock_messages
    
    def test_get_expert_statistics(self):
        """Test expert statistics calculation."""
        mock_db = MagicMock()
        mock_expert = MagicMock()
        mock_expert.name = "Dr. Test"
        mock_expert.role = "Test Expert"
        mock_expert.prompt_version = "1.0"
        
        # Mock various counts
        mock_db.query.return_value.filter.return_value.count.side_effect = [10, 5, 3]  # messages, knowledge, training
        
        # Mock feedback ratings
        mock_feedback1 = MagicMock()
        mock_feedback1.rating = 4
        mock_feedback2 = MagicMock()
        mock_feedback2.rating = 5
        mock_db.query.return_value.join.return_value.filter.return_value.all.return_value = [mock_feedback1, mock_feedback2]
        
        with patch.object(service, 'get_expert', return_value=mock_expert):
            result = service.get_expert_statistics(mock_db, uuid4())
            
            assert result["expert_name"] == "Dr. Test"
            assert result["message_count"] == 10
            assert result["knowledge_doc_count"] == 5
            assert result["training_example_count"] == 3
            assert result["average_rating"] == 4.5
            assert result["total_feedback"] == 2


class TestPanelOrchestration:
    """Test panel response orchestration."""
    
    def test_detect_disagreements(self):
        """Test disagreement detection function."""
        responses = [
            {"content": "I recommend approach A"},
            {"content": "I would not recommend approach A"},
        ]
        
        result = service.detect_disagreements(responses)
        
        # MVP implementation returns empty list
        assert isinstance(result, list)
    
    def test_synthesize_recommendations(self):
        """Test recommendation synthesis."""
        responses = [
            {"content": "First expert recommendation with details"},
            {"content": "Second expert recommendation"},
        ]
        
        result = service.synthesize_recommendations(responses)
        
        assert "Based on expert consensus" in result
        assert "First expert recommendation" in result
    
    def test_synthesize_recommendations_empty(self):
        """Test synthesis with no responses."""
        result = service.synthesize_recommendations([])
        
        assert result == "No recommendations available."
