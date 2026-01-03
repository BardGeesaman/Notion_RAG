"""Tests for AI Expert Agents service with CrewAI."""

import pytest
from uuid import uuid4
from unittest.mock import MagicMock, patch

from amprenta_rag.services import ai_expert_agents as service


class TestCrewAIToolInitialization:
    """Test CrewAI tool initialization."""
    
    def test_compound_lookup_tool_creation(self):
        """Test CompoundLookupTool can be created."""
        tool = service.CompoundLookupTool()
        assert tool.name == "compound_lookup"
        assert "compound" in tool.description.lower()
    
    def test_rag_search_tool_creation(self):
        """Test RAGSearchTool with expert ID."""
        expert_id = uuid4()
        tool = service.RAGSearchTool(expert_id)
        assert tool.name == "rag_search"
        assert tool.expert_id == expert_id
    
    def test_dataset_info_tool_creation(self):
        """Test DatasetInfoTool creation."""
        tool = service.DatasetInfoTool()
        assert tool.name == "dataset_info"
        assert "dataset" in tool.description.lower()
    
    def test_experiment_info_tool_creation(self):
        """Test ExperimentInfoTool creation."""
        tool = service.ExperimentInfoTool()
        assert tool.name == "experiment_info"
        assert "experiment" in tool.description.lower()


class TestCrewAIAgentCreation:
    """Test CrewAI agent creation from database models."""
    
    def test_create_expert_agent(self):
        """Test creating CrewAI agent from ExpertAgent model."""
        mock_expert = MagicMock()
        mock_expert.id = uuid4()
        mock_expert.role = "Medicinal Chemist"
        mock_expert.system_prompt = "You are an expert medicinal chemist"
        
        with patch('crewai.Agent') as mock_agent_class:
            mock_agent = MagicMock()
            mock_agent_class.return_value = mock_agent
            
            result = service.create_expert_agent(mock_expert)
            
            assert result == mock_agent
            mock_agent_class.assert_called_once()
            
            # Check agent was created with correct parameters
            call_args = mock_agent_class.call_args[1]
            assert call_args['role'] == "Medicinal Chemist"
            assert call_args['backstory'] == "You are an expert medicinal chemist"
            assert 'tools' in call_args
            assert call_args['memory'] is True
    
    def test_create_expert_crew(self):
        """Test creating CrewAI crew from multiple experts."""
        mock_experts = [MagicMock(), MagicMock()]
        
        with patch.object(service, 'create_expert_agent') as mock_create_agent:
            with patch('crewai.Crew') as mock_crew_class:
                mock_agents = [MagicMock(), MagicMock()]
                mock_create_agent.side_effect = mock_agents
                mock_crew = MagicMock()
                mock_crew_class.return_value = mock_crew
                
                result = service.create_expert_crew(mock_experts, process="hierarchical")
                
                assert result == mock_crew
                assert mock_create_agent.call_count == 2
                mock_crew_class.assert_called_once()


class TestResponseGeneration:
    """Test response generation with CrewAI."""
    
    @patch('crewai.Crew')
    @patch('crewai.Task')
    def test_get_expert_response(self, mock_task_class, mock_crew_class):
        """Test single expert response generation."""
        mock_db = MagicMock()
        
        # Mock expert and conversation
        mock_expert = MagicMock()
        mock_expert.id = uuid4()
        mock_expert.name = "Dr. Test"
        mock_expert.prompt_version = "1.0"
        
        mock_conversation = MagicMock()
        mock_conversation.context_entity_type = None
        
        mock_db.query.return_value.filter.return_value.first.side_effect = [mock_expert, mock_conversation]
        mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = []
        
        # Mock CrewAI execution
        mock_crew = MagicMock()
        mock_crew.kickoff.return_value = "Expert response content"
        mock_crew_class.return_value = mock_crew
        
        mock_task = MagicMock()
        mock_task_class.return_value = mock_task
        
        mock_message = MagicMock()
        with patch('amprenta_rag.services.ai_expert_agents.ExpertMessage', return_value=mock_message):
            with patch.object(service, 'create_expert_agent'):
                result = service.get_expert_response(
                    mock_db,
                    uuid4(),
                    uuid4(),
                    "Test question"
                )
                
                assert result == mock_message
                mock_crew.kickoff.assert_called_once()
                mock_db.add.assert_called()
                mock_db.commit.assert_called()
    
    @patch('crewai.Crew')
    @patch('crewai.Task')
    def test_get_panel_response(self, mock_task_class, mock_crew_class):
        """Test panel response generation."""
        mock_db = MagicMock()
        
        # Mock experts and conversation
        mock_experts = [MagicMock(), MagicMock()]
        for expert in mock_experts:
            expert.id = uuid4()
            expert.name = "Dr. Test"
            expert.prompt_version = "1.0"
        
        mock_conversation = MagicMock()
        mock_conversation.context_entity_type = None
        
        mock_db.query.return_value.filter.return_value.all.return_value = mock_experts
        mock_db.query.return_value.filter.return_value.first.return_value = mock_conversation
        
        # Mock CrewAI execution
        mock_crew = MagicMock()
        mock_crew.kickoff.return_value = "Panel consensus response"
        mock_crew_class.return_value = mock_crew
        
        with patch.object(service, 'create_expert_crew', return_value=mock_crew):
            with patch('amprenta_rag.services.ai_expert_agents.ExpertMessage'):
                result = service.get_panel_response(
                    mock_db,
                    uuid4(),
                    [uuid4(), uuid4()],
                    "Panel question"
                )
                
                assert "consensus" in result
                assert "expert_responses" in result
                assert result["primary_recommendation"] == "Panel consensus response"


class TestContextInjection:
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
    
    def test_inject_entity_context_unknown_type(self):
        """Test unknown entity type returns None."""
        mock_db = MagicMock()
        
        result = service.inject_entity_context(mock_db, "UnknownType", uuid4())
        
        assert result is None


class TestFeedbackAndTraining:
    """Test feedback and training functions."""
    
    def test_record_feedback_rating(self):
        """Test recording rating feedback."""
        mock_db = MagicMock()
        mock_feedback = MagicMock()
        
        with patch('amprenta_rag.services.ai_expert_agents.ExpertFeedback', return_value=mock_feedback):
            result = service.record_feedback(
                mock_db,
                uuid4(),
                feedback_type="rating",
                feedback_value=4,
                comments="Good response"
            )
            
            assert result == mock_feedback
            mock_db.add.assert_called_once()
    
    def test_record_feedback_thumbs(self):
        """Test recording thumbs up/down feedback."""
        mock_db = MagicMock()
        mock_feedback = MagicMock()
        
        with patch('amprenta_rag.services.ai_expert_agents.ExpertFeedback', return_value=mock_feedback):
            # Thumbs up
            service.record_feedback(mock_db, uuid4(), "thumbs", "up")
            
            # Thumbs down
            service.record_feedback(mock_db, uuid4(), "thumbs", "down")
            
            assert mock_db.add.call_count == 2
    
    def test_record_feedback_invalid_rating(self):
        """Test invalid rating raises ValueError."""
        mock_db = MagicMock()
        
        with pytest.raises(ValueError, match="Rating must be between 1 and 5"):
            service.record_feedback(mock_db, uuid4(), "rating", 6)
    
    def test_export_training_data_jsonl(self):
        """Test exporting training data in JSONL format."""
        mock_db = MagicMock()
        
        mock_example = MagicMock()
        mock_example.expert.system_prompt = "You are an expert"
        mock_example.question = "Test question?"
        mock_example.ideal_answer = "Test answer"
        
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_example]
        
        result = service.export_training_data(mock_db, format='jsonl')
        
        assert len(result) == 1
        assert "messages" in result[0]
        assert len(result[0]["messages"]) == 3
    
    def test_export_training_data_raw(self):
        """Test exporting training data in raw format."""
        mock_db = MagicMock()
        
        mock_example = MagicMock()
        mock_example.expert.name = "Dr. Test"
        mock_example.question = "Question"
        mock_example.ideal_answer = "Answer"
        mock_example.prompt_version = "1.0"
        
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_example]
        
        result = service.export_training_data(mock_db, format='raw')
        
        assert len(result) == 1
        assert result[0]["expert_name"] == "Dr. Test"
        assert result[0]["question"] == "Question"


class TestRAGIntegration:
    """Test RAG knowledge retrieval."""
    
    def test_retrieve_expert_knowledge_success(self):
        """Test successful knowledge retrieval."""
        mock_db = MagicMock()
        mock_docs = [MagicMock(), MagicMock()]
        
        with patch('amprenta_rag.services.ai_expert_agents.get_embedding', return_value=[0.1, 0.2]):
            mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_docs
            
            result = service.retrieve_expert_knowledge(mock_db, uuid4(), "test query")
            
            assert result == mock_docs
    
    def test_retrieve_expert_knowledge_fallback(self):
        """Test knowledge retrieval fallback on error."""
        mock_db = MagicMock()
        mock_docs = [MagicMock()]
        
        with patch('amprenta_rag.services.ai_expert_agents.get_embedding', side_effect=Exception("Embedding failed")):
            # Mock the fallback query
            mock_db.query.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = mock_docs
            
            result = service.retrieve_expert_knowledge(mock_db, uuid4(), "test query")
            
            assert result == mock_docs


class TestConversationManagement:
    """Test conversation management functions."""
    
    def test_create_conversation(self):
        """Test creating conversation."""
        mock_db = MagicMock()
        mock_conversation = MagicMock()
        mock_experts = [MagicMock(), MagicMock()]
        
        with patch('amprenta_rag.services.ai_expert_agents.ExpertConversation', return_value=mock_conversation):
            mock_db.query.return_value.filter.return_value.all.return_value = mock_experts
            
            result = service.create_conversation(
                mock_db,
                uuid4(),
                [uuid4(), uuid4()],
                title="Test Panel"
            )
            
            assert result == mock_conversation
            mock_db.add.assert_called_once()
    
    def test_get_conversation(self):
        """Test getting conversation by ID."""
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
    
    def test_add_user_message(self):
        """Test adding user message."""
        mock_db = MagicMock()
        mock_message = MagicMock()
        
        with patch('amprenta_rag.services.ai_expert_agents.ExpertMessage', return_value=mock_message):
            result = service.add_user_message(mock_db, uuid4(), "Test message")
            
            assert result == mock_message
            mock_db.add.assert_called_once()


class TestExpertUtilities:
    """Test expert utility functions."""
    
    def test_get_experts(self):
        """Test listing experts."""
        mock_db = MagicMock()
        mock_experts = [MagicMock(), MagicMock()]
        mock_db.query.return_value.filter.return_value.order_by.return_value.all.return_value = mock_experts
        
        result = service.get_experts(mock_db, active_only=True)
        
        assert result == mock_experts
    
    def test_get_expert(self):
        """Test getting single expert."""
        mock_db = MagicMock()
        mock_expert = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_expert
        
        result = service.get_expert(mock_db, uuid4())
        
        assert result == mock_expert
    
    def test_get_expert_statistics(self):
        """Test expert statistics calculation."""
        mock_db = MagicMock()
        mock_expert = MagicMock()
        mock_expert.name = "Dr. Test"
        mock_expert.role = "Test Expert"
        mock_expert.prompt_version = "1.0"
        mock_expert.specializations = ["testing"]
        
        # Mock counts
        mock_db.query.return_value.filter.return_value.count.side_effect = [10, 5]  # messages, knowledge
        
        # Mock feedback
        mock_feedback = MagicMock()
        mock_feedback.rating = 4
        mock_db.query.return_value.join.return_value.filter.return_value.all.return_value = [mock_feedback]
        
        with patch.object(service, 'get_expert', return_value=mock_expert):
            result = service.get_expert_statistics(mock_db, uuid4())
            
            assert result["expert_name"] == "Dr. Test"
            assert result["message_count"] == 10
            assert result["knowledge_doc_count"] == 5
            assert result["average_rating"] == 4.0


class TestServiceIntegration:
    """Test service integration and imports."""
    
    def test_crewai_imports(self):
        """Test that CrewAI imports successfully."""
        from crewai import Agent, Task, Crew
        assert Agent is not None
        assert Task is not None
        assert Crew is not None
    
    def test_service_functions_exist(self):
        """Test all expected service functions exist."""
        expected_functions = [
            'create_expert_agent',
            'create_expert_crew',
            'get_expert_response',
            'get_panel_response',
            'inject_entity_context',
            'record_feedback',
            'retrieve_expert_knowledge',
            'export_training_data',
            'create_conversation',
            'get_conversation',
            'get_user_conversations',
            'add_user_message',
            'get_experts',
            'get_expert',
            'get_expert_statistics',
        ]
        
        for func_name in expected_functions:
            assert hasattr(service, func_name)
            assert callable(getattr(service, func_name))
    
    def test_tool_classes_exist(self):
        """Test all tool classes exist."""
        tool_classes = [
            'CompoundLookupTool',
            'RAGSearchTool',
            'DatasetInfoTool',
            'ExperimentInfoTool',
        ]
        
        for tool_name in tool_classes:
            assert hasattr(service, tool_name)
            tool_class = getattr(service, tool_name)
            assert issubclass(tool_class, service.BaseTool)


class TestToolExecution:
    """Test tool execution (with mocking)."""
    
    def test_compound_lookup_tool_execution(self):
        """Test compound lookup tool execution."""
        tool = service.CompoundLookupTool()
        
        with patch('amprenta_rag.services.ai_expert_agents.db_session') as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            mock_compound = MagicMock()
            mock_compound.compound_id = "COMP-001"
            mock_compound.smiles = "CCO"
            mock_compound.molecular_weight = 46.07
            mock_db.query.return_value.filter.return_value.first.return_value = mock_compound
            
            result = tool._run(compound_id="COMP-001")
            
            assert "COMP-001" in result
            assert "CCO" in result
    
    def test_dataset_info_tool_execution(self):
        """Test dataset info tool execution."""
        tool = service.DatasetInfoTool()
        
        with patch('amprenta_rag.services.ai_expert_agents.db_session') as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            mock_dataset = MagicMock()
            mock_dataset.name = "Test Dataset"
            mock_dataset.omics_type = "transcriptomics"
            mock_db.query.return_value.filter.return_value.first.return_value = mock_dataset
            mock_db.query.return_value.filter.return_value.count.return_value = 1000
            
            result = tool._run(dataset_id=str(uuid4()))
            
            assert "Test Dataset" in result
            assert "transcriptomics" in result
    
    def test_tool_error_handling(self):
        """Test tool error handling."""
        tool = service.CompoundLookupTool()
        
        result = tool._run()  # No parameters
        
        assert "Error" in result
        assert "Must provide" in result
