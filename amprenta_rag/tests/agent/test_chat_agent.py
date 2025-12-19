"""
Tests for the chat agent module.

Coverage:
- Intent routing (chat_router.route_intent)
- Chat session state management (chat_types)
- Chat turn execution (chat_agent.run_chat_turn)
"""

from datetime import datetime
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.agent.chat_router import route_intent
from amprenta_rag.agent.chat_types import ChatMessage, ChatSessionState, ChatTurn


class TestChatRouter:
    """Tests for intent routing."""

    def test_route_similar_datasets(self):
        """'similar dataset' should route to similar_datasets intent."""
        assert route_intent("Find similar datasets to this one") == "similar_datasets"
        assert route_intent("Show me similar study data") == "similar_datasets"

    def test_route_dataset_summary(self):
        """Dataset summary keywords should route correctly."""
        assert route_intent("Summarize dataset ABC-123") == "dataset_summary"
        assert route_intent("Tell me about dataset xyz") == "dataset_summary"

    def test_route_program_summary(self):
        """Program summary keywords should route correctly."""
        assert route_intent("Give me a program summary") == "program_summary"
        assert route_intent("Program overview for ALS") == "program_summary"

    def test_route_signature_summary(self):
        """Signature keyword should route to signature_summary."""
        assert route_intent("What is this signature?") == "signature_summary"
        assert route_intent("Signature details for SIG-001") == "signature_summary"

    def test_route_help(self):
        """Help keywords should route to help intent."""
        assert route_intent("help") == "help"
        assert route_intent("What can you do?") == "help"

    def test_route_freeform_rag_default(self):
        """Unrecognized queries should default to freeform_rag."""
        assert route_intent("What genes are upregulated in ALS?") == "freeform_rag"
        assert route_intent("Show me the data") == "freeform_rag"


class TestChatTypes:
    """Tests for chat data types."""

    def test_chat_message_creation(self):
        """ChatMessage should be created with required fields."""
        msg = ChatMessage(
            role="user",
            content="Hello",
            timestamp=datetime.utcnow(),
        )
        assert msg.role == "user"
        assert msg.content == "Hello"
        assert isinstance(msg.timestamp, datetime)

    def test_chat_turn_creation(self):
        """ChatTurn should hold a list of messages."""
        user_msg = ChatMessage(role="user", content="Hi", timestamp=datetime.utcnow())
        assistant_msg = ChatMessage(role="assistant", content="Hello!", timestamp=datetime.utcnow())
        turn = ChatTurn(messages=[user_msg, assistant_msg])
        assert len(turn.messages) == 2
        assert turn.messages[0].role == "user"
        assert turn.messages[1].role == "assistant"

    def test_chat_session_state_creation(self):
        """ChatSessionState should initialize with empty turns."""
        session = ChatSessionState(
            id=uuid4(),
            created_at=datetime.utcnow(),
        )
        assert session.turns == []

    def test_chat_session_state_with_turns(self):
        """ChatSessionState should accept pre-populated turns."""
        turn = ChatTurn(
            messages=[
                ChatMessage(role="user", content="Test", timestamp=datetime.utcnow()),
            ]
        )
        session = ChatSessionState(
            id=uuid4(),
            created_at=datetime.utcnow(),
            turns=[turn],
        )
        assert len(session.turns) == 1


class TestRunChatTurn:
    """Tests for run_chat_turn function."""

    @pytest.fixture
    def empty_session(self):
        """Create an empty chat session."""
        return ChatSessionState(
            id=uuid4(),
            created_at=datetime.utcnow(),
        )

    def test_help_intent_returns_help_message(self, empty_session):
        """Help intent should return help text."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        session, answer = run_chat_turn(empty_session, "help")
        assert "I can help you with" in answer
        assert len(session.turns) == 1

    def test_dataset_summary_without_id_prompts_for_id(self, empty_session):
        """Dataset summary without ID should prompt user."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        session, answer = run_chat_turn(empty_session, "summarize dataset")
        assert "Please specify a dataset ID" in answer

    def test_program_summary_without_id_prompts_for_id(self, empty_session):
        """Program summary without ID should prompt user."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        session, answer = run_chat_turn(empty_session, "program summary")
        assert "Please specify a program ID" in answer

    def test_signature_summary_without_id_prompts_for_id(self, empty_session):
        """Signature summary without ID should prompt user."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        session, answer = run_chat_turn(empty_session, "signature details")
        assert "Please specify a signature ID" in answer

    @patch("amprenta_rag.agent.chat_agent.query_rag")
    def test_freeform_rag_calls_query_rag(self, mock_query_rag, empty_session):
        """Freeform RAG intent should call query_rag."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        mock_result = MagicMock()
        mock_result.answer = "Test answer from RAG"
        mock_query_rag.return_value = mock_result

        session, answer = run_chat_turn(empty_session, "What genes are upregulated?")

        mock_query_rag.assert_called_once()
        assert answer == "Test answer from RAG"

    @patch("amprenta_rag.agent.chat_agent.query_rag")
    def test_freeform_rag_handles_empty_answer(self, mock_query_rag, empty_session):
        """Freeform RAG should handle empty/None answers gracefully."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        mock_result = MagicMock()
        mock_result.answer = None
        mock_query_rag.return_value = mock_result

        session, answer = run_chat_turn(empty_session, "Some query")

        assert "couldn't find enough evidence" in answer

    @patch("amprenta_rag.agent.chat_agent.cross_omics_dataset_summary_postgres")
    def test_dataset_summary_with_valid_id(self, mock_summary, empty_session):
        """Dataset summary with valid ID should call summary function."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        mock_summary.return_value = "Dataset summary text"

        session, answer = run_chat_turn(
            empty_session,
            "summarize dataset: 12345678-1234-1234-1234-123456789abc"
        )

        mock_summary.assert_called_once_with("12345678-1234-1234-1234-123456789abc")
        assert answer == "Dataset summary text"

    def test_session_accumulates_turns(self, empty_session):
        """Multiple chat turns should accumulate in session."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        session, _ = run_chat_turn(empty_session, "help")
        session, _ = run_chat_turn(session, "what can you do?")

        assert len(session.turns) == 2

    def test_error_handling(self, empty_session):
        """Errors inside try block should be caught and returned as error messages."""
        from amprenta_rag.agent.chat_agent import run_chat_turn

        # Mock query_rag to raise an exception (inside the try block)
        with patch("amprenta_rag.agent.chat_agent.query_rag") as mock_rag:
            mock_rag.side_effect = Exception("Test error")

            # This routes to freeform_rag which calls query_rag
            session, answer = run_chat_turn(empty_session, "random question")

            assert "Error:" in answer
            assert "Test error" in answer

