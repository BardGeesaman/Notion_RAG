"""AI Expert Agent models for specialized domain assistance."""

from __future__ import annotations

from datetime import datetime
from typing import List, Optional, TYPE_CHECKING
from uuid import uuid4

from sqlalchemy import (
    Boolean, Column, DateTime, Float, ForeignKey, Integer, String, Table, Text, func
)
from sqlalchemy.dialects.postgresql import JSONB, UUID as PostgresUUID
from sqlalchemy.orm import Mapped, mapped_column, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import User


# Association table for expert-conversation many-to-many
expert_conversation_participants = Table(
    "expert_conversation_participants",
    Base.metadata,
    Column("expert_id", PostgresUUID(as_uuid=True), ForeignKey("expert_agents.id", ondelete="CASCADE"), primary_key=True),
    Column("conversation_id", PostgresUUID(as_uuid=True), ForeignKey("expert_conversations.id", ondelete="CASCADE"), primary_key=True),
    Column("joined_at", DateTime, default=func.now()),
)


class ExpertAgent(Base):
    """AI expert persona with specialized domain knowledge."""
    __tablename__ = "expert_agents"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    name: Mapped[str] = mapped_column(String(100), nullable=False)
    role: Mapped[str] = mapped_column(String(200), nullable=False)
    system_prompt: Mapped[str] = mapped_column(Text, nullable=False)
    prompt_version: Mapped[str] = mapped_column(String(20), nullable=False, default="1.0")  # P1 FIX
    specializations: Mapped[List[str]] = mapped_column(JSONB, nullable=False, default=list)
    is_active: Mapped[bool] = mapped_column(Boolean, default=True, nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    updated_at: Mapped[datetime] = mapped_column(DateTime, default=func.now(), onupdate=func.now())
    
    # Relationships
    conversations: Mapped[List["ExpertConversation"]] = relationship(
        secondary=expert_conversation_participants,
        back_populates="participants"
    )
    messages: Mapped[List["ExpertMessage"]] = relationship(back_populates="expert", cascade="all, delete-orphan")
    training_examples: Mapped[List["ExpertTrainingExample"]] = relationship(back_populates="expert", cascade="all, delete-orphan")
    knowledge_docs: Mapped[List["ExpertKnowledgeDoc"]] = relationship(back_populates="expert", cascade="all, delete-orphan")
    
    def __repr__(self) -> str:
        return f"<ExpertAgent(id={self.id}, name='{self.name}', role='{self.role}')>"


class ExpertConversation(Base):
    """Multi-expert conversation thread."""
    __tablename__ = "expert_conversations"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    user_id: Mapped[Optional[uuid4]] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"))
    title: Mapped[str] = mapped_column(String(500), nullable=False)
    context_entity_type: Mapped[Optional[str]] = mapped_column(String(100))
    context_entity_id: Mapped[Optional[uuid4]] = mapped_column(PostgresUUID(as_uuid=True))
    max_messages: Mapped[int] = mapped_column(Integer, default=50, nullable=False)  # P1 FIX
    is_panel: Mapped[bool] = mapped_column(Boolean, default=False, nullable=False)
    is_active: Mapped[bool] = mapped_column(Boolean, default=True, nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    last_message_at: Mapped[Optional[datetime]] = mapped_column(DateTime)
    
    # Relationships
    user: Mapped[Optional["User"]] = relationship(foreign_keys=[user_id])
    participants: Mapped[List["ExpertAgent"]] = relationship(
        secondary=expert_conversation_participants,
        back_populates="conversations"
    )
    messages: Mapped[List["ExpertMessage"]] = relationship(back_populates="conversation", cascade="all, delete-orphan")
    
    def __repr__(self) -> str:
        return f"<ExpertConversation(id={self.id}, title='{self.title}')>"


class ExpertMessage(Base):
    """Message in expert conversation."""
    __tablename__ = "expert_messages"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    conversation_id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("expert_conversations.id", ondelete="CASCADE"))
    expert_id: Mapped[Optional[uuid4]] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("expert_agents.id", ondelete="SET NULL"))
    role: Mapped[str] = mapped_column(String(20), nullable=False)  # user, assistant, system
    content: Mapped[str] = mapped_column(Text, nullable=False)
    prompt_version_used: Mapped[Optional[str]] = mapped_column(String(20))  # P1 FIX
    reasoning: Mapped[Optional[str]] = mapped_column(Text)
    citations: Mapped[Optional[List[str]]] = mapped_column(JSONB)
    token_count: Mapped[Optional[int]] = mapped_column(Integer)  # P1 FIX
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    
    # Relationships
    conversation: Mapped["ExpertConversation"] = relationship(back_populates="messages")
    expert: Mapped[Optional["ExpertAgent"]] = relationship(back_populates="messages")
    feedback: Mapped[List["ExpertFeedback"]] = relationship(back_populates="message", cascade="all, delete-orphan")
    
    def __repr__(self) -> str:
        return f"<ExpertMessage(id={self.id}, role='{self.role}')>"


class ExpertFeedback(Base):
    """User feedback on expert responses."""
    __tablename__ = "expert_feedback"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    message_id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("expert_messages.id", ondelete="CASCADE"))
    user_id: Mapped[Optional[uuid4]] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"))
    rating: Mapped[int] = mapped_column(Integer, nullable=False)  # P1 FIX: 1-5 integer rating
    correction: Mapped[Optional[str]] = mapped_column(Text)
    tags: Mapped[Optional[List[str]]] = mapped_column(JSONB)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    
    # Relationships
    message: Mapped["ExpertMessage"] = relationship(back_populates="feedback")
    user: Mapped[Optional["User"]] = relationship(foreign_keys=[user_id])
    
    def __repr__(self) -> str:
        return f"<ExpertFeedback(id={self.id}, rating={self.rating})>"


class ExpertTrainingExample(Base):
    """Training examples for expert fine-tuning."""
    __tablename__ = "expert_training_examples"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    expert_id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("expert_agents.id", ondelete="CASCADE"))
    question: Mapped[str] = mapped_column(Text, nullable=False)
    ideal_answer: Mapped[str] = mapped_column(Text, nullable=False)
    is_approved: Mapped[bool] = mapped_column(Boolean, default=False, nullable=False)  # P1 FIX
    prompt_version: Mapped[str] = mapped_column(String(20), nullable=False, default="1.0")  # P1 FIX
    source: Mapped[Optional[str]] = mapped_column(String(200))  # feedback, manual, generated
    created_by_id: Mapped[Optional[uuid4]] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"))
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    
    # Relationships
    expert: Mapped["ExpertAgent"] = relationship(back_populates="training_examples")
    created_by: Mapped[Optional["User"]] = relationship(foreign_keys=[created_by_id])
    
    def __repr__(self) -> str:
        return f"<ExpertTrainingExample(id={self.id}, expert='{self.expert.name if self.expert else 'None'}')>"


class ExpertKnowledgeDoc(Base):
    """RAG knowledge documents for expert context."""
    __tablename__ = "expert_knowledge_docs"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    expert_id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("expert_agents.id", ondelete="CASCADE"))
    namespace: Mapped[str] = mapped_column(String(100), nullable=False)  # P1 FIX
    title: Mapped[str] = mapped_column(String(500), nullable=False)
    content: Mapped[str] = mapped_column(Text, nullable=False)
    chunk_index: Mapped[int] = mapped_column(Integer, default=0, nullable=False)  # P1 FIX
    embedding: Mapped[Optional[List[float]]] = mapped_column(JSONB)
    embedding_model: Mapped[Optional[str]] = mapped_column(String(100))  # P1 FIX
    source_url: Mapped[Optional[str]] = mapped_column(String(1000))
    source_type: Mapped[Optional[str]] = mapped_column(String(50))  # paper, manual, web, internal
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    
    # Relationships
    expert: Mapped["ExpertAgent"] = relationship(back_populates="knowledge_docs")
    
    def __repr__(self) -> str:
        return f"<ExpertKnowledgeDoc(id={self.id}, title='{self.title[:50]}...')>"


__all__ = [
    "ExpertAgent",
    "ExpertConversation", 
    "ExpertMessage",
    "ExpertFeedback",
    "ExpertTrainingExample",
    "ExpertKnowledgeDoc",
    "expert_conversation_participants",
]
