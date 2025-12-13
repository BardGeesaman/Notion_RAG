"""Q&A persistence models."""

from __future__ import annotations

import uuid
from datetime import datetime

from sqlalchemy import Boolean, Column, DateTime, Float, ForeignKey, Integer, JSON, String, Text
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import relationship

from amprenta_rag.database.base import Base


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class SavedQuestion(Base):
    """Saved user question for tracking Q&A."""

    __tablename__ = "saved_questions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    question_text = Column(Text, nullable=False)
    tags = Column(JSON, nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    is_archived = Column(Boolean, default=False, nullable=False)

    answers = relationship("SavedAnswer", back_populates="question", cascade="all, delete-orphan")
    created_by = relationship("User")


class SavedAnswer(Base):
    """Saved answer linked to a question with evidence and model metadata."""

    __tablename__ = "saved_answers"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    question_id = Column(UUID(as_uuid=True), ForeignKey("saved_questions.id"), nullable=False)
    answer_text = Column(Text, nullable=False)
    evidence = Column(JSON, nullable=True)
    model_used = Column(String(200), nullable=True)
    version = Column(Integer, default=1, nullable=False)
    confidence_score = Column(Float, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    question = relationship("SavedQuestion", back_populates="answers")


__all__ = ["SavedQuestion", "SavedAnswer"]

