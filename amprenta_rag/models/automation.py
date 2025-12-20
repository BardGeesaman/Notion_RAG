"""Workflow automation models."""

from __future__ import annotations

import uuid
from typing import List, Optional, TYPE_CHECKING

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, JSON, String, Text, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import User


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class WorkflowRule(Base):
    """Workflow automation rules."""

    __tablename__ = "workflow_rules"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    trigger_type = Column(String(100), nullable=False)  # e.g., "experiment_created", "compound_registered"
    trigger_config = Column(JSON, nullable=True)  # Filter conditions
    action_type = Column(String(100), nullable=False)  # e.g., "send_notification", "add_tag"
    action_config = Column(JSON, nullable=True)  # Action parameters
    is_active = Column(Boolean, default=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])
    executions: Mapped[List["WorkflowExecution"]] = relationship(
        "WorkflowExecution", back_populates="rule", cascade="all, delete-orphan"
    )


class WorkflowExecution(Base):
    """Workflow rule execution history."""

    __tablename__ = "workflow_executions"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    rule_id = Column(UUID(as_uuid=True), ForeignKey("workflow_rules.id"), nullable=False)
    trigger_context = Column(JSON, nullable=True)  # What triggered it
    status = Column(String(50), nullable=False, default="pending")  # pending, success, failed
    result = Column(JSON, nullable=True)  # Execution result/error
    triggered_at = Column(DateTime(timezone=True), server_default=func.now())
    completed_at = Column(DateTime(timezone=True), nullable=True)

    rule: Mapped["WorkflowRule"] = relationship("WorkflowRule", back_populates="executions")


__all__ = ["WorkflowRule", "WorkflowExecution"]

